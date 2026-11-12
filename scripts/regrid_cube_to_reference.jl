using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using Mmap
using Printf
using WCS
using VGT

function parse_args(args)
    opts = Dict{String,Any}(
        "source" => nothing,
        "reference" => nothing,
        "output" => nothing,
        "row-block" => "64",
    )
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("Unexpected positional argument: $arg"))
        key = arg[3:end]
        haskey(opts, key) || throw(ArgumentError("Unknown option --$key"))
        i == length(args) && throw(ArgumentError("Missing value for --$key"))
        opts[key] = args[i + 1]
        i += 2
    end
    for key in ("source", "reference", "output")
        opts[key] === nothing && throw(ArgumentError("--$key is required"))
    end
    opts["row-block"] = parse(Int, opts["row-block"])
    opts["row-block"] >= 1 || throw(ArgumentError("--row-block must be positive"))
    return opts
end

function header_value(header, key::AbstractString)
    haskey(header, key) || throw(ArgumentError("Missing FITS keyword $key"))
    return Float64(header[key])
end

function linear_axis(header, n::Int)
    crval = header_value(header, "CRVAL3")
    crpix = header_value(header, "CRPIX3")
    cdelt = header_value(header, "CDELT3")
    return [crval + (k - crpix) * cdelt for k in 1:n]
end

function wcs_transform(header)
    transforms = WCS.from_header(VGT._wcs_header_string(header); ignore_rejected=true)
    isempty(transforms) && throw(ArgumentError("Could not construct a WCS transform from FITS header"))
    return transforms[1]
end

function spatial_source_coordinates(target_wcs, source_wcs, source_header,
                                    x_values::AbstractVector, y_values::AbstractVector)
    n = length(x_values) * length(y_values)
    target_pixels = ones(Float64, 3, n)
    idx = 0
    for y in y_values, x in x_values
        idx += 1
        target_pixels[1, idx] = x
        target_pixels[2, idx] = y
    end
    sky = WCS.pix_to_world(target_wcs, target_pixels)
    # Spatial WCS is separable here, but supplying a valid source velocity
    # avoids asking the source transform to invert an unrelated target velocity.
    sky[3, :] .= header_value(source_header, "CRVAL3")
    source_pixels = WCS.world_to_pix(source_wcs, sky)
    nx = length(x_values)
    ny = length(y_values)
    return reshape(source_pixels[1, :], nx, ny), reshape(source_pixels[2, :], nx, ny)
end

function read_channel(f, channel::Int)
    return Float32.(read(f[1], :, :, channel:channel)[:, :, 1])
end

function bilinear_stencil(source_x::AbstractMatrix, source_y::AbstractMatrix,
                          source_nx::Int, source_ny::Int)
    nx, ny = size(source_x)
    xlo = zeros(Int32, nx, ny)
    xhi = zeros(Int32, nx, ny)
    ylo = zeros(Int32, nx, ny)
    yhi = zeros(Int32, nx, ny)
    wx = zeros(Float32, nx, ny)
    wy = zeros(Float32, nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        bx = VGT.linear_bracket(source_x[i, j], source_nx)
        by = VGT.linear_bracket(source_y[i, j], source_ny)
        if bx !== nothing && by !== nothing
            xlo[i, j] = bx.lo
            xhi[i, j] = bx.hi
            ylo[i, j] = by.lo
            yhi[i, j] = by.hi
            wx[i, j] = bx.w
            wy[i, j] = by.w
        end
    end
    return (; xlo, xhi, ylo, yhi, wx, wy)
end

function regrid_channel!(destination::AbstractMatrix, source_low::AbstractMatrix,
                         source_high::AbstractMatrix, spectral_weight::Float64,
                         stencil)
    nx, ny = size(destination)
    @inbounds Threads.@threads for local_y in 1:ny
        for local_x in 1:nx
            xlo = stencil.xlo[local_x, local_y]
            ylo = stencil.ylo[local_x, local_y]
            if xlo == 0 || ylo == 0
                destination[local_x, local_y] = NaN32
                continue
            end
            xhi = stencil.xhi[local_x, local_y]
            yhi = stencil.yhi[local_x, local_y]
            tx = stencil.wx[local_x, local_y]
            ty = stencil.wy[local_x, local_y]

            low11 = source_low[xlo, ylo]
            low21 = source_low[xhi, ylo]
            low12 = source_low[xlo, yhi]
            low22 = source_low[xhi, yhi]
            high11 = source_high[xlo, ylo]
            high21 = source_high[xhi, ylo]
            high12 = source_high[xlo, yhi]
            high22 = source_high[xhi, yhi]
            if isfinite(low11) && isfinite(low21) && isfinite(low12) && isfinite(low22) &&
               isfinite(high11) && isfinite(high21) && isfinite(high12) && isfinite(high22)
                low = (1 - tx) * (1 - ty) * low11 + tx * (1 - ty) * low21 +
                      (1 - tx) * ty * low12 + tx * ty * low22
                high = (1 - tx) * (1 - ty) * high11 + tx * (1 - ty) * high21 +
                       (1 - tx) * ty * high12 + tx * ty * high22
                destination[local_x, local_y] = Float32((1 - spectral_weight) * low + spectral_weight * high)
            else
                destination[local_x, local_y] = NaN32
            end
        end
    end
    return destination
end

function make_output_header(reference_header, first_channel::Int, last_channel::Int)
    output_header = deepcopy(reference_header)
    old_crpix3 = header_value(reference_header, "CRPIX3")
    output_header["NAXIS3"] = last_channel - first_channel + 1
    output_header["CRPIX3"] = old_crpix3 - first_channel + 1
    output_header["BUNIT"] = "K"
    output_header["REGRIDX"] = "BILINEAR"
    output_header["REGRIDV"] = "LINEAR"
    return output_header
end

function regrid_cube(opts)
    source_path = opts["source"]
    reference_path = opts["reference"]
    output_path = opts["output"]
    isfile(source_path) || throw(ArgumentError("Source cube not found: $source_path"))
    isfile(reference_path) || throw(ArgumentError("Reference cube not found: $reference_path"))
    output_path = abspath(output_path)
    mkpath(dirname(output_path))

    FITS(source_path) do source_fits
        FITS(reference_path) do reference_fits
            source_hdu = source_fits[1]
            reference_hdu = reference_fits[1]
            source_dims = size(source_hdu)
            reference_dims = size(reference_hdu)
            length(source_dims) == 3 || throw(ArgumentError("Source must be a 3-D FITS image"))
            length(reference_dims) == 3 || throw(ArgumentError("Reference must be a 3-D FITS image"))

            source_header = read_header(source_hdu)
            reference_header = read_header(reference_hdu)
            source_velocity = linear_axis(source_header, source_dims[3])
            reference_velocity = linear_axis(reference_header, reference_dims[3])
            overlap = VGT.velocity_overlap_indices(
                reference_header, reference_dims[3], source_header, source_dims[3])
            first_channel, last_channel = first(overlap), last(overlap)
            target_velocity = reference_velocity[overlap]
            source_velocity_min, source_velocity_max = extrema(source_velocity)

            @printf("Source: %s  size=%s  velocity=%.3f..%.3f\n",
                    source_path, source_dims, source_velocity_min, source_velocity_max)
            @printf("Reference: %s  size=%s  velocity=%.3f..%.3f\n",
                    reference_path, reference_dims, minimum(reference_velocity), maximum(reference_velocity))
            @printf("Output velocity overlap: channels %d:%d, %.3f..%.3f\n",
                    first_channel, last_channel, first(target_velocity), last(target_velocity))

            source_wcs = wcs_transform(source_header)
            reference_wcs = wcs_transform(reference_header)
            nx, ny = reference_dims[1], reference_dims[2]
            output_nv = length(overlap)
            output_dims = (nx, ny, output_nv)
            raw_path = string(output_path, ".raw.partial")
            partial_fits_path = string(output_path, ".partial")
            isfile(output_path) && throw(ArgumentError("Output already exists: $output_path"))
            isfile(partial_fits_path) && throw(ArgumentError("Partial FITS already exists: $partial_fits_path"))
            isfile(raw_path) && throw(ArgumentError("Partial raw buffer already exists: $raw_path"))

            output_io = open(raw_path, "w+")
            try
                truncate(output_io, prod(output_dims) * sizeof(Float32))
            finally
                close(output_io)
            end
            output_io = open(raw_path, "r+")
            output = Mmap.mmap(output_io, Array{Float32,3}, output_dims, 0; shared=true)
            fill!(output, NaN32)

            row_block = opts["row-block"]
            try
                for y_start in 1:row_block:ny
                    y_end = min(ny, y_start + row_block - 1)
                    x_values = collect(1:nx)
                    y_values = collect(y_start:y_end)
                    source_x, source_y = spatial_source_coordinates(
                        reference_wcs, source_wcs, source_header, x_values, y_values)
                    stencil = bilinear_stencil(source_x, source_y, source_dims[1], source_dims[2])
                    destination = Matrix{Float32}(undef, nx, length(y_values))

                    for (output_channel, target_channel) in enumerate(overlap)
                        target_v = reference_velocity[target_channel]
                        source_position = VGT.source_axis_position(
                            target_v,
                            header_value(source_header, "CRVAL3"),
                            header_value(source_header, "CRPIX3"),
                            header_value(source_header, "CDELT3"),
                        )
                        bracket = VGT.linear_bracket(source_position, source_dims[3])
                        bracket === nothing && continue
                        source_low = read_channel(source_fits, bracket.lo)
                        source_high = bracket.hi == bracket.lo ? source_low : read_channel(source_fits, bracket.hi)
                        regrid_channel!(destination, source_low, source_high, bracket.w,
                                        stencil)
                        @views output[:, y_start:y_end, output_channel] .= destination
                    end
                    @printf("Processed rows %d:%d/%d\n", y_start, y_end, ny)
                end
                Mmap.sync!(output)
                output_header = make_output_header(reference_header, first_channel, last_channel)
                output_header["SRCFILE"] = basename(source_path)
                output_header["REFFILE"] = basename(reference_path)
                FITS(partial_fits_path, "w") do output_fits
                    write(output_fits, output; header=output_header)
                end
            finally
                Mmap.sync!(output)
                close(output_io)
            end

            FITS(partial_fits_path) do check_fits
                check_dims = size(check_fits[1])
                check_dims == output_dims || throw(ArgumentError("Output size $check_dims does not match $output_dims"))
                check_header = read_header(check_fits[1])
                @printf("Validated output: size=%s CRPIX3=%s CRVAL3=%s CDELT3=%s\n",
                        check_dims, check_header["CRPIX3"], check_header["CRVAL3"], check_header["CDELT3"])
            end
            mv(partial_fits_path, output_path; force=true)
            rm(raw_path; force=true)
            println("Wrote ", output_path)
        end
    end
    return output_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        regrid_cube(parse_args(ARGS))
    catch err
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end

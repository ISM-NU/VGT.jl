using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using Printf
using TOML
using VGT

function default_options()
    return Dict(
        "input" => nothing,
        "output-dir" => rawdatadir("velocity_slabs"),
        "ranges" => "-50:-20,-20:0,0:30,30:50,50:90",
        "prefix" => "",
        "force" => "false",
    )
end

function parse_args(args)
    opts = default_options()
    i = 1
    while i <= length(args)
        arg = args[i]
        startswith(arg, "--") || throw(ArgumentError("Unexpected positional argument: $arg"))
        key = arg[3:end]
        haskey(opts, key) || throw(ArgumentError("Unknown option --$key"))
        if key == "force"
            opts[key] = "true"
            i += 1
        else
            i == length(args) && throw(ArgumentError("Missing value for --$key"))
            opts[key] = args[i + 1]
            i += 2
        end
    end
    opts["input"] === nothing && throw(ArgumentError("--input is required"))
    return opts
end

function parse_bool(raw::AbstractString)
    raw in ("true", "yes", "1") && return true
    raw in ("false", "no", "0") && return false
    throw(ArgumentError("Expected boolean string, got $raw"))
end

function parse_ranges(raw::AbstractString)
    ranges = Tuple{Float64, Float64}[]
    for part in split(raw, ',')
        text = strip(part)
        isempty(text) && continue
        pieces = split(text, ':'; limit=2)
        length(pieces) == 2 || throw(ArgumentError("Velocity range must be lo:hi, got $text"))
        lo = parse(Float64, strip(pieces[1]))
        hi = parse(Float64, strip(pieces[2]))
        lo <= hi || throw(ArgumentError("Velocity range must have lo <= hi, got $text"))
        push!(ranges, (lo, hi))
    end
    isempty(ranges) && throw(ArgumentError("--ranges must contain at least one lo:hi range"))
    return ranges
end

function velocity_axis_kms(header, nv::Int)
    crval = Float64(header["CRVAL3"])
    crpix = Float64(header["CRPIX3"])
    cdelt = Float64(header["CDELT3"])
    axis_native = [crval + (k - crpix) * cdelt for k in 1:nv]
    scale = maximum(abs.(axis_native)) > 1_000 ? 1e-3 : 1.0
    return axis_native .* scale, scale
end

function safe_label(lo::Real, hi::Real)
    fmt(v) = replace(@sprintf("%.1f", v), "-" => "m", "." => "p")
    return string("v", fmt(lo), "_to_", fmt(hi), "kms")
end

function update_slab_header!(header, first_channel::Int)
    try
        old_crpix3 = Float64(header["CRPIX3"])
        old_crval3 = Float64(header["CRVAL3"])
        cdelt3 = Float64(header["CDELT3"])
        header["CRPIX3"] = 1.0
        header["CRVAL3"] = old_crval3 + (first_channel - old_crpix3) * cdelt3
    catch err
        @warn "Could not update slab velocity WCS keywords; writing data with original header" exception=(err, catch_backtrace())
    end
    return header
end

function write_slab_metadata(path::AbstractString, metadata)
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
end

function make_velocity_slabs(opts)
    ranges = parse_ranges(opts["ranges"])
    force = parse_bool(opts["force"])
    input = opts["input"]
    output_dir = opts["output-dir"]
    mkpath(output_dir)
    f = FITS(input)
    try
        nx, ny, nv = size(f[1])
        header = read_header(f[1])
        velocities_kms, native_to_kms = velocity_axis_kms(header, nv)
        stem = splitext(basename(input))[1]
        written = String[]
        for (lo, hi) in ranges
            indices = findall(v -> lo <= v <= hi, velocities_kms)
            if isempty(indices)
                @warn "No channels in requested velocity range" lo hi
                continue
            end
            label = safe_label(lo, hi)
            out_stem = isempty(opts["prefix"]) ? string(stem, "_", label) :
                       string(opts["prefix"], "_", label)
            fits_path = joinpath(output_dir, string(out_stem, ".fits"))
            if isfile(fits_path) && !force
                @warn "Output exists; use --force to overwrite" path=fits_path
                push!(written, fits_path)
                continue
            end
            slab = read(f[1], :, :, minimum(indices):maximum(indices))
            slab_header = deepcopy(header)
            update_slab_header!(slab_header, minimum(indices))
            VGT.write_fits(fits_path, slab; header=slab_header)
            metadata_path = joinpath(output_dir, string(out_stem, ".toml"))
            write_slab_metadata(metadata_path, Dict(
                "input" => input,
                "input_absolute_path" => abspath(input),
                "output" => fits_path,
                "velocity_range_kms" => Dict("lo" => lo, "hi" => hi),
                "channel_range" => Dict("first" => minimum(indices), "last" => maximum(indices)),
                "n_channels" => length(indices),
                "native_velocity_to_kms_scale" => native_to_kms,
                "velocity_min_kms" => minimum(velocities_kms[indices]),
                "velocity_max_kms" => maximum(velocities_kms[indices]),
                "cube_size" => Dict("nx" => nx, "ny" => ny, "nv" => nv),
            ))
            push!(written, fits_path)
            @printf "Wrote %s (%d channels, %.2f..%.2f km/s)\n" fits_path length(indices) minimum(velocities_kms[indices]) maximum(velocities_kms[indices])
        end
        return written
    finally
        close(f)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        written = make_velocity_slabs(parse_args(ARGS))
        println("Wrote ", length(written), " slab FITS files.")
    catch err
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end

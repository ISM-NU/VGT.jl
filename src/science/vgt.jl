# Core Velocity Gradient Technique (VGT) computation.
#
# The input cube is interpreted as intensity I(x, y, v): two sky axes and one
# velocity/channel axis. For each channel, the code computes local spatial
# gradients, estimates the dominant gradient orientation in a sliding window,
# and accumulates those orientations into doubled-angle Stokes-like Q/U maps.

function read_cube_metadata(path::AbstractString)
    # Read only dimensions and header metadata. The full cube can be many GB, so
    # the pipeline avoids loading it all into memory at this stage.
    f = FITS(path)
    try
        dims = size(f[1])
        header = read_header(f[1])
        return dims[1], dims[2], dims[3], header
    finally
        close(f)
    end
end

function _fits_card_value(value)
    if value isa AbstractString
        escaped = replace(value, "'" => "''")
        return "'$(rpad(escaped, 8))'"
    elseif value isa Bool
        return value ? "T" : "F"
    elseif value isa Integer
        return string(value)
    elseif value isa Real
        return string(Float64(value))
    end
    return string(value)
end

function _fits_header_card(key::AbstractString, value)
    length(key) <= 8 || throw(ArgumentError("FITS keyword is longer than 8 characters: $key"))
    body = rpad(key, 8) * "= " * lpad(_fits_card_value(value), 20)
    return rpad(body, 80)
end

function _wcs_header_string(header)
    cards = String[]
    for key in keys(header)
        push!(cards, _fits_header_card(String(key), header[key]))
    end
    push!(cards, rpad("END", 80))
    return join(cards)
end

function _linear_axis_coordinates(header, axis::Int, n::Int)
    crpix = Float64(header["CRPIX$(axis)"])
    crval = Float64(header["CRVAL$(axis)"])
    cdelt = Float64(header["CDELT$(axis)"])
    return [crval + (idx - crpix) * cdelt for idx in 1:n]
end

function _linear_getcoor(header, nx, ny, nv)
    RA = zeros(nx)
    DEC = zeros(ny)
    vz = zeros(nv)
    RA .= _linear_axis_coordinates(header, 1, nx)
    DEC .= _linear_axis_coordinates(header, 2, ny)
    vz .= _linear_axis_coordinates(header, 3, nv)
    xi = Float64(header["CRVAL1"]) + (0 - Float64(header["CRPIX1"])) * Float64(header["CDELT1"])
    xe = Float64(header["CRVAL1"]) + (nx - Float64(header["CRPIX1"])) * Float64(header["CDELT1"])
    yi = Float64(header["CRVAL2"]) + (0 - Float64(header["CRPIX2"])) * Float64(header["CDELT2"])
    ye = Float64(header["CRVAL2"]) + (ny - Float64(header["CRPIX2"])) * Float64(header["CDELT2"])
    @debug "map extent" xi xe yi ye
    return RA, DEC, vz
end

function _wcs_axis_coordinates(transform, axis::Int, n::Int, naxis::Int)
    pixcoords = ones(Float64, naxis, n)
    @inbounds for idx in 1:n
        pixcoords[axis, idx] = idx
    end
    world = WCS.pix_to_world(transform, pixcoords)
    return vec(world[axis, :])
end

function getcoor(header, nx, ny, nv)
    # FITS WCS metadata converts pixel indices into sky/velocity coordinates.
    # The pipeline can still run without this metadata; it is used for reporting
    # and sanity checks, not for the VGT numerical result.
    try
        transforms = WCS.from_header(_wcs_header_string(header); ignore_rejected=true)
        if isempty(transforms)
            return _linear_getcoor(header, nx, ny, nv)
        end

        transform = transforms[1]
        naxis = max(3, transform.naxis)
        RA = _wcs_axis_coordinates(transform, 1, nx, naxis)
        DEC = _wcs_axis_coordinates(transform, 2, ny, naxis)
        vz = _wcs_axis_coordinates(transform, 3, nv, naxis)

        corners = ones(Float64, naxis, 4)
        corners[1:2, :] .= [0.0 nx 0.0 nx;
                            0.0 0.0 ny ny]
        world_corners = WCS.pix_to_world(transform, corners)
        xi, xe = extrema(world_corners[1, :])
        yi, ye = extrema(world_corners[2, :])
        @debug "map extent" xi xe yi ye
        return RA, DEC, vz
    catch err
        is_interrupt_exception(err) && rethrow()
        @warn "WCS coordinate extraction failed; falling back to linear CRPIX/CRVAL/CDELT coordinates" exception=(err, catch_backtrace())
        return _linear_getcoor(header, nx, ny, nv)
    end
end

@inline function _accumulate_weighted_diagnostic!(num::AbstractMatrix{Float64},
                                                 weight_sum::AbstractMatrix{Float64},
                                                 i::Int, j::Int, value::Real, weight::Real)
    if isfinite(value) && weight > 0.0
        num[i, j] = muladd(Float64(weight), Float64(value), num[i, j])
        weight_sum[i, j] += Float64(weight)
    end
    return nothing
end

@inline function _weighted_diagnostic_value(num::AbstractArray{Float64},
                                            weight_sum::AbstractArray{Float64}, idx)
    return weight_sum[idx] > 0.0 ? num[idx] / weight_sum[idx] : NaN
end


# Gradient stencils at interior pixel (i, j) of a 2D image; the sign convention
# (left minus right) is the historical one and does not matter for axial angles.
# :central_checked  (a[i-1] - a[i+1]) / 2.
# :sobel_checked    central difference times [1, 2, 1] / 4 across it (Sobel / 8).
# :isotropic_checked central difference times [1, 4, 1] / 6 across it.
# On a plane wave of wavenumber k (rad/pixel) at gradient angle theta, central
# differences return the angle theta + (k^2 / 24) sin(4 theta) + O(k^4): the
# square grid pulls gradients towards the diagonals. Smoothing across the
# difference with [a, 1 - 2a, a] adds -(a k^2 / 4) sin(4 theta); Sobel (a = 1/4)
# overshoots to -(k^2 / 48) sin(4 theta), and a = 1/6 cancels the k^2 term
# (Farid & Simoncelli 2004 treat the general rotation-invariance problem).
# test/angle_calibration.jl checks these responses; validation in
# docs/validation-angle-calibration.md.
@inline function _stencil_gradient(img::AbstractMatrix, i::Int, j::Int, mode::Symbol)
    @inbounds begin
        if mode === :central_checked
            return (img[i - 1, j] - img[i + 1, j]) / 2, (img[i, j - 1] - img[i, j + 1]) / 2
        elseif mode === :sobel_checked
            gx = (
                img[i - 1, j - 1] + 2 * img[i - 1, j] + img[i - 1, j + 1] -
                img[i + 1, j - 1] - 2 * img[i + 1, j] - img[i + 1, j + 1]
            ) / 8
            gy = (
                img[i - 1, j - 1] + 2 * img[i, j - 1] + img[i + 1, j - 1] -
                img[i - 1, j + 1] - 2 * img[i, j + 1] - img[i + 1, j + 1]
            ) / 8
            return gx, gy
        elseif mode === :isotropic_checked
            gx = (
                img[i - 1, j - 1] + 4 * img[i - 1, j] + img[i - 1, j + 1] -
                img[i + 1, j - 1] - 4 * img[i + 1, j] - img[i + 1, j + 1]
            ) / 12
            gy = (
                img[i - 1, j - 1] + 4 * img[i, j - 1] + img[i + 1, j - 1] -
                img[i - 1, j + 1] - 4 * img[i, j + 1] - img[i + 1, j + 1]
            ) / 12
            return gx, gy
        end
    end
    throw(ArgumentError("No stencil for gradient_mode=$mode"))
end

# Input pixels read by the gradient at (i, j), as (di, dj) offsets. These are
# exactly the pixels _stencil_gradient (or the backward difference) touches;
# the corners of the 3x3 block are not read by the central difference.
const _CROSS_STENCIL = ((0, 0), (-1, 0), (1, 0), (0, -1), (0, 1))
const _BLOCK_STENCIL = Tuple((di, dj) for dj in -1:1 for di in -1:1)
const _BACKWARD_STENCIL = ((0, 0), (-1, 0), (0, -1))

function _gradient_stencil_offsets(gradient_mode::Symbol)
    gradient_mode === :central_checked && return _CROSS_STENCIL
    gradient_mode === :nan_edge_backward && return _BACKWARD_STENCIL
    return _BLOCK_STENCIL
end

"""
    _edge_footprint_mask(cube, weight_cube, dn, first_channel, last_channel, gradient_mode)

`BitMatrix`, true where the footprint of the VGT window at that pixel reaches a
non-finite input pixel or the array border. The footprint is what the kernel
reads: the window of pixels `i-dn/2+1 : i+dn/2` (and the same in j) plus the
gradient stencil of every pixel in it, so it spans `i-dn/2 : i+dn/2+1` for the
central, Sobel and isotropic stencils and `i-dn/2 : i+dn/2` for the backward
difference. A pixel is non-finite when it is NaN or Inf in any channel of
`first_channel:last_channel` of `cube` or `weight_cube`. A fully blank channel
therefore masks every pixel; trim such channels (`velocity_range_kms`) first. The
rule depends only on the data, so chunked and one-piece runs agree. Pixels below
the noise threshold are not non-finite and do not enter.
"""
function _edge_footprint_mask(cube::AbstractArray{<:Any, 3}, weight_cube, dn::Int,
                              first_channel::Int, last_channel::Int, gradient_mode::Symbol)
    nx, ny, _ = size(cube)
    bad = falses(nx, ny)
    for src in (cube, weight_cube)
        src === nothing && continue
        for channel in first_channel:last_channel
            img = @view src[:, :, channel]
            @inbounds for idx in eachindex(img)
                isfinite(img[idx]) || (bad[idx] = true)
            end
        end
    end
    offsets = _gradient_stencil_offsets(gradient_mode)
    # grad_bad: the gradient at (i, j) is NaN (stencil leaves the array or hits a bad pixel).
    grad_bad = falses(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        for (di, dj) in offsets
            ii, jj = i + di, j + dj
            if ii < 1 || ii > nx || jj < 1 || jj > ny || bad[ii, jj]
                grad_bad[i, j] = true
                break
            end
        end
    end
    # Summed-area table of grad_bad, so each window is counted in O(1).
    sat = zeros(Int, nx + 1, ny + 1)
    @inbounds for j in 1:ny, i in 1:nx
        sat[i + 1, j + 1] = grad_bad[i, j] + sat[i, j + 1] + sat[i + 1, j] - sat[i, j]
    end
    half = div(dn, 2)
    masked = trues(nx, ny)
    @inbounds for j in 1:ny, i in 1:nx
        is, ie, js, je = i - half + 1, i + half, j - half + 1, j + half
        (is < 1 || ie > nx || js < 1 || je > ny) && continue
        count = sat[ie + 1, je + 1] - sat[is, je + 1] - sat[ie + 1, js] + sat[is, js]
        count == 0 && (masked[i, j] = false)
    end
    return masked
end

function _vgt_qu_mv_exact_impl(cube::AbstractArray{<:Any, 3}, dn, first_channel, last_channel, noise;
                               min_valid::Int=100,
                               fitter_mode::Symbol=:axial_vonmises,
                               collect_profile::Bool=false,
                               rbar2_weight::Bool=false,
                               rbar2_gate::Float64=0.0,
                               gradient_mode::Symbol=:central_checked,
                               gradient_whiten::Bool=false,
                               gradient_noise_cov::AbstractVector{<:Real}=Float64[],
                               gradient_noise_cov_map=nothing,
                               weight_cube=nothing,
                               edge_mask::Bool=false)
    # This is the science kernel. It returns Q/U maps accumulated over channels.
    # The final magnetic-field angle is computed later from these maps.
    nx, ny, nv = size(cube)
    iseven(dn) || throw(ArgumentError("dn must be even; got $dn"))
    min_valid >= 1 || throw(ArgumentError("min_valid must be >= 1"))
    min_valid <= dn * dn ||
        throw(ArgumentError("min_valid=$min_valid exceeds the maximum window area dn^2=$(dn * dn)"))
    dn <= nx || throw(ArgumentError("dn=$dn is too large for nx=$nx (need nx ≥ $dn)"))
    dn <= ny || throw(ArgumentError("dn=$dn is too large for ny=$ny (need ny ≥ $dn)"))
    gradient_mode in GRADIENT_MODES ||
        throw(ArgumentError("Unsupported gradient_mode=$gradient_mode"))
    0.0 <= rbar2_gate <= 1.0 ||
        throw(ArgumentError("rbar2_gate must be in [0, 1], got $rbar2_gate"))
    _validate_gradient_whiten(gradient_whiten, gradient_noise_cov,
                              gradient_noise_cov_map === nothing ? "" : "<map>", gradient_mode)
    gradient_noise_cov_map === nothing || size(gradient_noise_cov_map) == (nx, ny, 3) ||
        throw(ArgumentError("gradient_noise_cov_map size $(size(gradient_noise_cov_map)) must be ($nx, $ny, 3)"))
    # Whitening maps each gradient by W = C^(-1/2) (C scaled to unit determinant),
    # so noise gradients with covariance C become isotropic; the fitted window
    # angle is mapped back with S = C^(1/2). Skipped when off, so the default
    # output is unchanged.
    whiten_maps = (gradient_whiten && gradient_noise_cov_map !== nothing) ?
                  _whitening_maps(gradient_noise_cov_map) : nothing
    whiten_w = (gradient_whiten && gradient_noise_cov_map === nothing) ?
               _gradient_whitening(gradient_noise_cov) : nothing

    # q_accum and u_accum store orientation information. We use cos(2*angle) and
    # sin(2*angle) because field orientations are axial: angle and angle + pi
    # describe the same headless direction.
    q_accum = zeros(nx, ny)
    u_accum = zeros(nx, ny)
    # Opt-in: pixels whose window footprint reaches a non-finite pixel or the
    # array border get NaN Q/U ("no direction") instead of a fit on a partial
    # window. Off by default, so the default output is unchanged.
    footprint_masked = edge_mask ?
        _edge_footprint_mask(cube, weight_cube, dn, first_channel, last_channel, gradient_mode) : nothing
    # w_accum sums the same Stokes weights that enter Q/U, over fitted windows
    # only. It is the reliability denominator: with weights >= 0,
    # sqrt(Q^2 + U^2) <= W by the triangle inequality, so P = |Q, U| / W <= 1.
    w_accum = zeros(nx, ny)
    grad_x = zeros(nx, ny)
    grad_y = zeros(nx, ny)
    half_dn = div(dn, 2)
    # weight_cube (same size as cube) supplies the noise mask and the Stokes
    # weights while the gradients come from cube, e.g. gradients of the VDA
    # velocity-caustics cube (which can be negative) weighted by the observed
    # intensity. Without it both come from cube.
    if weight_cube !== nothing
        size(weight_cube) == size(cube) ||
            throw(ArgumentError("weight_cube size $(size(weight_cube)) must match cube size $(size(cube))"))
    end
    scalar_noise = noise isa Number
    if !scalar_noise
        size(noise) == size(cube) ||
            throw(ArgumentError("noise array size $(size(noise)) must match cube size $(size(cube))"))
    end

    profile = Dict{Symbol, Any}(
        :collect_profile => collect_profile,
        :gradient_ns => Int[],
        :window_ns => Int[],
        :channels => Int[],
        :fit_attempts => 0,
        :fit_successes => 0,
        :fit_failures => 0,
        :fit_ns => 0,
        :valid_windows => 0,
        :fit_low_confidence_windows => 0,
        :fit_rbar2 => nothing,
        :fit_attempt_fraction => nothing,
        :fit_support => nothing,
        :fitted_weight => w_accum,
    )

    fit_rbar2_num = collect_profile ? zeros(nx, ny) : Matrix{Float64}(undef, 0, 0)
    fit_rbar2_w = collect_profile ? zeros(nx, ny) : Matrix{Float64}(undef, 0, 0)
    fit_attempt_count = collect_profile ? zeros(nx, ny) : Matrix{Float64}(undef, 0, 0)
    fit_success_count = collect_profile ? zeros(nx, ny) : Matrix{Float64}(undef, 0, 0)

    for channel in first_channel:last_channel
        channel_image = @view cube[:, :, channel]
        mask_image = weight_cube === nothing ? channel_image : @view weight_cube[:, :, channel]

        # Compute backward spatial differences within one velocity channel.
        # The production mode marks the first row/column invalid instead of
        # wrapping around the image edge, which avoids inventing a gradient
        # between unrelated sky pixels.
        t_grad = collect_profile ? time_ns() : Int64(0)
        @inbounds for j in 1:ny, i in 1:nx
            if gradient_mode === :nan_edge_backward
                grad_x[i, j] = i == 1 ? NaN : channel_image[i - 1, j] - channel_image[i, j]
                grad_y[i, j] = j == 1 ? NaN : channel_image[i, j - 1] - channel_image[i, j]
                # The backward compatibility path only masks the center pixel,
                # preserving the historical production semantics.
                threshold = scalar_noise ? noise : noise[i, j, channel]
                if mask_image[i, j] < threshold
                    grad_x[i, j] = NaN
                    grad_y[i, j] = NaN
                end
            elseif gradient_mode === :central_checked
                if i == 1 || i == nx || j == 1 || j == ny
                    grad_x[i, j] = NaN
                    grad_y[i, j] = NaN
                    continue
                end
                center_ok = mask_image[i, j] >= (scalar_noise ? noise : noise[i, j, channel])
                left_ok = mask_image[i - 1, j] >= (scalar_noise ? noise : noise[i - 1, j, channel])
                right_ok = mask_image[i + 1, j] >= (scalar_noise ? noise : noise[i + 1, j, channel])
                down_ok = mask_image[i, j - 1] >= (scalar_noise ? noise : noise[i, j - 1, channel])
                up_ok = mask_image[i, j + 1] >= (scalar_noise ? noise : noise[i, j + 1, channel])
                if center_ok && left_ok && right_ok && down_ok && up_ok
                    gx, gy = _stencil_gradient(channel_image, i, j, :central_checked)
                    if whiten_maps !== nothing
                        gx, gy = _whiten_gradient_at(whiten_maps, i, j, gx, gy)
                    elseif whiten_w !== nothing
                        gx, gy = _whiten_gradient(whiten_w, gx, gy)
                    end
                    grad_x[i, j] = gx
                    grad_y[i, j] = gy
                else
                    grad_x[i, j] = NaN
                    grad_y[i, j] = NaN
                end
            else # :sobel_checked or :isotropic_checked (3x3 stencils)
                if i == 1 || i == nx || j == 1 || j == ny
                    grad_x[i, j] = NaN
                    grad_y[i, j] = NaN
                    continue
                end
                valid_stencil = true
                for jj in (j - 1):(j + 1), ii in (i - 1):(i + 1)
                    threshold = scalar_noise ? noise : noise[ii, jj, channel]
                    if mask_image[ii, jj] < threshold
                        valid_stencil = false
                        break
                    end
                end
                if valid_stencil
                    gx, gy = _stencil_gradient(channel_image, i, j, gradient_mode)
                    if whiten_maps !== nothing
                        gx, gy = _whiten_gradient_at(whiten_maps, i, j, gx, gy)
                    elseif whiten_w !== nothing
                        gx, gy = _whiten_gradient(whiten_w, gx, gy)
                    end
                    grad_x[i, j] = gx
                    grad_y[i, j] = gy
                else
                    grad_x[i, j] = NaN
                    grad_y[i, j] = NaN
                end
            end
        end
        if collect_profile
            push!(profile[:channels], channel)
            push!(profile[:gradient_ns], time_ns() - t_grad)
        end

        fit_attempts_tls = collect_profile ? zeros(Int, Threads.maxthreadid()) : Int[]
        fit_successes_tls = collect_profile ? zeros(Int, Threads.maxthreadid()) : Int[]
        fit_failures_tls = collect_profile ? zeros(Int, Threads.maxthreadid()) : Int[]
        fit_ns_tls = collect_profile ? zeros(Int, Threads.maxthreadid()) : Int[]
        valid_windows_tls = collect_profile ? zeros(Int, Threads.maxthreadid()) : Int[]
        fit_low_confidence_tls = collect_profile ? zeros(Int, Threads.maxthreadid()) : Int[]

        t_fit_stage = collect_profile ? time_ns() : Int64(0)
        Threads.@threads for idx in 1:(nx * ny)
            ii = mod1(idx, nx)
            jj = cld(idx, nx)
            if ii < half_dn || ii > nx - half_dn || jj < half_dn || jj > ny - half_dn
                continue
            end
            footprint_masked !== nothing && footprint_masked[ii, jj] && continue

            tid = collect_profile ? Threads.threadid() : 1
            is = ii - half_dn + 1
            ie = ii + half_dn
            js = jj - half_dn + 1
            je = jj + half_dn
            grad_x_window = @view grad_x[is:ie, js:je]
            grad_y_window = @view grad_y[is:ie, js:je]

            # A local window must contain enough valid gradient samples before
            # fitting an orientation. This prevents isolated noisy pixels from
            # producing a misleading field estimate.
            valid = 0
            for v in grad_x_window
                if !isnan(v)
                    valid += 1
                end
            end
            valid >= min_valid || continue

            if collect_profile
                valid_windows_tls[tid] += 1
                fit_attempts_tls[tid] += 1
                fit_attempt_count[ii, jj] += 1
            end

            try
                t_fit = collect_profile ? time_ns() : Int64(0)
                result = fit_window_result(fitter_mode, grad_x_window, grad_y_window)
                if collect_profile
                    fit_ns_tls[tid] += time_ns() - t_fit
                end

                fit_angle = result.angle
                if whiten_maps !== nothing
                    fit_angle = _unwhiten_angle_at(whiten_maps, ii, jj, fit_angle)
                elseif whiten_w !== nothing
                    fit_angle = _unwhiten_angle(whiten_w, fit_angle)
                end
                # Clamped at zero: a negative centre pixel (noise around zero
                # intensity) would otherwise enter Q/U with a negative weight,
                # i.e. vote for the orientation rotated by 90 degrees. NaN still
                # propagates (max(NaN, 0.0) is NaN) and is masked downstream.
                intensity_weight = max(Float64(mask_image[ii, jj]), 0.0)
                if collect_profile
                    _accumulate_weighted_diagnostic!(fit_rbar2_num, fit_rbar2_w, ii, jj, result.rbar2, intensity_weight)
                    if result.low_confidence
                        fit_low_confidence_tls[tid] += 1
                    end
                end

                passes_rbar2_gate = rbar2_gate <= 0.0 || (isfinite(result.rbar2) && result.rbar2 >= rbar2_gate)
                # isfinite: a window centre without a noise covariance (whitening map NaN there) is not used.
                if result.fit_ok && isfinite(fit_angle) && passes_rbar2_gate
                    stokes_weight = (rbar2_weight && isfinite(result.rbar2)) ? intensity_weight * result.rbar2 : intensity_weight
                    # Accumulate the fitted gradient orientation into Q/U.
                    # The raw gradient is later rotated by 90 degrees to become
                    # the magnetic-field orientation proxy.
                    q_accum[ii, jj] = muladd(stokes_weight, cos(2.0 * fit_angle), q_accum[ii, jj])
                    u_accum[ii, jj] = muladd(stokes_weight, sin(2.0 * fit_angle), u_accum[ii, jj])
                    w_accum[ii, jj] += stokes_weight

                    if collect_profile
                        fit_successes_tls[tid] += 1
                        fit_success_count[ii, jj] += 1
                    end
                elseif collect_profile
                    fit_failures_tls[tid] += 1
                end
            catch err
                is_interrupt_exception(err) && rethrow()
                if collect_profile
                    fit_failures_tls[tid] += 1
                end
            end
        end

        if collect_profile
            profile[:fit_attempts] += sum(fit_attempts_tls)
            profile[:fit_successes] += sum(fit_successes_tls)
            profile[:fit_failures] += sum(fit_failures_tls)
            profile[:valid_windows] += sum(valid_windows_tls)
            profile[:fit_low_confidence_windows] += sum(fit_low_confidence_tls)
            profile[:fit_ns] += sum(fit_ns_tls)
            push!(profile[:window_ns], time_ns() - t_fit_stage)
        end

    end

    if footprint_masked !== nothing
        @inbounds for idx in eachindex(footprint_masked)
            if footprint_masked[idx]
                q_accum[idx] = NaN
                u_accum[idx] = NaN
                # Keep W on the same mask as Q/U so the reliability stays <= 1.
                w_accum[idx] = NaN
            end
        end
    end

    if collect_profile
        fit_rbar2_map = fill(NaN, nx, ny)
        @inbounds for idx in eachindex(fit_rbar2_map)
            fit_rbar2_map[idx] = _weighted_diagnostic_value(fit_rbar2_num, fit_rbar2_w, idx)
        end
        profile[:fit_rbar2] = fit_rbar2_map
        channel_count = max(last_channel - first_channel + 1, 1)
        profile[:fit_attempt_fraction] = fit_attempt_count ./ channel_count
        profile[:fit_support] = fit_success_count ./ channel_count
    end

    return q_accum, u_accum, profile
end

function VGT_QU_mv_2d(cube::AbstractArray{<:Any, 3}, dn, first_channel, last_channel, noise;
                      min_valid::Int=100, fitter_mode::Symbol=:axial_vonmises,
                      rbar2_weight::Bool=false,
                      rbar2_gate::Float64=0.0,
                      gradient_mode::Symbol=:central_checked,
                      gradient_whiten::Bool=false,
                      gradient_noise_cov::AbstractVector{<:Real}=Float64[],
                      gradient_noise_cov_map=nothing,
                      weight_cube=nothing,
                      edge_mask::Bool=false)
    # Public non-profiling VGT entry point used by tests and analysis scripts.
    q_accum, u_accum, _ = _vgt_qu_mv_exact_impl(cube, dn, first_channel, last_channel, noise; min_valid=min_valid,
                                                fitter_mode=fitter_mode, collect_profile=false,
                                                rbar2_weight=rbar2_weight,
                                                rbar2_gate=rbar2_gate,
                                                gradient_mode=gradient_mode,
                                                gradient_whiten=gradient_whiten,
                                                gradient_noise_cov=gradient_noise_cov,
                                                gradient_noise_cov_map=gradient_noise_cov_map,
                                                weight_cube=weight_cube,
                                                edge_mask=edge_mask)
    return q_accum, u_accum
end

function VGT_QU_mv_2d_profiled(cube::AbstractArray{<:Any, 3}, dn, first_channel, last_channel, noise;
                               min_valid::Int=100, fitter_mode::Symbol=:axial_vonmises,
                               rbar2_weight::Bool=false,
                               rbar2_gate::Float64=0.0,
                               gradient_mode::Symbol=:central_checked,
                               gradient_whiten::Bool=false,
                               gradient_noise_cov::AbstractVector{<:Real}=Float64[],
                               gradient_noise_cov_map=nothing,
                               weight_cube=nothing,
                               edge_mask::Bool=false)
    # Same VGT computation as VGT_QU_mv_2d, but also returns diagnostic maps and
    # counters that help evaluate fit quality across the sky, and
    # profile[:fitted_weight], the reliability denominator for
    # finalize_vgt_products.
    return _vgt_qu_mv_exact_impl(cube, dn, first_channel, last_channel, noise; min_valid=min_valid,
                                 fitter_mode=fitter_mode, collect_profile=true,
                                 rbar2_weight=rbar2_weight,
                                 rbar2_gate=rbar2_gate,
                                 gradient_mode=gradient_mode,
                                 gradient_whiten=gradient_whiten,
                                 gradient_noise_cov=gradient_noise_cov,
                                                gradient_noise_cov_map=gradient_noise_cov_map,
                                 weight_cube=weight_cube,
                                 edge_mask=edge_mask)
end

function chunk_y_bounds(j::Int, ny::Int, ny_chunks::Int)
    chunk_y = ceil(Int, ny / ny_chunks)
    y_start = (j - 1) * chunk_y + 1
    y_start > ny && return nothing
    y_end = min(j * chunk_y, ny)
    return y_start, y_end
end

Base.@kwdef struct VGTDiagnostics
    # These counters are science diagnostics, not just logging. They tell us how
    # much of the map was actually supported by successful local fits.
    fit_attempts::Int = 0
    fit_successes::Int = 0
    fit_failures::Int = 0
    valid_windows::Int = 0
    fit_low_confidence_windows::Int = 0
end

function save_chunk_atomic(filename::String, q_accum, u_accum, intensity_accum, fitted_weight,
                           diagnostics::NamedTuple; fit_maps::NamedTuple=(;))
    # Write through a temporary file and then move it into place. If a run is
    # interrupted, this avoids leaving a half-written chunk with the final name.
    tmp = string(filename, ".tmp.", getpid(), ".", Threads.threadid())
    try
        jldsave(tmp;
            q=q_accum,
            u=u_accum,
            i=intensity_accum,
            w=fitted_weight,
            fit_attempts=diagnostics.fit_attempts,
            fit_successes=diagnostics.fit_successes,
            fit_failures=diagnostics.fit_failures,
            valid_windows=diagnostics.valid_windows,
            fit_low_confidence_windows=get(diagnostics, :fit_low_confidence_windows, 0),
            fit_rbar2=get(fit_maps, :fit_rbar2, nothing),
            fit_attempt_fraction=get(fit_maps, :fit_attempt_fraction, nothing),
            fit_support=get(fit_maps, :fit_support, nothing),
        )
        mv(tmp, filename; force=true)
    finally
        if isfile(tmp)
            rm(tmp; force=true)
        end
    end
end

function process_chunk(j::Int, cfg::VGTConfig, ny::Int, nv::Int; noise_model=nothing, channels=1:nv)
    # Large FITS cubes are processed in y chunks. Each chunk includes padding so
    # sliding windows near chunk boundaries still see the same neighbors they
    # would see in a full-cube run.
    y_bounds = chunk_y_bounds(j, ny, cfg.ny_chunks)
    y_bounds === nothing && return nothing
    y_start, y_end = y_bounds
    chunk_start_ns = time_ns()

    f = FITS(cfg.input_file)
    fw = isempty(cfg.weight_file) ? nothing : FITS(cfg.weight_file)
    fcov = (cfg.gradient_whiten && !isempty(cfg.gradient_noise_cov_file)) ? FITS(cfg.gradient_noise_cov_file) : nothing
    try
        pad = chunk_padding_pixels(cfg)
        py_start = max(1, y_start - pad)
        py_end = min(ny, y_end + pad)
        # Only the selected channels are read (velocity_range_kms); nv is their count.
        length(channels) == nv || throw(ArgumentError("channels ($(length(channels))) must have nv=$nv entries"))
        sub_cube = read(f[1], :, py_start:py_end, channels)
        sub_cov_map = fcov === nothing ? nothing : read(fcov[1], :, py_start:py_end, 1:3)
        sub_weight = fw === nothing ? nothing : read(fw[1], :, py_start:py_end, channels)
        resolved_noise_model = noise_model === nothing ? build_noise_model(
            weight_source(cfg),
            cfg.noise,
            cfg.noise_mode;
            tile_size=cfg.noise_tile_size,
            channels=channels,
        ) : noise_model
        noise = noise_for_chunk(
            resolved_noise_model,
            size(sub_cube, 1),
            py_start:py_end,
            nv,
            multiplier=cfg.noise_sigma_multiplier,
        )
        q_padded, u_padded, profile = VGT_QU_mv_2d_profiled(
            sub_cube, cfg.dn, 1, nv, noise;
            min_valid=cfg.min_valid,
            fitter_mode=cfg.fitter_mode,
            rbar2_weight=cfg.rbar2_weight,
            rbar2_gate=cfg.rbar2_gate,
            gradient_mode=cfg.gradient_mode,
            gradient_whiten=cfg.gradient_whiten,
            gradient_noise_cov=cfg.gradient_noise_cov,
            gradient_noise_cov_map=sub_cov_map,
            weight_cube=sub_weight,
            edge_mask=cfg.edge_mask,
        )
        diagnostics = (
            fit_attempts=profile[:fit_attempts],
            fit_successes=profile[:fit_successes],
            fit_failures=profile[:fit_failures],
            valid_windows=profile[:valid_windows],
            fit_low_confidence_windows=profile[:fit_low_confidence_windows],
        )
        fit_maps = (
            fit_rbar2=profile[:fit_rbar2],
            fit_attempt_fraction=profile[:fit_attempt_fraction],
            fit_support=profile[:fit_support],
        )

        trim_y_start = y_start - py_start + 1
        trim_y_end = trim_y_start + (y_end - y_start)
        # Remove padding before saving the chunk. Padding is only for correct
        # window context; it must not duplicate pixels in the stitched output.
        q_accum = q_padded[:, trim_y_start:trim_y_end]
        u_accum = u_padded[:, trim_y_start:trim_y_end]
        w_accum = profile[:fitted_weight][:, trim_y_start:trim_y_end]
        intensity_cube = sub_weight === nothing ? sub_cube : sub_weight
        @views intensity_accum = sum(intensity_cube[:, trim_y_start:trim_y_end, :]; dims=3)[:, :, 1]
        fit_maps_trimmed = (
            fit_rbar2=fit_maps.fit_rbar2[:, trim_y_start:trim_y_end],
            fit_attempt_fraction=fit_maps.fit_attempt_fraction[:, trim_y_start:trim_y_end],
            fit_support=fit_maps.fit_support[:, trim_y_start:trim_y_end],
        )

        filename = joinpath(cfg.output_dir, "chunk_2d_$(j).jld2")
        save_chunk_atomic(filename, q_accum, u_accum, intensity_accum, w_accum, diagnostics; fit_maps=fit_maps_trimmed)
        elapsed_s = round((time_ns() - chunk_start_ns) / 1e9; digits=2)
        @info "Done chunk" chunk=j chunks=cfg.ny_chunks y_start=y_start y_end=y_end elapsed_s=elapsed_s
        return filename
    finally
        close(f)
        fw === nothing || close(fw)
        fcov === nothing || close(fcov)
    end
end

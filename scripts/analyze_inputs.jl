using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO, Statistics, Printf, StatsBase

const NORMAL_MAD = 0.6744897501960817

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

# Estimate 1σ noise from negative-valued pixels. In emission-line cubes,
# negative brightness is normally baseline noise rather than real emission.
#
# Use the median absolute negative amplitude scaled by the standard normal
# median absolute deviation constant. This avoids the low bias from taking
# `std` over only the negative half of a Gaussian noise distribution.
function estimate_noise(sl::AbstractArray)
    vals = filter(v -> isfinite(v) && v < 0.0, vec(sl))
    length(vals) < 100 && return NaN
    return median(abs.(vals)) / NORMAL_MAD
end

# Fallback for cubes with no negative samples. This is less authoritative than
# estimate_noise because clipped cubes have already destroyed half the baseline
# distribution. It is still useful for line-free channels or blank regions whose
# low-intensity values are dominated by noise rather than emission.
function estimate_fallback_noise(sl::AbstractArray)
    vals = Float64[v for v in vec(sl) if isfinite(v)]
    length(vals) < 100 && return NaN

    candidates = Float64[]

    centered_mad = median(abs.(vals .- median(vals))) / NORMAL_MAD
    if isfinite(centered_mad) && centered_mad > 0.0
        push!(candidates, centered_mad)
    end

    positive_vals = Float64[v for v in vals if v > 0.0]
    if length(positive_vals) >= 100
        positive_half = median(positive_vals) / NORMAL_MAD
        if isfinite(positive_half) && positive_half > 0.0
            push!(candidates, positive_half)
        end
    end

    isempty(candidates) && return NaN
    return minimum(candidates)
end

function safe_pct(spread, ref)
    return (isfinite(spread) && isfinite(ref) && ref != 0.0) ? 100 * spread / ref : NaN
end

function collect_patch_ranges(n::Int, patch_size::Int)
    patch_size = max(1, min(n, patch_size))
    if patch_size == n
        return [(1, n)]
    end
    centers = unique(round.(Int, range(cld(patch_size, 2), n - fld(patch_size, 2), length=3)))
    half_lo = fld(patch_size - 1, 2)
    half_hi = cld(patch_size - 1, 2)
    ranges = Tuple{Int,Int}[]
    for c in centers
        lo = max(1, c - half_lo)
        hi = min(n, c + half_hi)
        lo = max(1, hi - patch_size + 1)
        hi = min(n, lo + patch_size - 1)
        push!(ranges, (lo, hi))
    end
    return unique(ranges)
end

function summarize_noise(values::AbstractVector{<:Real})
    clean = sort!(collect(filter(isfinite, values)))
    isempty(clean) && return nothing
    med = median(clean)
    q16 = quantile(clean, 0.16)
    q84 = quantile(clean, 0.84)
    return (
        min=first(clean),
        q16=q16,
        median=med,
        q84=q84,
        max=last(clean),
        spread=q84 - q16,
        spread_pct=safe_pct(q84 - q16, med),
    )
end

function header_str(h, key, default="N/A")
    try
        ;
        return string(h[key]);
    catch
        ;
        return default;
    end
end

# ─────────────────────────────────────────────────────────────────────────────
# Per-file analysis
# ─────────────────────────────────────────────────────────────────────────────

function analyze_file(filepath::String; label::String="")
    println("\n" * "="^72)
    println("FILE: $filepath")
    println("  Label  : $label")
    @printf "  Size   : %.1f MB  (%.2f GB)\n" filesize(filepath)/1024^2 filesize(filepath)/1024^3
    println("="^72)

    f = FITS(filepath)
    h = read_header(f[1])
    dims = size(f[1])
    ndim = length(dims)
    println("\n── Header ────────────────────────────────────────────────────────────")
    println("  NAXIS   : $ndim")
    for ax in 1:ndim
        @printf "  NAXIS%d  : %d\n" ax dims[ax]
    end
    for key in ["BUNIT", "OBJECT", "TELESCOP", "INSTRUME", "RESTFRQ", "RESTFREQ",
        "CTYPE1", "CTYPE2", "CTYPE3", "CDELT1", "CDELT2", "CDELT3",
        "CRVAL1", "CRVAL2", "CRVAL3", "BMAJ", "BMIN", "BPA"]
        v = header_str(h, key)
        v != "N/A" && @printf "  %-10s: %s\n" key v
    end

    if ndim < 3
        println("\n  (2D image — skipping VGT analysis)")
        close(f)
        return
    end

    nx, ny, nv = dims[1], dims[2], dims[3]
    println("\n── Cube dimensions ───────────────────────────────────────────────────")
    @printf "  Spatial : %d × %d  pixels\n" nx ny
    @printf "  Velocity: %d  channels\n" nv

    try
        cdelt3 = h["CDELT3"];
        crval3 = h["CRVAL3"];
        crpix3 = h["CRPIX3"]
        v_start = crval3 + (1 - crpix3) * cdelt3
        v_end = crval3 + (nv - crpix3) * cdelt3
        ctype3 = header_str(h, "CTYPE3", "?")
        @printf "  Vel axis: %.2f → %.2f  (CTYPE3=%s, CDELT3=%g)\n" v_start v_end ctype3 cdelt3
    catch
        ;
    end

    try
        cdelt1 = h["CDELT1"];
        cdelt2 = h["CDELT2"]
        @printf "  Pixel scale: CDELT1=%.6f  CDELT2=%.6f\n" cdelt1 cdelt2
    catch
        ;
    end

    # ── Noise estimation ─────────────────────────────────────────────────────
    println("\n── Noise estimation ──────────────────────────────────────────────────")
    sample_chans = unique(vcat(
        collect(1:max(1, nv÷10)),
        collect(max(1, nv÷2-5):min(nv, nv÷2+5)),
        collect(max(1, nv-nv÷10):nv)
    ))

    noise_ests = Float64[]
    patch_noise_ests = Float64[]
    fallback_noise_ests = Float64[]
    fallback_patch_noise_ests = Float64[]
    global_min = Inf;
    global_max = -Inf
    n_finite = 0;
    n_total = 0

    patch_size_x = max(32, min(128, nx ÷ 4))
    patch_size_y = max(32, min(128, ny ÷ 4))
    xranges = collect_patch_ranges(nx, patch_size_x)
    yranges = collect_patch_ranges(ny, patch_size_y)

    for k in sample_chans
        sl = read(f[1],:,:,(k:k))[:, :, 1]
        ne = estimate_noise(sl)
        if !isnan(ne)
            ;
            push!(noise_ests, ne);
        end
        fne = estimate_fallback_noise(sl)
        if !isnan(fne)
            ;
            push!(fallback_noise_ests, fne);
        end
        for (x1, x2) in xranges, (y1, y2) in yranges
            patch_view = @view sl[x1:x2, y1:y2]
            pne = estimate_noise(patch_view)
            if !isnan(pne)
                ;
                push!(patch_noise_ests, pne);
            end
            fpne = estimate_fallback_noise(patch_view)
            if !isnan(fpne)
                ;
                push!(fallback_patch_noise_ests, fpne);
            end
        end
        n_finite += count(isfinite, sl)
        n_total += length(sl)
        mn, mx = extrema(filter(isfinite, vec(sl)))
        global_min = min(global_min, mn)
        global_max = max(global_max, mx)
    end

    used_noise_ests = noise_ests
    used_patch_noise_ests = patch_noise_ests
    used_fallback_noise = false

    if !isempty(noise_ests)
        noise_est = median(noise_ests)
        @printf "  Noise threshold (robust, from negative pixels): %.5f  (sampled %d channels)\n" noise_est length(noise_ests)
        @printf "  Suggested config: noise = %.2f\n" noise_est
    elseif !isempty(fallback_noise_ests)
        used_noise_ests = fallback_noise_ests
        used_patch_noise_ests = fallback_patch_noise_ests
        used_fallback_noise = true
        noise_est = median(fallback_noise_ests)
        @printf "  Noise threshold fallback (line-free/blank robust scatter): %.5f  (sampled %d channels)\n" noise_est length(fallback_noise_ests)
        @printf "  Suggested config: noise = %.2f\n" noise_est
        println("  WARNING: no negative pixels were available; this fallback is weaker than")
        println("           negative-pixel noise. Prefer an unclipped calibrated cube when possible.")
    else
        noise_est = NaN
        println("  Noise: could not estimate (no negative pixels in sample)")
    end
    @printf "  Data range (sampled): %.4f  to  %.4f\n" global_min global_max
    @printf "  Finite values in sampled channels: %.1f%%\n" 100*n_finite/n_total
    if used_fallback_noise && isfinite(global_min) && global_min >= 0.0
        println("  Sampled minimum is nonnegative; cube may be clipped or baseline-shifted.")
    end

    chan_summary = summarize_noise(used_noise_ests)
    if chan_summary !== nothing
        label = used_fallback_noise ? "fallback noise spread" : "noise spread"
        println("  Channel-to-channel $label:")
        @printf "    min / p16 / median / p84 / max = %.5f / %.5f / %.5f / %.5f / %.5f\n" chan_summary.min chan_summary.q16 chan_summary.median chan_summary.q84 chan_summary.max
        @printf "    central spread (p84-p16)      = %.5f  (%.1f%% of median)\n" chan_summary.spread chan_summary.spread_pct
    end

    patch_summary = summarize_noise(used_patch_noise_ests)
    if patch_summary !== nothing
        label = used_fallback_noise ? "fallback noise spread" : "noise spread"
        println("  Spatial patch $label:")
        @printf "    sampled %d patches (%d×%d over %d channels)\n" length(used_patch_noise_ests) length(xranges) length(yranges) length(sample_chans)
        @printf "    min / p16 / median / p84 / max = %.5f / %.5f / %.5f / %.5f / %.5f\n" patch_summary.min patch_summary.q16 patch_summary.median patch_summary.q84 patch_summary.max
        @printf "    central spread (p84-p16)      = %.5f  (%.1f%% of median)\n" patch_summary.spread patch_summary.spread_pct
    end

    recommend_noise_map = false
    if chan_summary !== nothing && isfinite(chan_summary.spread_pct) && chan_summary.spread_pct > 35
        recommend_noise_map = true
    end
    if patch_summary !== nothing && isfinite(patch_summary.spread_pct) && patch_summary.spread_pct > 35
        recommend_noise_map = true
    end
    println("  Recommendation:")
    if recommend_noise_map
        println("    Noise varies materially across sampled channels/patches.")
        println("    A scalar noise value may mis-mask parts of the cube; investigate noise_mode=\"auto\" or \"separable\".")
    else
        println("    Noise looks stable enough for a scalar config value.")
        println("    Keep `noise` constant unless science validation shows edge/coverage artifacts.")
    end

    thresh_use = isnan(noise_est) ? 0.0 : noise_est

    # ── Channel analysis ─────────────────────────────────────────────────────
    println("\n── Channel analysis ──────────────────────────────────────────────────")
    print("  Scanning all $nv channels ... ")
    chan_sum = zeros(nv)
    chan_filled = zeros(nv)

    for k in 1:nv
        sl = read(f[1],:,:,(k:k))[:, :, 1]
        above = count(v -> isfinite(v) && v > thresh_use, sl)
        chan_filled[k] = above / length(sl)
        chan_sum[k] = sum(v -> (isfinite(v) && v > thresh_use) ? v : 0.0, sl)
    end
    k_peak = argmax(chan_sum)
    k_lo = max(1, k_peak - 20);
    k_hi = min(nv, k_peak + 20)
    @printf "done\n  Peak emission channel: %d  (analysis band %d:%d)\n" k_peak k_lo k_hi
    @printf "  Fill at peak (>noise): %.4f%%\n" 100*chan_filled[k_peak]

    nonzero_chans = count(>(0), chan_filled)
    @printf "  Channels with any emission: %d / %d  (%.1f%%)\n" nonzero_chans nv 100*nonzero_chans/nv
    if nonzero_chans > 0
        cmin = 100*minimum(filter(>(0), chan_filled))
        cmed = 100*median(filter(>(0), chan_filled))
        cmax = 100*maximum(chan_filled)
        @printf "  Fill stats across emission channels: min=%.4f%%  median=%.4f%%  max=%.4f%%\n" cmin cmed cmax
    end

    # ── Integrated intensity map ──────────────────────────────────────────────
    println("\n── Integrated intensity map ──────────────────────────────────────────")
    print("  Building imap over channels $k_lo:$k_hi ... ")
    imap = zeros(Float32, nx, ny)
    for k in k_lo:k_hi
        sl = read(f[1],:,:,(k:k))[:, :, 1]
        @inbounds for i in eachindex(imap)
            v = sl[i]
            if isfinite(v) && v > thresh_use
                ;
                imap[i] += v;
            end
        end
    end
    n_nonzero = count(>(0), imap)
    @printf "done\n  Non-zero pixels: %d / %d  (%.2f%%)\n" n_nonzero length(imap) 100*n_nonzero/length(imap)

    # ── Spatial patch: density table for window fit viability ──────────────
    patch_size = min(400, nx÷2, ny÷2)
    half_p = patch_size ÷ 2
    println("\n── window fit viability (densest patch, size=$(patch_size)×$(patch_size)) ──────────")

    step_p = max(1, patch_size ÷ 4)
    best_sum = -1.0;
    best_cx = half_p;
    best_cy = half_p
    for cx in half_p:step_p:(nx-half_p), cy in half_p:step_p:(ny-half_p)
        s = sum(imap[(cx-half_p+1):(cx+half_p), (cy-half_p+1):(cy+half_p)])
        if s > best_sum
            ;
            best_sum = s;
            best_cx = cx;
            best_cy = cy;
        end
    end
    x1 = best_cx-half_p+1;
    x2 = min(nx, best_cx+half_p)
    y1 = best_cy-half_p+1;
    y2 = min(ny, best_cy+half_p)
    @printf "  Best patch: [%d:%d, %d:%d]\n" x1 x2 y1 y2

    print("  Reading patch ... ")
    patch = read(f[1], x1:x2, y1:y2, k_lo:k_hi)
    close(f)
    pnx, pny, pnv = size(patch)
    println("done ($(pnx)×$(pny)×$(pnv))")

    above = count(v -> isfinite(v) && v > thresh_use, patch)
    @printf "  Voxels above noise in patch: %.3f%%\n" 100*above/length(patch)

    k_local_peak = argmax([count(v -> isfinite(v) && v > thresh_use, @view patch[:, :, k]) for k in 1:pnv])
    peak_slice = @view patch[:, :, k_local_peak]
    @printf "  Peak channel in patch: local=%d  global=%d\n" k_local_peak (k_lo + k_local_peak - 1)

    println("\n  Valid pixels (>noise) in densest dn×dn window at peak channel:")
    println("  (need at least min_valid (default 100) valid gradients per window; dn in config.toml controls spatial resolution)")
    @printf "  %-6s  %-8s  %-8s  %-6s  %s\n" "dn" "window" "n_valid" "fill%" "verdict"
    dn_ok = Int[]
    for test_dn in [5, 10, 15, 20, 30, 50, 75, 100]
        h_dn = div(test_dn, 2)
        h_dn < 1 && continue
        (h_dn > pnx-1 || h_dn > pny-1) && continue
        best_v = 0
        for j in h_dn:(pny-h_dn), i in h_dn:(pnx-h_dn)
            v = count(x -> isfinite(x) && x > thresh_use, @view peak_slice[(i-h_dn+1):(i+h_dn), (j-h_dn+1):(j+h_dn)])
            if v > best_v
                ;
                best_v = v;
            end
        end
        win = (2*h_dn)^2
        fill_pct = 100.0 * best_v / win
        verdict = best_v >= 100 ? "✓ window fit OK" : (best_v >= 10 ? "~ marginal" : "✗ too sparse")
        if best_v >= 100
            ;
            push!(dn_ok, test_dn);
        end
        @printf "  dn=%-3d  win=%-6d  n=%-5d  %5.1f%%  %s\n" test_dn win best_v fill_pct verdict
    end

    # ── Parameter recommendations ─────────────────────────────────────────────
    println("\n── Parameter recommendations ─────────────────────────────────────────")
    if !isnan(noise_est)
        @printf "  noise       = %.2f\n" noise_est
    end
    if !isempty(dn_ok)
        @printf "  dn        = %d  (smallest giving ≥100 valid px; larger = smoother output)\n" minimum(dn_ok)
        @printf "  min_valid = 100\n"
    else
        println("  WARNING: no tested dn reaches 100 valid pixels in the densest window.")
        println("  The noise threshold may need lowering, or the data is genuinely too sparse for VGT.")
    end
    println()
end

function usage()
    println("Usage:")
    println("  julia scripts/analyze_inputs.jl [fits_path ...]")
    println()
    println("If no FITS paths are provided, the script analyzes the built-in default file list.")
end

function default_files()
    return [
        (datadir("exp_raw", "Tb_combined_mw_full_hnv_res1_3pperb_reprojected_on_gass.fits"),
        "CO cube — primary VGT input"),
    ]
end

function parse_args(args)
    if any(arg -> arg in ("-h", "--help"), args)
        usage()
        return nothing
    end

    isempty(args) && return default_files()

    files = Tuple{String,String}[]
    for path in args
        startswith(path, "--") && error("Unknown option: $path")
        push!(files, (path, basename(path)))
    end
    return files
end

function main(args)
    files = parse_args(args)
    files === nothing && return 0

    for (path, lbl) in files
        if isfile(path)
            analyze_file(path; label=lbl)
        else
            println("\nFILE NOT FOUND: $path")
        end
    end

    println("\n" * "="^72)
    println("Analysis complete.")
    println("="^72)
    return 0
end

if abspath(PROGRAM_FILE) == @__FILE__
    exit(main(ARGS))
end

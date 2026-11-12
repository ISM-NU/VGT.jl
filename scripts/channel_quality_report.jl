using DrWatson
quickactivate(@__DIR__, "VGT")

using FITSIO
using Printf
using Statistics
using VGT

const MAD_NORMAL = 0.6744897501960817

function default_options()
    return Dict(
        "input" => nothing,
        "noise" => nothing,
        "output-dir" => experiment_default_dir("channel_quality"),
        "min-fill-3sigma" => "0.001",
        "min-positive-negative-ratio" => "1.25",
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
        i == length(args) && throw(ArgumentError("Missing value for --$key"))
        opts[key] = args[i + 1]
        i += 2
    end
    opts["input"] === nothing && throw(ArgumentError("--input is required"))
    return opts
end

function robust_negative_noise(values; min_count::Int=100)
    vals = Float64[]
    for v in values
        if isfinite(v) && v < 0
            push!(vals, abs(Float64(v)))
        end
    end
    length(vals) < min_count && return NaN
    return median(vals) / MAD_NORMAL
end

function velocity_axis_kms(header, nv::Int)
    crval = Float64(header["CRVAL3"])
    crpix = Float64(header["CRPIX3"])
    cdelt = Float64(header["CDELT3"])
    native = [crval + (k - crpix) * cdelt for k in 1:nv]
    scale = maximum(abs.(native)) > 1_000 ? 1e-3 : 1.0
    return native .* scale
end

function finite_positive_median(values, fallback::Real)
    clean = [Float64(v) for v in values if isfinite(v) && v > 0]
    isempty(clean) && return Float64(fallback)
    return median(clean)
end

function write_csv(path::AbstractString, rows, columns)
    mkpath(dirname(path))
    open(path, "w") do io
        println(io, join(columns, ","))
        for row in rows
            println(io, join([string(get(row, col, "")) for col in columns], ","))
        end
    end
end

function channel_quality_report(opts)
    input = opts["input"]
    min_fill_3sigma = parse(Float64, opts["min-fill-3sigma"])
    min_pos_neg_ratio = parse(Float64, opts["min-positive-negative-ratio"])
    f = FITS(input)
    try
        nx, ny, nv = size(f[1])
        header = read_header(f[1])
        velocities = velocity_axis_kms(header, nv)
        fallback_noise = opts["noise"] === nothing ? NaN : parse(Float64, opts["noise"])
        channel_noise = fill(NaN, nv)
        rows = Vector{Dict{String, Any}}()

        for k in 1:nv
            sl = Float64.(read(f[1], :, :, k:k)[:, :, 1])
            noise = robust_negative_noise(sl)
            if !(isfinite(noise) && noise > 0) && isfinite(fallback_noise)
                noise = fallback_noise
            end
            channel_noise[k] = noise
            finite_count = count(isfinite, sl)
            fill1 = isfinite(noise) ? count(v -> isfinite(v) && v > noise, sl) / length(sl) : NaN
            fill2 = isfinite(noise) ? count(v -> isfinite(v) && v > 2noise, sl) / length(sl) : NaN
            fill3 = isfinite(noise) ? count(v -> isfinite(v) && v > 3noise, sl) / length(sl) : NaN
            positive_sum = sum(v -> isfinite(v) && v > 0 ? Float64(v) : 0.0, sl)
            negative_abs_sum = sum(v -> isfinite(v) && v < 0 ? abs(Float64(v)) : 0.0, sl)
            pos_neg_ratio = negative_abs_sum > 0 ? positive_sum / negative_abs_sum : Inf
            vals = [Float64(v) for v in vec(sl) if isfinite(v)]
            p99 = isempty(vals) ? NaN : quantile(vals, 0.99)
            max_value = isempty(vals) ? NaN : maximum(vals)
            no_signal = positive_sum == 0.0 && negative_abs_sum == 0.0 && max_value == 0.0
            likely_noise_only = !no_signal &&
                                isfinite(fill3) &&
                                fill3 < min_fill_3sigma &&
                                pos_neg_ratio < min_pos_neg_ratio
            push!(rows, Dict{String, Any}(
                "channel" => k,
                "velocity_kms" => velocities[k],
                "noise" => noise,
                "finite_fraction" => finite_count / length(sl),
                "fill_1sigma" => fill1,
                "fill_2sigma" => fill2,
                "fill_3sigma" => fill3,
                "positive_sum" => positive_sum,
                "negative_abs_sum" => negative_abs_sum,
                "positive_negative_ratio" => pos_neg_ratio,
                "p99" => p99,
                "max" => max_value,
                "no_signal" => no_signal,
                "likely_noise_only" => likely_noise_only,
            ))
        end

        scalar_noise = isfinite(fallback_noise) ? fallback_noise :
                       finite_positive_median(channel_noise, NaN)
        output_dir = opts["output-dir"]
        mkpath(output_dir)
        stem = splitext(basename(input))[1]
        csv_path = joinpath(output_dir, string("channel_quality_", stem, ".csv"))
        columns = [
            "channel", "velocity_kms", "noise", "finite_fraction",
            "fill_1sigma", "fill_2sigma", "fill_3sigma",
            "positive_sum", "negative_abs_sum", "positive_negative_ratio",
            "p99", "max", "no_signal", "likely_noise_only",
        ]
        write_csv(csv_path, rows, columns)

        likely = [row for row in rows if row["likely_noise_only"]]
        no_signal_rows = [row for row in rows if row["no_signal"]]
        low_fill = [row for row in rows if isfinite(row["fill_3sigma"]) &&
                    row["fill_3sigma"] < min_fill_3sigma]
        @printf "Wrote %s\n" csv_path
        @printf "Channels: %d, scalar reference noise %.5g\n" nv scalar_noise
        @printf "Likely noise-only channels: %d\n" length(likely)
        if !isempty(likely)
            println("  ", join([@sprintf("%d(%.2f km/s)", row["channel"], row["velocity_kms"]) for row in likely], ", "))
        end
        @printf "No-signal zero channels: %d\n" length(no_signal_rows)
        if !isempty(no_signal_rows)
            println("  ", join([@sprintf("%d(%.2f km/s)", row["channel"], row["velocity_kms"]) for row in no_signal_rows], ", "))
        end
        @printf "Channels with fill_3sigma < %.4f: %d\n" min_fill_3sigma length(low_fill)
        return csv_path
    finally
        close(f)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    try
        channel_quality_report(parse_args(ARGS))
    catch err
        showerror(stderr, err)
        println(stderr)
        exit(1)
    end
end

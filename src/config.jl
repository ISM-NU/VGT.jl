# Configuration and execution settings for a VGT run.
#
# Most fields are not generic software options; they encode scientific choices
# such as window size, noise model, fitter mode, and masking policy.

const ANGLE_FRAMES = (:pixel, :galactic, :icrs)

function normalize_angle_frame(value)
    frame = _symbol_config_value(value)
    frame in ANGLE_FRAMES && return frame
    throw(ArgumentError("angle_frame must be one of $(join(ANGLE_FRAMES, ", ")); got $(value)"))
end

struct VGTConfig
    input_file::String
    output_dir::String
    keep_intermediate_chunks::Bool
    noise::Float64
    dn::Int
    ker::Int
    ny_chunks::Union{Int, Symbol}
    min_valid::Int
    fitter_mode::Symbol
    rbar2_weight::Bool
    noise_mode::Symbol
    noise_tile_size::Int
    zero_stokes_policy::Symbol
    output_support_mode::Symbol
    gradient_mode::Symbol
    noise_sigma_multiplier::Float64
    rbar2_gate::Float64
    # Optional cube for the noise mask, Stokes weights and intensity map while the
    # gradients come from input_file (VDA: gradients of the velocity caustics,
    # weighted by the observed intensity). Empty = input_file.
    weight_file::String
    # Frame of the written psi map: :pixel (counter-clockwise from +FITS axis 1,
    # the frame VGT computes in) or an IAU position angle, east of north, in
    # :galactic or :icrs (converted per pixel through the WCS after smoothing).
    angle_frame::Symbol
    # Output grid: 1 keeps the input grid; b > 1 averages the smoothed Q/U/I over
    # b x b input-pixel blocks before psi and reliability are computed, and
    # scales the celestial WCS to match (see `bin_products`).
    output_bin::Int
    # Velocity window [lo, hi] in km/s: only channels whose centres lie inside it
    # (inclusive) are used, as if the cube were cut to them. Empty = all channels.
    velocity_range_kms::Vector{Float64}

    function VGTConfig(input_file::String, output_dir::String,
                       keep_intermediate_chunks::Bool, noise::Float64, dn::Int, ker::Int,
                       ny_chunks::Union{Int, Symbol},
                       min_valid::Int,
                       fitter_mode::Symbol,
                       rbar2_weight::Bool,
                       noise_mode::Symbol, noise_tile_size::Int,
                       zero_stokes_policy::Symbol, output_support_mode::Symbol,
                       gradient_mode::Symbol,
                       noise_sigma_multiplier::Float64,
                       rbar2_gate::Float64,
                       weight_file::String="",
                       angle_frame::Symbol=:pixel,
                       output_bin::Int=1,
                       velocity_range_kms::Vector{Float64}=Float64[])
        # Validate here so invalid science settings fail before a large FITS
        # cube run starts.
        _validate_config_fields(noise, dn, ker, min_valid, fitter_mode,
                                noise_mode, noise_tile_size,
                                noise_sigma_multiplier, zero_stokes_policy,
                                output_support_mode, gradient_mode, rbar2_gate)
        angle_frame in ANGLE_FRAMES ||
            throw(ArgumentError("angle_frame must be one of $(join(ANGLE_FRAMES, ", ")); got $(angle_frame)"))
        output_bin >= 1 || throw(ArgumentError("output_bin must be a positive integer; got $(output_bin)"))
        isempty(velocity_range_kms) || (length(velocity_range_kms) == 2 && all(isfinite, velocity_range_kms) &&
                                        velocity_range_kms[1] < velocity_range_kms[2]) ||
            throw(ArgumentError("velocity_range_kms must be empty or [lo, hi] with lo < hi; got $(velocity_range_kms)"))
        return new(input_file, output_dir, keep_intermediate_chunks, noise, dn, ker,
                   ny_chunks, min_valid, fitter_mode,
                   rbar2_weight, noise_mode,
                   noise_tile_size, zero_stokes_policy, output_support_mode,
                   gradient_mode, noise_sigma_multiplier,
                   rbar2_gate, weight_file, angle_frame, output_bin, velocity_range_kms)
    end
end

function VGTConfig(input_file::AbstractString, output_dir::AbstractString,
                   keep_intermediate_chunks::Bool, noise, dn, ker,
                   ny_chunks=1, min_valid=100,
                   fitter_mode=:axial_vonmises,
                   rbar2_weight=false, noise_mode=:scalar,
                   noise_tile_size=128, zero_stokes_policy=:nan_before_smoothing,
                   output_support_mode=:smoothed_stokes_support,
                   gradient_mode=:central_checked,
                   noise_sigma_multiplier=1.0,
                   rbar2_gate=0.0,
                   weight_file="",
                   angle_frame=:pixel,
                   output_bin=1,
                   velocity_range_kms=Float64[])
    return VGTConfig(
        String(input_file),
        String(output_dir),
        keep_intermediate_chunks,
        Float64(noise),
        Int(dn),
        Int(ker),
        normalize_execution_setting(ny_chunks, :ny_chunks),
        Int(min_valid),
        normalize_fitter_mode(fitter_mode),
        Bool(rbar2_weight),
        normalize_noise_mode(noise_mode),
        Int(noise_tile_size),
        normalize_zero_stokes_policy(zero_stokes_policy),
        normalize_output_support_mode(output_support_mode),
        normalize_gradient_mode(gradient_mode),
        Float64(noise_sigma_multiplier),
        Float64(rbar2_gate),
        String(weight_file),
        normalize_angle_frame(angle_frame),
        Int(output_bin),
        collect(Float64, velocity_range_kms),  # Float64.(x) of an empty Vector{Any} stays Vector{Any}
    )
end

function VGTConfig(; input_file::AbstractString, output_dir::AbstractString="subcube",
                   keep_intermediate_chunks::Bool=true, noise=nothing, dn, ker,
                   ny_chunks=1, n_threads=nothing, min_valid=100,
                   fitter_mode=:axial_vonmises,
                   rbar2_weight=false, noise_mode=:scalar,
                   noise_tile_size=128,
                   zero_stokes_policy=:nan_before_smoothing,
                   output_support_mode=:smoothed_stokes_support,
                   gradient_mode=:central_checked,
                   noise_sigma_multiplier=1.0,
                   rbar2_gate=0.0,
                   weight_file="",
                   angle_frame=:pixel,
                   output_bin=1,
                   velocity_range_kms=Float64[],
                   # Removed in 0.1.0; only the old default value is tolerated.
                   finalize_mask_mode=nothing,
                   histogram_edge_mode=nothing,
                   gradient_sample_mode=nothing,
                   noise_spatial_exponent=nothing,
                   noise_channel_exponent=nothing)
    _check_obsolete_keys(; finalize_mask_mode, histogram_edge_mode, gradient_sample_mode,
                         noise_spatial_exponent, noise_channel_exponent)
    noise === nothing && throw(ArgumentError("Missing required keyword `noise`"))
    n_threads === nothing ||
        @warn "`n_threads` no longer controls execution; start Julia with `--threads=N` or set JULIA_NUM_THREADS" maxlog=1
    return VGTConfig(
        String(input_file),
        String(output_dir),
        keep_intermediate_chunks,
        Float64(noise),
        Int(dn),
        Int(ker),
        normalize_execution_setting(ny_chunks, :ny_chunks),
        Int(min_valid),
        normalize_fitter_mode(fitter_mode),
        Bool(rbar2_weight),
        normalize_noise_mode(noise_mode),
        Int(noise_tile_size),
        normalize_zero_stokes_policy(zero_stokes_policy),
        normalize_output_support_mode(output_support_mode),
        normalize_gradient_mode(gradient_mode),
        Float64(noise_sigma_multiplier),
        Float64(rbar2_gate),
        String(weight_file),
        normalize_angle_frame(angle_frame),
        Int(output_bin),
        collect(Float64, velocity_range_kms),  # Float64.(x) of an empty Vector{Any} stays Vector{Any}
    )
end

# Copy of `cfg` with some fields replaced. Copies every field by name, so a new
# field cannot be dropped by a hand-written field list.
function config_with(cfg::VGTConfig; changes...)
    values = (haskey(changes, f) ? changes[f] : getfield(cfg, f) for f in fieldnames(VGTConfig))
    return VGTConfig(values...)
end

function _symbol_config_value(value)
    value isa Symbol && return value
    if value isa AbstractString
        text = String(value)
        startswith(text, ":") && (text = text[2:end])
        return Symbol(text)
    end
    return Symbol(value)
end

function normalize_execution_setting(value, name::Symbol)
    (value === :auto || value == "auto" || value == ":auto") && return :auto
    value isa Integer && value > 0 && return Int(value)
    throw(ArgumentError("$(name) must be a positive integer or :auto, got $(repr(value))"))
end

const REMOVED_IN_0_1_0 = "removed in VGT 0.1.0 (public release) because validation found it worse than the default; it remains in the archived private repository"

# Removed values of kept options, so a stale config gets a clear message
# instead of a generic "unknown value" error.
const REMOVED_OPTION_VALUES = Dict(
    :fitter_mode => (:lsqfit_reference, :legacy_lsqfit_reference),
    :gradient_mode => (:legacy_circshift,),
    :noise_mode => (:separable_shrinkage,),
)

function _normalize_mode(value, name::Symbol, allowed::Tuple)
    mode = _symbol_config_value(value)
    mode in allowed && return mode
    mode in get(REMOVED_OPTION_VALUES, name, ()) &&
        throw(ArgumentError("$(name)=$(repr(mode)) was $REMOVED_IN_0_1_0"))
    throw(ArgumentError("$(name) must be one of $(join(repr.(allowed), ", ", ", or ")); got $(repr(value))"))
end

normalize_fitter_mode(value) = _normalize_mode(value, :fitter_mode, (:axial_vonmises,))
normalize_noise_mode(value) =
    _normalize_mode(value, :noise_mode, (:scalar, :channel, :spatial, :separable, :auto))
normalize_gradient_mode(value) =
    _normalize_mode(value, :gradient_mode, (:central_checked, :sobel_checked, :nan_edge_backward))

# Config keys removed in 0.1.0 and the only value each still tolerates (its old
# default, which changes nothing). Any other value would have changed the science.
const OBSOLETE_CONFIG_KEYS = Dict(
    :finalize_mask_mode => :strict_nonpositive,
    :histogram_edge_mode => :modern_range,
    :gradient_sample_mode => :paired_finite,
    :noise_spatial_exponent => 1.0,
    :noise_channel_exponent => 1.0,
)

function _check_obsolete_keys(; kwargs...)
    for (key, value) in kwargs
        value === nothing && continue
        default = OBSOLETE_CONFIG_KEYS[key]
        given = default isa Symbol ? _obsolete_symbol(value) : value
        if given == default
            @warn "Config key `$key` is obsolete and ignored (value $(repr(value)) is the old default)"
        else
            throw(ArgumentError("$key=$(repr(value)) was $REMOVED_IN_0_1_0"))
        end
    end
    return nothing
end

_obsolete_symbol(value) = value isa Union{Symbol, AbstractString} ? _symbol_config_value(value) : value

function normalize_zero_stokes_policy(value)
    mode = _symbol_config_value(value)
    mode in (:nan_before_smoothing, :keep_zero_before_smoothing) && return mode
    throw(ArgumentError("zero_stokes_policy must be one of :nan_before_smoothing or :keep_zero_before_smoothing; got $(repr(value))"))
end

function normalize_output_support_mode(value)
    mode = _symbol_config_value(value)
    mode in (:smoothed_stokes_support, :direct_fit_support) && return mode
    throw(ArgumentError("output_support_mode must be one of :smoothed_stokes_support or :direct_fit_support; got $(repr(value))"))
end

function _validate_config_fields(noise, dn, ker, min_valid, fitter_mode,
                                 noise_mode, noise_tile_size,
                                 noise_sigma_multiplier, zero_stokes_policy,
                                 output_support_mode, gradient_mode,
                                 rbar2_gate=0.0)
    noise > 0 || throw(ArgumentError("noise must be > 0"))
    dn >= 2 || throw(ArgumentError("dn must be >= 2"))
    iseven(dn) || throw(ArgumentError("dn must be even; got $dn"))
    ker >= 1 || throw(ArgumentError("ker must be >= 1"))
    min_valid >= 1 || throw(ArgumentError("min_valid must be >= 1"))
    min_valid <= dn * dn ||
        throw(ArgumentError("min_valid=$min_valid exceeds the maximum window area dn^2=$(dn * dn)"))
    fitter_mode === :axial_vonmises ||
        throw(ArgumentError("Unsupported fitter_mode=$fitter_mode"))
    noise_mode in (:scalar, :channel, :spatial, :separable, :auto) ||
        throw(ArgumentError("Unsupported noise_mode=$noise_mode"))
    noise_tile_size >= 1 || throw(ArgumentError("noise_tile_size must be >= 1"))
    noise_sigma_multiplier > 0 ||
        throw(ArgumentError("noise_sigma_multiplier must be > 0"))
    zero_stokes_policy in (:nan_before_smoothing, :keep_zero_before_smoothing) ||
        throw(ArgumentError("Unsupported zero_stokes_policy=$zero_stokes_policy"))
    output_support_mode in (:smoothed_stokes_support, :direct_fit_support) ||
        throw(ArgumentError("Unsupported output_support_mode=$output_support_mode"))
    gradient_mode in (:central_checked, :sobel_checked, :nan_edge_backward) ||
        throw(ArgumentError("Unsupported gradient_mode=$gradient_mode"))
    0.0 <= rbar2_gate <= 1.0 ||
        throw(ArgumentError("rbar2_gate must be in [0, 1], got $rbar2_gate"))
    return nothing
end

function validate_static_config(cfg::VGTConfig)
    _validate_config_fields(cfg.noise, cfg.dn, cfg.ker, cfg.min_valid,
                            cfg.fitter_mode, cfg.noise_mode, cfg.noise_tile_size,
                            cfg.noise_sigma_multiplier, cfg.zero_stokes_policy,
                            cfg.output_support_mode, cfg.gradient_mode,
                            cfg.rbar2_gate)
    return cfg
end

function chunk_padding_pixels(cfg::VGTConfig)
    # Half-window padding is the exact context required by the current dn x dn
    # local VGT window. Final smoothing is applied only after chunk stitching.
    return div(cfg.dn, 2)
end

function choose_auto_threads(cpu_threads::Int)
    return max(1, cpu_threads)
end

function choose_auto_chunks(ny::Int, dn::Int, n_threads::Int)
    target = min(max(1, 2 * n_threads), 32)
    # dn is the correctness floor (process_chunk pads by dn/2 on each side);
    # The quality floor of 32 avoids pathologically thin strips.
    min_chunk_height = max(dn, 32)
    max_reasonable = max(1, cld(ny, min_chunk_height))
    return max(1, min(target, max_reasonable, ny))
end

# Estimate peak RAM per chunk: sub-cube (Float32) + working channel arrays (Float64).
# Sub-cube: nx * strip_height * nv * 4 bytes
# Channel working arrays (dsx, dsy, q_full, u_full, v_full): 5 * nx * strip_height * 8 bytes
function estimate_peak_chunk_bytes(nx::Int, nv::Int, strip_height::Int)
    # This is a conservative memory estimate used only to choose chunk counts.
    subcube = nx * strip_height * nv * 4
    channel_arrays = 5 * nx * strip_height * 8
    return subcube + channel_arrays
end

function choose_auto_execution(cfg::VGTConfig, nx::Int, ny::Int, nv::Int,
                               cpu_threads::Int, budget_bytes::Int)
    # Auto execution chooses threads and y chunks without changing the science
    # parameters. It only decides how to split the same computation in memory.
    n_threads = choose_auto_threads(cpu_threads)
    ny_chunks = cfg.ny_chunks === :auto ? choose_auto_chunks(ny, cfg.dn, n_threads) : cfg.ny_chunks
    ny_chunks = max(1, min(ny_chunks, ny))

    while ny_chunks < ny
        strip_height = cld(ny, ny_chunks) + 2 * chunk_padding_pixels(cfg)
        estimate_peak_chunk_bytes(nx, nv, strip_height) <= budget_bytes && break
        ny_chunks += 1
    end

    strip_height = cld(ny, ny_chunks) + 2 * chunk_padding_pixels(cfg)
    estimate_peak_chunk_bytes(nx, nv, strip_height) <= budget_bytes ||
        throw(ArgumentError("Unable to resolve a feasible chunk configuration for cube size ($nx, $ny, $nv)"))
    return (; julia_threads=n_threads, ny_chunks)
end

# Half of the memory this process may use. Under Slurm that is the job's
# allocation (SLURM_MEM_PER_NODE, in MB), not the node's RAM: a budget based on
# the whole node lets chunks outgrow the allocation and the job gets killed.
function memory_budget_bytes(env=ENV)
    total = Int(Sys.total_memory())
    slurm_mb = tryparse(Int, get(env, "SLURM_MEM_PER_NODE", ""))
    slurm_mb === nothing || (total = min(total, slurm_mb * 1024^2))
    return max(1, fld(total, 2))
end

# Cube that sets the noise mask, Stokes weights and intensity map.
weight_source(cfg::VGTConfig) = isempty(cfg.weight_file) ? cfg.input_file : cfg.weight_file

function resolve_execution_config(cfg::VGTConfig, nx::Int, ny::Int, nv::Int=0)
    validate_static_config(cfg)

    cfg.dn <= nx || throw(ArgumentError("dn=$(cfg.dn) is too large for nx=$nx (need nx ≥ $(cfg.dn))"))
    cfg.dn <= ny || throw(ArgumentError("dn=$(cfg.dn) is too large for ny=$ny (need ny ≥ $(cfg.dn))"))
    cfg.output_bin <= min(nx, ny) ||
        throw(ArgumentError("output_bin=$(cfg.output_bin) is larger than the map ($nx x $ny pixels)"))

    cpu_threads = max(1, Threads.nthreads())
    budget = memory_budget_bytes()
    auto_execution = choose_auto_execution(cfg, nx, ny, nv, cpu_threads, budget)
    ny_chunks = auto_execution.ny_chunks

    return config_with(cfg; ny_chunks=ny_chunks)
end

const VGT_CONFIG_TABLE_KEYS = Dict(
    "input" => Set(["input_file", "output_dir", "keep_intermediate_chunks", "weight_file", "velocity_range_kms"]),
    "noise" => Set([
        "noise", "noise_mode", "noise_tile_size", "noise_sigma_multiplier",
        "noise_spatial_exponent", "noise_channel_exponent",  # obsolete, see OBSOLETE_CONFIG_KEYS
    ]),
    "vgt" => Set([
        "dn", "ker", "min_valid", "fitter_mode",
        "zero_stokes_policy", "output_support_mode", "gradient_mode",
        "rbar2_weight", "rbar2_gate", "angle_frame", "output_bin",
        "finalize_mask_mode", "histogram_edge_mode", "gradient_sample_mode",  # obsolete
    ]),
    "execution" => Set(["ny_chunks"]),
)

function _flatten_toml_config(raw, abspath_path::AbstractString)
    flat = Dict{String, Any}()
    for (table_name, table_value) in raw
        haskey(VGT_CONFIG_TABLE_KEYS, table_name) ||
            throw(ArgumentError("Unknown VGT config table `$table_name` in $abspath_path"))
        table_value isa AbstractDict ||
            throw(ArgumentError("Top-level VGT config entry `$table_name` must be a TOML table"))
        allowed = VGT_CONFIG_TABLE_KEYS[table_name]
        for (key, value) in table_value
            key in allowed ||
                throw(ArgumentError("Unknown VGT config key `$table_name.$key` in $abspath_path"))
            haskey(flat, key) &&
                throw(ArgumentError("Duplicate VGT config key `$key` in $abspath_path"))
            flat[key] = value
        end
    end
    return flat
end

function _config_key(cfg::AbstractDict, key::AbstractString, T::Type, abspath_path::AbstractString)
    haskey(cfg, key) || throw(ArgumentError("Missing required key '$key' in $abspath_path"))
    val = cfg[key]
    try
        return T === Any ? val : T(val)
    catch
        throw(ArgumentError("Key '$key' in $abspath_path: cannot convert $(repr(val)) to $T"))
    end
end

function _config_key_or(cfg::AbstractDict, key::AbstractString, default)
    return get(cfg, key, default)
end

function load_config(path::AbstractString="config.toml")
    abspath_path = abspath(path)
    isfile(abspath_path) || throw(ArgumentError("Config file not found: $abspath_path"))
    cfg = _flatten_toml_config(TOML.parsefile(abspath_path), abspath_path)
    noise = _config_key(cfg, "noise", Float64, abspath_path)
    noise_mode = _config_key_or(cfg, "noise_mode", :scalar)
    noise_tile_size = _config_key_or(cfg, "noise_tile_size", 128)
    noise_sigma_multiplier = _config_key_or(cfg, "noise_sigma_multiplier", 1.0)

    return validate_static_config(VGTConfig(
        input_file=_config_key(cfg, "input_file", String, abspath_path),
        output_dir=_config_key(cfg, "output_dir", String, abspath_path),
        keep_intermediate_chunks=_config_key_or(cfg, "keep_intermediate_chunks", true),
        noise=noise,
        dn=_config_key(cfg, "dn", Int, abspath_path),
        ker=_config_key(cfg, "ker", Int, abspath_path),
        ny_chunks=_config_key(cfg, "ny_chunks", Any, abspath_path),
        min_valid=_config_key(cfg, "min_valid", Int, abspath_path),
        fitter_mode=_config_key_or(cfg, "fitter_mode", :axial_vonmises),
        rbar2_weight=_config_key_or(cfg, "rbar2_weight", false),
        rbar2_gate=Float64(_config_key_or(cfg, "rbar2_gate", 0.0)),
        noise_mode=noise_mode,
        noise_tile_size=Int(noise_tile_size),
        noise_sigma_multiplier=Float64(noise_sigma_multiplier),
        zero_stokes_policy=_config_key_or(cfg, "zero_stokes_policy", :nan_before_smoothing),
        output_support_mode=_config_key_or(cfg, "output_support_mode", :smoothed_stokes_support),
        gradient_mode=_config_key_or(cfg, "gradient_mode", :central_checked),
        weight_file=String(_config_key_or(cfg, "weight_file", "")),
        angle_frame=_config_key_or(cfg, "angle_frame", :pixel),
        output_bin=_config_key_or(cfg, "output_bin", 1),
        velocity_range_kms=_config_key_or(cfg, "velocity_range_kms", Float64[]),
        # Obsolete keys: the constructor accepts the old default with a
        # warning and rejects anything else.
        finalize_mask_mode=get(cfg, "finalize_mask_mode", nothing),
        histogram_edge_mode=get(cfg, "histogram_edge_mode", nothing),
        gradient_sample_mode=get(cfg, "gradient_sample_mode", nothing),
        noise_spatial_exponent=get(cfg, "noise_spatial_exponent", nothing),
        noise_channel_exponent=get(cfg, "noise_channel_exponent", nothing),
    ))
end

"""
    selected_channels(cfg, header, nv) -> UnitRange{Int}

Channels (1-based) of the input cube used by the run: all of them, or those whose
centres lie inside `cfg.velocity_range_kms` (inclusive, 1e-6 km/s tolerance).
The spectral axis comes from CRVAL3/CRPIX3/CDELT3 (`velocity_axis_kms`).
"""
function selected_channels(cfg::VGTConfig, header, nv::Int)
    isempty(cfg.velocity_range_kms) && return 1:nv
    lo, hi = cfg.velocity_range_kms
    v, _ = velocity_axis_kms(header, nv)
    inside = findall(x -> lo - 1e-6 <= x <= hi + 1e-6, v)
    isempty(inside) && throw(ArgumentError(
        "velocity_range_kms = [$lo, $hi] contains no channel centre; the cube spans $(minimum(v))..$(maximum(v)) km/s"))
    return first(inside):last(inside)   # contiguous: the axis is linear
end

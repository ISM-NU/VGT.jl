# Configuration and execution settings for a VGT run.
#
# Most fields are not generic software options; they encode scientific choices
# such as window size, noise model, fitter mode, and masking policy.

const ANGLE_FRAMES = (:pixel, :galactic, :icrs)
# :isotropic_checked is opt-in (docs/validation-angle-calibration.md); :central_checked is the default.
const GRADIENT_MODES = (:central_checked, :sobel_checked, :isotropic_checked, :nan_edge_backward)
const FITTER_MODES = (:axial_vonmises,)
const NOISE_MODES = (:scalar, :channel, :spatial, :separable, :auto)
const ZERO_STOKES_POLICIES = (:nan_before_smoothing, :keep_zero_before_smoothing)
const OUTPUT_SUPPORT_MODES = (:smoothed_stokes_support, :direct_fit_support)

# Adding a field takes two edits: the struct field and its CONFIG_SPEC line below.
# Range checks go in _validate_config. Keyword, positional and TOML construction follow.
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
    # Output grid: 1 keeps the input grid; b > 1 averages the smoothed Q/U/W over
    # b x b input-pixel blocks before psi and reliability are computed, and
    # scales the celestial WCS to match (see `bin_products`).
    output_bin::Int
    # Velocity window [lo, hi] in km/s: only channels whose centres lie inside it
    # (inclusive) are used, as if the cube were cut to them. Empty = all channels.
    velocity_range_kms::Vector{Float64}
    # Optional whitening of the gradient components (off by default). With
    # gradient_whiten every gradient is mapped by C^(-1/2), C = gradient_noise_cov
    # = [sxx, syy, sxy] scaled to unit determinant, so anisotropic noise gradients
    # have no preferred angle; each window's fitted angle is mapped back with
    # C^(1/2). C is the covariance of one gradient (data units per pixel) in
    # emission-free channels with the stencil of gradient_mode (see
    # `measure_gradient_noise`).
    gradient_whiten::Bool
    gradient_noise_cov::Vector{Float64}
    # Alternative to gradient_noise_cov for noise whose anisotropy varies across
    # the map: FITS cube (nx, ny, 3) of (sxx, syy, sxy) on the input's sky grid,
    # used per pixel (see `gradient_noise_cov_map`). Exclusive with the constant.
    gradient_noise_cov_file::String
    # Opt-in: NaN Q/U where the window footprint reaches a non-finite pixel or the
    # array border (see `_edge_footprint_mask`). Off by default.
    edge_mask::Bool

    # Takes every field, in field order. Every construction path ends here, so
    # normalizing and validating here means invalid science settings fail
    # before a large FITS cube run starts, whichever constructor was used.
    function VGTConfig(values::Tuple)
        fields = map(_normalize_field, fieldnames(VGTConfig), values)
        _validate_config(NamedTuple{fieldnames(VGTConfig)}(fields))
        return new(fields...)
    end
end

# One entry per VGTConfig field, in field order: the TOML table the key lives
# in, its default (none = required), whether a TOML file must set it even when
# it has a default, and an optional normalizer. Without a normalizer the value
# is converted to the field's type.
const CONFIG_SPEC = (
    input_file = (table="input",),
    output_dir = (table="input", default="subcube", toml_required=true),
    keep_intermediate_chunks = (table="input", default=true),
    noise = (table="noise",),
    dn = (table="vgt",),
    ker = (table="vgt",),
    ny_chunks = (table="execution", default=1, toml_required=true,
                 normalize=v -> normalize_execution_setting(v, :ny_chunks)),
    min_valid = (table="vgt", default=100, toml_required=true),
    fitter_mode = (table="vgt", default=:axial_vonmises,
                   normalize=v -> _normalize_mode(v, :fitter_mode, FITTER_MODES)),
    rbar2_weight = (table="vgt", default=false),
    noise_mode = (table="noise", default=:scalar,
                  normalize=v -> _normalize_mode(v, :noise_mode, NOISE_MODES)),
    noise_tile_size = (table="noise", default=128),
    zero_stokes_policy = (table="vgt", default=:nan_before_smoothing,
                          normalize=v -> _normalize_mode(v, :zero_stokes_policy, ZERO_STOKES_POLICIES)),
    output_support_mode = (table="vgt", default=:smoothed_stokes_support,
                           normalize=v -> _normalize_mode(v, :output_support_mode, OUTPUT_SUPPORT_MODES)),
    gradient_mode = (table="vgt", default=:central_checked,
                     normalize=v -> _normalize_mode(v, :gradient_mode, GRADIENT_MODES)),
    noise_sigma_multiplier = (table="noise", default=1.0),
    rbar2_gate = (table="vgt", default=0.0),
    weight_file = (table="input", default=""),
    angle_frame = (table="vgt", default=:pixel,
                   normalize=v -> _normalize_mode(v, :angle_frame, ANGLE_FRAMES)),
    output_bin = (table="vgt", default=1),
    velocity_range_kms = (table="input", default=Float64[]),
    gradient_whiten = (table="noise", default=false),
    gradient_noise_cov = (table="noise", default=Float64[]),
    gradient_noise_cov_file = (table="noise", default=""),
    edge_mask = (table="vgt", default=false),
)
keys(CONFIG_SPEC) == fieldnames(VGTConfig) ||
    error("CONFIG_SPEC must list every VGTConfig field, in field order")

# Config keys removed in 0.1.0 and the only value each still tolerates (its old
# default, which changes nothing). Any other value would have changed the science.
const OBSOLETE_CONFIG_KEYS = (
    finalize_mask_mode = (table="vgt", old_default=:strict_nonpositive),
    histogram_edge_mode = (table="vgt", old_default=:modern_range),
    gradient_sample_mode = (table="vgt", old_default=:paired_finite),
    noise_spatial_exponent = (table="noise", old_default=1.0),
    noise_channel_exponent = (table="noise", old_default=1.0),
)

_coerce(::Type{Vector{Float64}}, value) = collect(Float64, value)  # Float64.(x) of an empty Vector{Any} stays Vector{Any}
_coerce(T::Type, value) = T(value)

function _normalize_field(name::Symbol, value)
    spec = CONFIG_SPEC[name]
    haskey(spec, :normalize) && return spec.normalize(value)
    return _coerce(fieldtype(VGTConfig, name), value)
end

# VGTConfig from fields given by name; the rest take their defaults.
function _build_config(given::AbstractDict{Symbol})
    unknown = setdiff(keys(given), fieldnames(VGTConfig))
    isempty(unknown) ||
        throw(ArgumentError("Unknown VGTConfig field(s): $(join(sort!(collect(unknown)), ", "))"))
    values = map(fieldnames(VGTConfig)) do name
        spec = CONFIG_SPEC[name]
        haskey(given, name) && return given[name]
        haskey(spec, :default) && return spec.default
        throw(ArgumentError("Missing required keyword `$name`"))
    end
    return VGTConfig(values)
end

# Positional form: fields in declaration order, the first six required.
function VGTConfig(input_file::AbstractString, output_dir::AbstractString,
                   keep_intermediate_chunks::Bool, noise, dn, ker, rest...)
    args = (input_file, output_dir, keep_intermediate_chunks, noise, dn, ker, rest...)
    length(args) <= fieldcount(VGTConfig) ||
        throw(ArgumentError("VGTConfig takes at most $(fieldcount(VGTConfig)) positional arguments; got $(length(args))"))
    return _build_config(Dict{Symbol, Any}(zip(fieldnames(VGTConfig), args)))
end

function VGTConfig(; n_threads=nothing, kwargs...)
    given = Dict{Symbol, Any}(kwargs)
    _check_obsolete_keys(given)
    n_threads === nothing ||
        @warn "`n_threads` no longer controls execution; start Julia with `--threads=N` or set JULIA_NUM_THREADS" maxlog=1
    foreach(key -> delete!(given, key), keys(OBSOLETE_CONFIG_KEYS))
    return _build_config(given)
end

# Copy of `cfg` with some fields replaced. Copies every field by name, so a new
# field cannot be dropped by a hand-written field list.
function config_with(cfg::VGTConfig; changes...)
    return _build_config(Dict{Symbol, Any}(
        f => (haskey(changes, f) ? changes[f] : getfield(cfg, f)) for f in fieldnames(VGTConfig)))
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

const REMOVED_IN_0_1_0 = "removed in VGT 0.1.0 (public release); it remains in the archived private repository (see docs/provenance.md)"

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
    choices = length(allowed) == 2 ? join(repr.(allowed), " or ") : join(repr.(allowed), ", ", ", or ")
    throw(ArgumentError("$(name) must be one of $choices; got $(repr(value))"))
end

function _check_obsolete_keys(given::AbstractDict{Symbol})
    for (key, spec) in pairs(OBSOLETE_CONFIG_KEYS)
        value = get(given, key, nothing)
        value === nothing && continue
        seen = spec.old_default isa Symbol ? _obsolete_symbol(value) : value
        if seen == spec.old_default
            @warn "Config key `$key` is obsolete and ignored (value $(repr(value)) is the old default)"
        else
            throw(ArgumentError("$key=$(repr(value)) was $REMOVED_IN_0_1_0"))
        end
    end
    return nothing
end

_obsolete_symbol(value) = value isa Union{Symbol, AbstractString} ? _symbol_config_value(value) : value

# Settings the kernel also takes as loose keywords (VGT_QU_mv_2d is a public
# entry point), so the kernel and the config share these checks.
function _check_kernel_settings(dn, min_valid, rbar2_gate)
    iseven(dn) || throw(ArgumentError("dn must be even; got $dn"))
    min_valid >= 1 || throw(ArgumentError("min_valid must be >= 1"))
    min_valid <= dn * dn ||
        throw(ArgumentError("min_valid=$min_valid exceeds the maximum window area dn^2=$(dn * dn)"))
    0.0 <= rbar2_gate <= 1.0 ||
        throw(ArgumentError("rbar2_gate must be in [0, 1], got $rbar2_gate"))
    return nothing
end

# Checks the normalized fields (a NamedTuple in field order) at construction.
function _validate_config(c)
    c.noise > 0 || throw(ArgumentError("noise must be > 0"))
    c.dn >= 2 || throw(ArgumentError("dn must be >= 2"))
    c.ker >= 1 || throw(ArgumentError("ker must be >= 1"))
    _check_kernel_settings(c.dn, c.min_valid, c.rbar2_gate)
    c.noise_tile_size >= 1 || throw(ArgumentError("noise_tile_size must be >= 1"))
    c.noise_sigma_multiplier > 0 ||
        throw(ArgumentError("noise_sigma_multiplier must be > 0"))
    c.output_bin >= 1 || throw(ArgumentError("output_bin must be a positive integer; got $(c.output_bin)"))
    isempty(c.velocity_range_kms) || (length(c.velocity_range_kms) == 2 && all(isfinite, c.velocity_range_kms) &&
                                      c.velocity_range_kms[1] < c.velocity_range_kms[2]) ||
        throw(ArgumentError("velocity_range_kms must be empty or [lo, hi] with lo < hi; got $(c.velocity_range_kms)"))
    _validate_gradient_whiten(c.gradient_whiten, c.gradient_noise_cov, c.gradient_noise_cov_file, c.gradient_mode)
    return nothing
end

function _validate_gradient_whiten(gradient_whiten, gradient_noise_cov, gradient_noise_cov_file, gradient_mode)
    gradient_whiten || return nothing
    has_const = !isempty(gradient_noise_cov)
    has_file = !isempty(gradient_noise_cov_file)
    (has_const && has_file) &&
        throw(ArgumentError("set either gradient_noise_cov or gradient_noise_cov_file, not both"))
    (has_const || has_file) ||
        throw(ArgumentError("gradient_whiten needs gradient_noise_cov = [sxx, syy, sxy] or gradient_noise_cov_file (covariance of one gradient in emission-free channels; see scripts/analyze_inputs.jl --noise-vrange)"))
    if has_const
        length(gradient_noise_cov) == 3 ||
            throw(ArgumentError("gradient_noise_cov must be [sxx, syy, sxy]; got $(gradient_noise_cov)"))
        sxx, syy, sxy = gradient_noise_cov
        (all(isfinite, gradient_noise_cov) && sxx > 0 && syy > 0 && sxx * syy - sxy^2 > 0) ||
            throw(ArgumentError("gradient_noise_cov = [sxx, syy, sxy] must be finite and positive definite (sxx > 0, syy > 0, sxx*syy > sxy^2); got $(gradient_noise_cov)"))
    end
    gradient_mode === :nan_edge_backward &&
        throw(ArgumentError("gradient_whiten is not supported with gradient_mode=:nan_edge_backward; use :central_checked, :sobel_checked or :isotropic_checked"))
    return nothing
end

function chunk_padding_pixels(cfg::VGTConfig)
    # Rows a chunk reads beyond its own: a window centred on row j covers rows
    # j - dn/2 + 1 .. j + dn/2, and the central and Sobel stencils need one more
    # row on each side of a gradient (the backward difference only the row
    # before, which the window's first row already has). One row less made the
    # last window rows of every chunk NaN, so results depended on ny_chunks
    # (test/chunk_seams.jl). Final smoothing is applied only after stitching.
    return div(cfg.dn, 2) + (cfg.gradient_mode === :nan_edge_backward ? 0 : 1)
end

function choose_auto_chunks(ny::Int, dn::Int, n_threads::Int)
    target = min(max(1, 2 * n_threads), 32)
    # dn is the correctness floor (process_chunk pads by chunk_padding_pixels on each side);
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
    n_threads = max(1, cpu_threads)
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
    cfg.dn <= nx || throw(ArgumentError("dn=$(cfg.dn) is too large for nx=$nx (need nx ≥ $(cfg.dn))"))
    cfg.dn <= ny || throw(ArgumentError("dn=$(cfg.dn) is too large for ny=$ny (need ny ≥ $(cfg.dn))"))
    cfg.output_bin <= min(nx, ny) ||
        throw(ArgumentError("output_bin=$(cfg.output_bin) is larger than the map ($nx x $ny pixels)"))

    if cfg.gradient_whiten && !isempty(cfg.gradient_noise_cov_file)
        isfile(cfg.gradient_noise_cov_file) ||
            throw(ArgumentError("gradient_noise_cov_file not found: $(cfg.gradient_noise_cov_file)"))
        size(read_gradient_noise_cov_map(cfg.gradient_noise_cov_file)) == (nx, ny, 3) ||
            throw(ArgumentError("gradient_noise_cov_file must be a ($nx, $ny, 3) cube (sxx, syy, sxy) on the input grid"))
    end

    cpu_threads = max(1, Threads.nthreads())
    budget = memory_budget_bytes()
    auto_execution = choose_auto_execution(cfg, nx, ny, nv, cpu_threads, budget)
    ny_chunks = auto_execution.ny_chunks

    return config_with(cfg; ny_chunks=ny_chunks)
end

# TOML table -> allowed keys, from the field table and the obsolete keys.
const VGT_CONFIG_TABLE_KEYS = let tables = Dict{String, Set{String}}()
    for specs in (CONFIG_SPEC, OBSOLETE_CONFIG_KEYS), (name, spec) in pairs(specs)
        push!(get!(Set{String}, tables, spec.table), String(name))
    end
    tables
end

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

function load_config(path::AbstractString="config.toml")
    abspath_path = abspath(path)
    isfile(abspath_path) || throw(ArgumentError("Config file not found: $abspath_path"))
    raw = _flatten_toml_config(TOML.parsefile(abspath_path), abspath_path)
    given = Dict{Symbol, Any}(Symbol(key) => value for (key, value) in raw)
    for (name, spec) in pairs(CONFIG_SPEC)
        (haskey(spec, :default) && !get(spec, :toml_required, false)) && continue
        T = fieldtype(VGTConfig, name)
        given[name] = _config_key(raw, String(name), T isa Union ? Any : T, abspath_path)
    end
    # Obsolete keys go through unchanged: the constructor accepts the old
    # default with a warning and rejects anything else.
    return VGTConfig(; given...)
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

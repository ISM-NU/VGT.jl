function _fitsio_read_header(path::AbstractString)
    f = FITS(path)
    try
        return read_header(f[_fitsio_image_hdu_index(f)])
    finally
        close(f)
    end
end

function _header_bool(value)
    value isa Bool && return value
    lowercase(strip(string(value))) in ("t", "true", "1")
end

function _header_int(header, key::AbstractString, default::Int=-1)
    haskey(header, key) || return default
    try
        return Int(header[key])
    catch
        return default
    end
end

function _is_compressed_image_header(header)
    haskey(header, "ZIMAGE") && _header_bool(header["ZIMAGE"])
end

function _is_plain_image_header(header)
    xtension = haskey(header, "XTENSION") ? uppercase(strip(string(header["XTENSION"]))) : ""
    return xtension in ("", "IMAGE") && _header_int(header, "NAXIS") > 0
end

function _fitsio_image_hdu_index(f::FITS)
    plain_idx = nothing
    for i in 1:length(f)
        header = read_header(f[i])
        _is_compressed_image_header(header) && return i
        if plain_idx === nothing && _is_plain_image_header(header)
            plain_idx = i
        end
    end
    plain_idx !== nothing && return plain_idx
    return 1
end

function _plain_header(header)
    out = Dict{String, Any}()
    for card in header
        if hasproperty(card, :name) && hasproperty(card, :value)
            out[String(card.name)] = card.value
        end
    end
    return out
end

function read_fits_header(path::AbstractString)
    try
        header = _fitsio_read_header(path)
        return header
    catch err
        @warn "FITSIO header read failed; falling back to AstroFITS" path=path exception=(err, catch_backtrace())
        return _plain_header(read(path, AstroFITS.FitsHeader))
    end
end

function _fitsio_read_image(path::AbstractString, ranges...)
    f = FITS(path)
    try
        hdu = f[_fitsio_image_hdu_index(f)]
        return isempty(ranges) ? read(hdu) : read(hdu, ranges...)
    finally
        close(f)
    end
end

function read_fits_image(path::AbstractString)
    _is_fz_path(path) && return _fitsio_read_image(path)
    try
        return AstroFITS.readfits(path)
    catch err
        @warn "AstroFITS image read failed; falling back to FITSIO" path=path exception=(err, catch_backtrace())
        return _fitsio_read_image(path)
    end
end

function read_fits_image(path::AbstractString, ranges...)
    return _fitsio_read_image(path, ranges...)
end

# (key, value, comment) for every card, by position: repeated COMMENT/HISTORY
# cards must not collapse into copies of the first one, as key lookup does.
function _header_cards(header)
    if header isa FITSHeader
        return [(header.keys[i], header.values[i], header.comments[i]) for i in eachindex(header.keys)]
    end
    return [(String(k), header[k], "") for k in keys(header)]
end

function _header_pairs(header)
    header === nothing && return nothing
    pairs = Pair{String, Any}[]
    for (key, value, comment) in _header_cards(header)
        keystr = String(key)
        keystr in ("SIMPLE", "BITPIX", "NAXIS", "EXTEND") && continue
        startswith(keystr, "NAXIS") && continue
        if keystr in ("COMMENT", "HISTORY", "")
            # FITSIO keeps commentary text in the comment slot.
            push!(pairs, keystr => (value === nothing ? comment : value))
        elseif isempty(comment)
            push!(pairs, keystr => value)
        else
            push!(pairs, keystr => (value, comment))
        end
    end
    return pairs
end

function _is_fz_path(path::AbstractString)
    lowercase(path) == ".fz" && return true
    return endswith(lowercase(path), ".fz")
end

function normalize_fits_extension(ext::AbstractString)
    normalized = lowercase(strip(String(ext)))
    isempty(normalized) && return ".fits"
    startswith(normalized, ".") || (normalized = ".$normalized")
    normalized in (".fits", ".fit") && return normalized
    normalized in (".fz", ".fits.fz", ".fit.fz") && return ".fits.fz"
    throw(ArgumentError("FITS extension must be .fits, .fit, .fits.fz, or .fz; got $(repr(ext))"))
end

function _write_plain_fits(path::AbstractString, data; header=nothing)
    header_pairs = header === nothing ? Pair{String, Any}[] : _header_pairs(header)
    try
        AstroFITS.writefits(path, header_pairs, data; overwrite=true)
    catch err
        @warn "AstroFITS write failed; falling back to FITSIO" path=path exception=(err, catch_backtrace())
        rm(path; force=true)
        f = FITS(path, "w")
        try
            header === nothing ? write(f, data) : write(f, data, header=header)
        finally
            close(f)
        end
    end
    return path
end

function write_fits(path::AbstractString, data; header=nothing)
    mkpath(dirname(path))
    if _is_fz_path(path)
        Sys.which("fpack") === nothing &&
            throw(ArgumentError("Writing .fz FITS output requires fpack on PATH"))
        mktempdir(dirname(path)) do dir
            plain_path = joinpath(dir, replace(basename(path), r"\.fz$"i => ""))
            _write_plain_fits(plain_path, data; header=header)
            run(`fpack -g -Y -O $path $plain_path`)
        end
    else
        _write_plain_fits(path, data; header=header)
    end
    return path
end

# Spectral axis in km/s from a cube header (CRVAL3/CRPIX3/CDELT3, CUNIT3 or value scale).
function velocity_axis_kms(header, nv::Int)
    crval = Float64(header["CRVAL3"])
    crpix = Float64(header["CRPIX3"])
    cdelt = Float64(header["CDELT3"])
    # FITS convention: channel k (1-based) is at CRVAL3 + (k - CRPIX3) * CDELT3.
    native = [crval + (k - crpix) * cdelt for k in 1:nv]
    unit = try
        lowercase(strip(String(header["CUNIT3"])))
    catch
        ""
    end
    # Some survey cubes omit or vary CUNIT3 spelling; fall back to value scale
    # so native m/s axes are still converted before physical products are made.
    scale = if unit in ("m/s", "ms-1", "m s-1", "m/s ")
        1e-3
    elseif unit in ("km/s", "kms-1", "km s-1")
        1.0
    else
        maximum(abs.(native)) > 1_000 ? 1e-3 : 1.0
    end
    return native .* scale, abs(cdelt * scale)
end

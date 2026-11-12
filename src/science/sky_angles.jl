# Map headers and sky-frame angles for 2D VGT products.
#
# VGT angles (`psi`) are measured counter-clockwise from +FITS axis 1 in the
# map's own pixel grid; that is the frame in which gradients and Q/U smoothing
# are computed. Comparisons with Planck or starlight need an IAU position angle
# (east of north) in a named celestial frame, and on LMC-sized maps the local
# north direction moves by many degrees across the grid, so the conversion is
# done per pixel through the WCS.

# Cards describing a third (spectral) axis or the input data scaling. A 2D map
# that keeps them describes the input cube, not itself.
const _NON_CELESTIAL_KEY_PATTERNS = (
    r"^NAXIS\d*$",
    r"^(CTYPE|CRVAL|CRPIX|CDELT|CUNIT|CROTA|CNAME|CRDER|CSYER)[3-9]$",
    r"^(PC|CD)0*([3-9]_0*\d+|\d+_0*[3-9])$",
    r"^(PV|PS)[3-9]_\d+$",
    r"^(SPECSYS|SSYSOBS|SSYSSRC|RESTFRQ|RESTFREQ|RESTWAV|VELREF|VELO-LSR|VLSR|ZSOURCE|ALTRVAL|ALTRPIX)$",
    r"^(BUNIT|BSCALE|BZERO|BLANK|DATAMIN|DATAMAX|CHECKSUM|DATASUM|WCSAXES)$",
    r"^(SIMPLE|BITPIX|EXTEND)$",
    # The input beam is not the product's resolution; it is re-recorded as INBMAJ etc.
    r"^(BMAJ|BMIN|BPA)$",
)

function _keep_celestial_card(key::AbstractString)
    return !any(p -> occursin(p, key), _NON_CELESTIAL_KEY_PATTERNS)
end

"""
    celestial_header(header) -> FITSHeader

Header of a 2D map on the same sky grid as `header` (a 2D or 3D cube header):
celestial WCS and provenance cards are kept, spectral-axis, data-scaling and
beam cards are dropped, and `WCSAXES = 2` is set. The input beam, if present,
is kept as `INBMAJ`/`INBMIN`/`INBPA`.
"""
function celestial_header(header)
    keys_out = String[]
    values_out = Any[]
    comments_out = String[]
    push_card!(k, v, c) = (push!(keys_out, k); push!(values_out, v); push!(comments_out, c))
    push_card!("WCSAXES", 2, "number of celestial WCS axes")
    for (key, value, comment) in _header_cards(header)
        k = String(key)
        _keep_celestial_card(k) || continue
        push_card!(k, value, comment)
    end
    for k in ("BMAJ", "BMIN", "BPA")
        haskey(header, k) || continue
        push_card!("IN" * k, header[k], "input cube $(k), deg (not this map's resolution)")
    end
    return FITSHeader(keys_out, values_out, comments_out)
end

function _celestial_wcs(header)
    transforms = WCS.from_header(_wcs_header_string(celestial_header(header)); ignore_rejected=true)
    isempty(transforms) && throw(ArgumentError("Header has no celestial WCS"))
    wcs = transforms[1]
    wcs.naxis == 2 || throw(ArgumentError("Expected a 2-axis celestial WCS, got $(wcs.naxis) axes"))
    return wcs
end

function _map_frame(header)
    ctype1 = uppercase(String(header["CTYPE1"]))
    startswith(ctype1, "RA--") && return :icrs
    startswith(ctype1, "GLON") && return :galactic
    throw(ArgumentError("Unsupported celestial frame CTYPE1=$(ctype1); expected RA--- or GLON-"))
end

# ICRS -> Galactic rotation (Hipparcos, ESA 1997, Vol. 1 Sec. 1.5.3); rows are
# the Galactic x, y, z axes expressed in ICRS. FK5 J2000 differs by < 0.1",
# far below what matters for orientation angles.
const _ICRS_TO_GALACTIC = [
    -0.0548755604162154 -0.8734370902348850 -0.4838350155487132
     0.4941094278755837 -0.4448296299600112  0.7469822444972189
    -0.8676661490190047 -0.1980763734312015  0.4559837761750669
]

_unit_vector(lon, lat) = (cosd(lat) * cosd(lon), cosd(lat) * sind(lon), sind(lat))

function _lonlat(v)
    x, y, z = v
    return mod(atand(y, x), 360.0), atand(z, hypot(x, y))
end

function _convert_frame(lon, lat, from::Symbol, to::Symbol)
    from === to && return lon, lat
    v = collect(_unit_vector(lon, lat))
    if from === :icrs && to === :galactic
        return _lonlat(_ICRS_TO_GALACTIC * v)
    elseif from === :galactic && to === :icrs
        return _lonlat(transpose(_ICRS_TO_GALACTIC) * v)
    end
    throw(ArgumentError("Unsupported frame conversion $(from) -> $(to)"))
end

# Small sky step for the local north/east pixel vectors. Projection curvature
# over 1e-4 deg is negligible for any survey pixel scale used here.
const _SKY_STEP_DEG = 1e-4

"""
    pixel_to_position_angle(psi, header; frame=:galactic) -> Matrix{Float64}

Convert axial angles `psi` (radians, counter-clockwise from +FITS axis 1 in the
map's pixel grid, as written by the VGT pipeline) into IAU position angles
(radians, east of north) in `frame` (`:galactic` or `:icrs`), wrapped to
[-pi/2, pi/2). At every pixel the local north and east steps are projected
through the WCS and `psi`'s direction is decomposed on them, so grid rotation,
axis parity and non-square pixels are all accounted for. NaN stays NaN.
"""
function pixel_to_position_angle(psi::AbstractMatrix, header; frame::Symbol=:galactic)
    frame in (:galactic, :icrs) || throw(ArgumentError("frame must be :galactic or :icrs, got $(frame)"))
    wcs = _celestial_wcs(header)
    native = _map_frame(header)
    nx, ny = size(psi)
    out = fill(NaN, nx, ny)
    finite = findall(isfinite, psi)
    isempty(finite) && return out
    n = length(finite)
    # FITS pixel (i, j) is Julia array element [i, j] (both 1-based).
    pix = Matrix{Float64}(undef, 2, n)
    for (c, idx) in enumerate(finite)
        pix[1, c] = idx[1]
        pix[2, c] = idx[2]
    end
    world = pix_to_world(wcs, pix)
    north = Matrix{Float64}(undef, 2, n)
    east = Matrix{Float64}(undef, 2, n)
    for c in 1:n
        lon, lat = _convert_frame(world[1, c], world[2, c], native, frame)
        north[:, c] .= _convert_frame(lon, lat + _SKY_STEP_DEG, frame, native)
        east[:, c] .= _convert_frame(lon + _SKY_STEP_DEG / cosd(lat), lat, frame, native)
    end
    pix_n = world_to_pix(wcs, north)
    pix_e = world_to_pix(wcs, east)
    for (c, idx) in enumerate(finite)
        nxv, nyv = pix_n[1, c] - pix[1, c], pix_n[2, c] - pix[2, c]
        exv, eyv = pix_e[1, c] - pix[1, c], pix_e[2, c] - pix[2, c]
        dx, dy = cos(psi[idx]), sin(psi[idx])
        # Solve (dx, dy) = a * north + b * east; position angle = atan(b, a).
        det = nxv * eyv - exv * nyv
        a = (dx * eyv - exv * dy) / det
        b = (nxv * dy - dx * nyv) / det
        chi = atan(b, a)
        out[idx] = mod(chi + pi / 2, pi) - pi / 2
    end
    return out
end

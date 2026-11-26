module VGT

using AstroFITS
using CairoMakie
using Dates
using DrWatson
using FFTW
using FITSIO
using Images
using JLD2
using LinearAlgebra: mul!
using Random
using SHA
using Statistics
using TOML
using WCS

const FIT_LOW_CONFIDENCE_RBAR2_THRESHOLD = 0.5

include("config.jl")
include("paths.jl")
include("fits_io.jl")
include("regrid.jl")
include("science/noise.jl")
include("science/fitters.jl")
include("science/vgt.jl")
include("science/products.jl")
include("science/sky_angles.jl")
include("pipeline/chunks.jl")
include("pipeline/run.jl")
include("pipeline/outputs.jl")

export VGTConfig, load_config, run_pipeline, write_pipeline_outputs
export velocity_axis_kms
export block_mean, VGT_QU_mv_2d, VGT_QU_mv_2d_profiled
export measure_gradient_noise, gradient_anisotropy
export finalize_vgt_products
export write_fits, read_fits_header, read_fits_image
export celestial_header, pixel_to_position_angle
export experiment_default_dir

end

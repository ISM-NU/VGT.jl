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

export VGTConfig, load_config, run_pipeline, write_pipeline_outputs, is_interrupt_exception
export read_cube_metadata, getcoor, velocity_axis_kms
export WindowFitResult, block_mean, fit_window_result, VGT_QU_mv_2d, VGT_QU_mv_2d_profiled
export robust_negative_noise, build_noise_cube
export measure_gradient_noise, gradient_anisotropy, gradient_noise_cov_map
export write_gradient_noise_cov_map
export imfilter_gaussian, compute_stokes_intensity, finalize_vgt_products, headless_angle_difference
export write_fits, normalize_fits_extension, read_fits_header, read_fits_image
export celestial_header, pixel_to_position_angle
export source_axis_position, linear_bracket, bilinear_sample, velocity_overlap_indices
export rawdatadir, experiment_default_dir, parameterized_output_dir, experiment_summary_path
export collect_experiment_summaries

end

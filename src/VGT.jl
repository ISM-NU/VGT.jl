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

export VGTConfig, VGTDiagnostics
export load_config, validate_static_config, resolve_execution_config, chunk_padding_pixels
export read_cube_metadata, getcoor
export WindowFitResult, fit_window_result
export fit_axial_vonmises, fit_window_angle
export robust_negative_noise, estimate_channel_noise, estimate_spatial_noise
export NoiseModel, separable_noise_cube, build_noise_model
export measure_gradient_noise, gradient_anisotropy, gradient_noise_cov_map
export read_gradient_noise_cov_map, write_gradient_noise_cov_map
export noise_for_chunk, build_noise_cube
export velocity_axis_kms, selected_channels
export VGT_QU_mv_2d, VGT_QU_mv_2d_profiled
export process_chunk, chunk_y_bounds, save_chunk_atomic
export NEGATIVE_MAD_NORMAL
export choose_auto_threads, choose_auto_chunks, choose_auto_execution
export estimate_peak_chunk_bytes

export PipelineRun
export imfilter_gaussian
export compute_stokes_intensity, compute_reliability
export finalize_vgt_products, headless_angle_difference, block_mean, bin_header, bin_products
export coherence_phasor, coherence_floor2, debias_coherence, cross_coherence, split_coherence
export write_fits, normalize_fits_extension, next_available_output_path, prefixed_input_filename
export read_fits_header, read_fits_image
export celestial_header, pixel_to_position_angle, vgt_product_header
export source_axis_position, linear_bracket, bilinear_sample, velocity_overlap_indices
export project_config_path, rawdatadir, processeddatadir, simulationdir
export run_parameters, run_parameter_hash, run_savename, default_output_dir, experiment_default_dir
export parameterized_output_dir, experiment_output_dir, experiment_summary_dir, experiment_summary_path
export collect_experiment_summaries, resolved_output_dir
export run_metadata, write_run_metadata, write_pipeline_outputs, stitch_chunk_outputs
export cleanup_interrupted_run!
export is_interrupt_exception, run_pipeline

end

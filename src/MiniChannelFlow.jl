module MiniChannelFlow

using LinearAlgebra
using SparseArrays

include("fluids.jl")
include("geometry.jl")
include("correlations.jl")
include("rectduct.jl")
include("graetz.jl")
include("ns2d.jl")
include("heatsink.jl")

# fluids
export Fluid, water, WATER_20C, kinematic_viscosity, thermal_diffusivity, prandtl
# geometry
export RectChannel, hydraulic_diameter, aspect_ratio, cross_section_area,
       wetted_perimeter, classify_channel, reynolds_number
# correlations
export fRe_fanning, fRe_darcy, nusselt_T, nusselt_H1, apparent_fRe,
       hydrodynamic_entrance_length, thermal_entrance_length, pressure_drop
# rectangular duct
export exact_fRe, exact_velocity, DuctCrossSection, solve_duct,
       normalized_velocity, fRe_numeric, nusselt_H1_numeric
# thermal entrance
export ThermalEntrance, thermal_entrance, thermal_entrance_length_numeric
# 2-D Navier-Stokes
export PlaneChannelConfig, PlaneChannelSolution, solve_plane_channel,
       centerline_velocity, pressure_profile, fully_developed_fRe,
       entrance_length, incremental_pressure_drop
# heat sink
export HeatSink, HeatSinkResult, evaluate_heatsink, flow_rate_for_pressure_drop,
       optimize_channel_width, modify

end # module

using MiniChannelFlow
using Printf

function bench(f; n = 3)
    f()                                  # compile
    return minimum(@elapsed(f()) for _ in 1:n)
end

println("MiniChannelFlow.jl benchmarks (best of 3, seconds)")
@printf("  exact_fRe (α*=0.25, 100×100 terms)        %8.4f\n", bench(() -> exact_fRe(0.25)))
@printf("  solve_duct nz=64 (α*=1)                    %8.4f\n", bench(() -> solve_duct(1.0; nz = 64)))
@printf("  solve_duct nz=32 (α*=0.1, 320×32 cells)    %8.4f\n", bench(() -> solve_duct(0.1; nz = 32)))
@printf("  thermal_entrance α*=1 nz=32                %8.4f\n", bench(() -> thermal_entrance(1.0; nz = 32)))
@printf("  solve_plane_channel 120×16 Re=100          %8.4f\n",
        bench(() -> solve_plane_channel(PlaneChannelConfig(Re = 100.0, Nx = 120, Ny = 16, L = 12.0)); n = 1))
hs = HeatSink(channel_width = 50e-6, wall_width = 50e-6, channel_height = 302e-6,
              base_thickness = 100e-6, length = 0.01, width = 0.01, k_solid = 148.0)
@printf("  optimize_channel_width                     %8.4f\n",
        bench(() -> optimize_channel_width(hs; pressure_drop = 214e3, heat_flux = 7.9e6); n = 1))

using MiniChannelFlow
using Printf
using CairoMakie

mkpath(joinpath(@__DIR__, "output"))
hs = HeatSink(channel_width = 50e-6, wall_width = 50e-6, channel_height = 302e-6,
              base_thickness = 100e-6, length = 0.01, width = 0.01, k_solid = 148.0)
q = 7.9e6                                        # W/m² (790 W/cm²)

fig = Figure(size = (900, 400))
ax = Axis(fig[1, 1]; xscale = log10, xlabel = "channel width wc [µm]", ylabel = "R_total [K/W]",
          title = "Thermal resistance at fixed ΔP")
println("  ΔP [kPa]   wc_opt [µm]   R_opt [K/W]   Q [cm³/s]   Re")
for Δp in (50e3, 100e3, 214e3, 400e3)
    wcs = exp.(range(log(20e-6), log(300e-6); length = 40))
    Rs = Float64[]
    for wc in wcs
        h2 = modify(hs; channel_width = wc)
        Q = flow_rate_for_pressure_drop(h2, Δp; heat_flux = q)
        r = evaluate_heatsink(h2; flow_rate = Q, heat_flux = q)
        push!(Rs, r.laminar ? r.R_total : NaN)     # hide non-laminar designs
    end
    lines!(ax, wcs .* 1e6, Rs; label = "ΔP = $(Int(Δp / 1e3)) kPa")
    opt = optimize_channel_width(hs; pressure_drop = Δp, heat_flux = q)
    scatter!(ax, [opt.channel_width * 1e6], [opt.R_total]; color = :black)
    @printf("%9.0f   %11.1f   %11.4f   %9.2f   %5.0f\n", Δp / 1e3, opt.channel_width * 1e6,
            opt.R_total, opt.flow_rate * 1e6, opt.result.Re)
end
axislegend(ax; position = :rt)
save(joinpath(@__DIR__, "output", "04_heatsink_optimization.png"), fig)

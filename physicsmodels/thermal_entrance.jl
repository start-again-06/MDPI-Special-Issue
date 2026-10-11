using MiniChannelFlow
using Printf
using CairoMakie

mkpath(joinpath(@__DIR__, "output"))
αs = [1.0, 0.5, 0.25, 0.125]
ref = Dict(1.0 => 2.976, 0.5 => 3.391, 0.25 => 4.439, 0.125 => 5.597)   # Shah & London Nu_T

fig = Figure(size = (900, 420))
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = "x* = x/(Dh Re Pr)",
          ylabel = "local Nu_T", title = "Thermal entrance, constant wall temperature")
println("   α*   Nu_T(∞) solver   literature   error %   x*_t (2 %)")
open(joinpath(@__DIR__, "output", "02_graetz.csv"), "w") do io
    println(io, "alpha,xstar,Nu_local,Nu_mean")
    for α in αs
        r = thermal_entrance(α; nz = 32, x_end = 0.5)
        xt = thermal_entrance_length_numeric(r)
        @printf("%6.3f   %12.4f   %10.3f   %7.2f   %8.4f\n",
                α, r.Nu[end], ref[α], 100 * (r.Nu[end] / ref[α] - 1), xt)
        lines!(ax, r.xmid[r.xmid .> 1e-5], r.Nu[r.xmid .> 1e-5]; label = "α* = $α")
        hlines!(ax, [ref[α]]; linestyle = :dot, color = :gray)
        for (x, nl, nm) in zip(r.xmid[2:end], r.Nu[2:end], r.Nu_mean[2:end])
            println(io, α, ",", x, ",", nl, ",", nm)
        end
    end
end
axislegend(ax; position = :lb)
save(joinpath(@__DIR__, "output", "02_thermal_entrance.png"), fig)

# What does this mean for a real mini-channel?  Water, Re = 200, Dh = 1 mm
Pr = prandtl(water(25.0))
println("\nWater (Pr = $(round(Pr; digits = 2))), Re = 200, Dh = 1 mm:")
r = thermal_entrance(1.0; nz = 32, x_end = 0.5)
xt = thermal_entrance_length_numeric(r)
@printf("  square duct thermal entrance length ≈ %.0f mm (hydrodynamic: %.1f mm)\n",
        xt * 1e-3 * 200 * Pr * 1e3, hydrodynamic_entrance_length(200.0) * 1.0)

using MiniChannelFlow
using Printf
using CairoMakie

mkpath(joinpath(@__DIR__, "output"))
fig = Figure(size = (900, 400))
ax = Axis(fig[1, 1]; xlabel = "x / Dh", ylabel = "u_centre / U", title = "Centre-line velocity development")

println("   Re   L_h/Dh (solver)   L_h/Dh (Chen)   f·Re (x>0.6L)   iterations")
for Re in (50.0, 100.0, 200.0, 400.0)
    L = max(30.0, 0.08 * Re * 2)              # domain in units of H
    cfg = PlaneChannelConfig(Re = Re, Nx = round(Int, 10L), Ny = 24, L = L)
    sol = solve_plane_channel(cfg)
    ucl = centerline_velocity(sol)
    x = (0:cfg.Nx) .* sol.dx ./ 2             # x / Dh
    lines!(ax, collect(x), ucl; label = "Re = $(Int(Re))")
    @printf("%5.0f   %14.3f   %13.3f   %13.3f   %10d\n", Re, entrance_length(sol),
            hydrodynamic_entrance_length(Re), fully_developed_fRe(sol), sol.iterations)
end
axislegend(ax; position = :rb)

# K(Re) = K∞ + c/Re  — two-point extrapolation to Re → ∞
Ks = Dict{Float64,Float64}()
for Re in (200.0, 400.0)
    sol = solve_plane_channel(PlaneChannelConfig(Re = Re, Nx = 400, Ny = 32, L = 40.0))
    _, K = incremental_pressure_drop(sol)
    Ks[Re] = K[length(K)÷2]
end
c = (Ks[200.0] - Ks[400.0]) / (1 / 200 - 1 / 400)
Kinf = Ks[400.0] - c / 400
@printf("\nK(200) = %.3f, K(400) = %.3f  →  K∞ = %.3f (Shah & London: 0.674), c = %.1f\n",
        Ks[200.0], Ks[400.0], Kinf, c)

ax2 = Axis(fig[1, 2]; xlabel = "1/Re", ylabel = "K", title = "Incremental pressure drop K")
scatter!(ax2, [1 / 200, 1 / 400], [Ks[200.0], Ks[400.0]]; markersize = 12, label = "solver")
lines!(ax2, [0.0, 1 / 150], [Kinf, Kinf + c / 150]; linestyle = :dash, label = "linear fit")
hlines!(ax2, [0.674]; color = :gray, linestyle = :dot, label = "K∞ = 0.674")
axislegend(ax2; position = :lt)
save(joinpath(@__DIR__, "output", "03_developing_flow.png"), fig)

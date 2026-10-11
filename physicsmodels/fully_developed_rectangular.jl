using MiniChannelFlow
using Printf
using CairoMakie

mkpath(joinpath(@__DIR__, "output"))
αs = [1.0, 0.5, 0.25, 0.125, 0.0625]

println("   α*    fRe(exact)  fRe(poly)  fRe(FV,nz=32)   Nu_H1(poly)  Nu_H1(FV,nz=48)")
rows = Tuple[]
for α in αs
    ex = exact_fRe(α)
    fv = fRe_numeric(solve_duct(α; nz = 32))
    nuh = nusselt_H1_numeric(solve_duct(α; nz = 48))
    @printf("%6.4f   %9.4f  %9.4f  %12.4f   %11.4f  %14.4f\n",
            α, ex, fRe_fanning(α), fv, nusselt_H1(α), nuh)
    push!(rows, (α, ex, fRe_fanning(α), fv, nusselt_H1(α), nuh))
end
open(joinpath(@__DIR__, "output", "01_fre_nu.csv"), "w") do io
    println(io, "alpha,fRe_exact,fRe_poly,fRe_fv,NuH1_poly,NuH1_fv")
    for r in rows
        println(io, join(r, ","))
    end
end

# --- grid convergence of fRe (square duct and α* = 0.25)
nzs = [8, 12, 16, 24, 32, 48, 64]
fig = Figure(size = (900, 380))
ax = Axis(fig[1, 1]; xlabel = "cells across short side, nz", ylabel = "relative error in fRe",
          xscale = log10, yscale = log10, title = "Second-order convergence")
for α in (1.0, 0.25)
    ex = exact_fRe(α)
    err = [abs(fRe_numeric(solve_duct(α; nz = n)) / ex - 1) for n in nzs]
    scatterlines!(ax, nzs, err; label = "α* = $α")
end
lines!(ax, nzs, 2.0 .* (nzs ./ nzs[1]) .^ (-2) .* 0.05; linestyle = :dash, color = :gray, label = "∝ nz⁻²")
axislegend(ax)

# --- velocity field of a flat mini-channel
sec = solve_duct(0.25; nz = 32)
U = normalized_velocity(sec)
ax2 = Axis(fig[1, 2]; title = "u/ū, α* = 0.25 (finite volume)", xlabel = "y / w", ylabel = "z / w",
           aspect = DataAspect())
hm = heatmap!(ax2, range(0, 1; length = sec.ny), range(0, 0.25; length = sec.nz), U)
Colorbar(fig[1, 2, Right()], hm)
save(joinpath(@__DIR__, "output", "01_fully_developed.png"), fig)
println("figure written to examples/output/01_fully_developed.png")

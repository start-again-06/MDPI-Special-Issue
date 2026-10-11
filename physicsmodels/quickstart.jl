using MiniChannelFlow
using Printf

fluid = water(25.0)
ch = RectChannel(1.0e-3, 2.0e-3, 0.05)       # 1 mm × 2 mm × 50 mm
Dh = hydraulic_diameter(ch)
α = aspect_ratio(ch)
u = 0.5                                       # mean velocity [m/s]
Re = reynolds_number(ch, fluid, u)

@printf("fluid       : %s\n", fluid)
@printf("channel     : Dh = %.3f mm, α* = %.2f, class = %s\n", 1e3Dh, α, classify_channel(ch))
@printf("Reynolds    : Re = %.0f (laminar if < ~2300)\n", Re)
@printf("fully dev.  : fRe = %.3f (exact %.3f), Nu_H1 = %.3f, Nu_T = %.3f\n",
        fRe_fanning(α), exact_fRe(α), nusselt_H1(α), nusselt_T(α))

Lh = hydrodynamic_entrance_length(Re) * Dh
Lt = thermal_entrance_length(Re, prandtl(fluid)) * Dh
@printf("entrance    : hydrodynamic ≈ %.1f mm, thermal ≈ %.0f mm (rule of thumb 0.05 Re Pr Dh)\n",
        1e3Lh, 1e3Lt)

Δp = pressure_drop(ch, fluid, u)
Δp_fd = 2 * fRe_fanning(α) * fluid.μ * u * ch.length / Dh^2
@printf("pressure    : Δp = %.0f Pa with developing-flow correction (fully developed: %.0f Pa)\n", Δp, Δp_fd)

# Thermal entrance length for this duct from the Graetz solver
te = thermal_entrance(α; nz = 24, x_end = 0.3)
xt = thermal_entrance_length_numeric(te)
@printf("graetz      : x*_t = %.3f  →  L_t = %.0f mm;  Nu_T(∞) = %.3f\n",
        xt, 1e3 * xt * Dh * Re * prandtl(fluid), te.Nu[end])

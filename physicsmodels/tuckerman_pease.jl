using MiniChannelFlow
using Printf

cases = [
    # name,   wc,    ww,    depth, flow [cm³/s], power [W], R_meas [K/W], ΔP_meas [kPa]
    ("sample 1", 56e-6, 44e-6, 320e-6, 4.7, 181.0, 0.110, 103.4),   # 15 psi
    ("sample 3", 50e-6, 50e-6, 302e-6, 8.6, 790.0, 0.090, 213.7),   # 31 psi
]
println("case        R_model  R_meas   err%    ΔP_model  ΔP_meas   err%    Re    fin η")
for (name, wc, ww, H, Q, P, Rm, dpm) in cases
    hs = HeatSink(channel_width = wc, wall_width = ww, channel_height = H,
                  base_thickness = 100e-6, length = 0.01, width = 0.01, k_solid = 148.0)
    r = evaluate_heatsink(hs; flow_rate = Q * 1e-6, heat_flux = P / 1e-4)
    @printf("%-9s  %7.4f  %6.3f  %6.1f  %8.1f  %8.1f  %6.1f  %5.0f  %5.2f\n", name,
            r.R_total, Rm, 100 * (r.R_total / Rm - 1), r.dp_total / 1e3, dpm,
            100 * (r.dp_total / 1e3 / dpm - 1), r.Re, r.fin_efficiency)
end
println("\nThe model uses bulk-mean fluid properties; the measured ΔP of sample 3 is lower")
println("because the near-wall fluid is hotter (less viscous) than the bulk.")

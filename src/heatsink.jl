Base.@kwdef struct HeatSink
    channel_width::Float64
    wall_width::Float64
    channel_height::Float64
    base_thickness::Float64
    length::Float64
    width::Float64
    k_solid::Float64
    K_contraction::Float64 = 0.5
    K_expansion::Float64 = 0.0
end

"""
    modify(hs::HeatSink; kwargs...) -> HeatSink

Copy of `hs` with the given fields replaced, e.g.
`modify(hs; channel_width = 60e-6)`.
"""
function modify(hs::HeatSink; kwargs...)
    d = Dict{Symbol,Any}(f => getfield(hs, f) for f in fieldnames(HeatSink))
    for (k, v) in kwargs
        d[k] = v
    end
    return HeatSink(; d...)
end

struct HeatSinkResult
    n_channels::Int
    velocity::Float64
    Re::Float64
    laminar::Bool
    h::Float64
    fin_efficiency::Float64
    T_fluid_mean::Float64
    dT_caloric::Float64
    dT_convective::Float64
    dT_base::Float64
    dT_max::Float64
    R_total::Float64
    dp_channel::Float64
    dp_total::Float64
end

function evaluate_heatsink(hs::HeatSink; flow_rate::Real, heat_flux::Real,
                           T_in::Real = 20.0, nusselt::Symbol = :H1,
                           fluid_model = water)
    flow_rate > 0 || throw(ArgumentError("flow_rate must be positive"))
    heat_flux >= 0 || throw(ArgumentError("heat_flux must be non-negative"))
    nusselt in (:H1, :T) || throw(ArgumentError("nusselt must be :H1 or :T"))

    wc, ww, H = hs.channel_width, hs.wall_width, hs.channel_height
    pitch = wc + ww
    N = floor(Int, hs.width / pitch + 1e-9)
    N >= 1 || throw(ArgumentError("heat sink width holds no channel"))
    Q = float(flow_rate)
    u = Q / (N * wc * H)
    qtot = heat_flux * hs.width * hs.length

    # mean-temperature fluid properties (a few fixed-point passes)
    Tm = T_in + 0.5 * qtot / (4.18e6 * Q)
    fl = fluid_model(clamp(Tm, 0.0, 100.0))
    dTc = qtot / (fl.ρ * fl.cp * Q)
    for _ in 1:3
        Tm = T_in + 0.5 * dTc
        fl = fluid_model(clamp(Tm, 0.0, 100.0))
        dTc = qtot / (fl.ρ * fl.cp * Q)
    end

    α = aspect_ratio(wc, H)
    Dh = hydraulic_diameter(wc, H)
    Nu = nusselt === :H1 ? nusselt_H1(α) : nusselt_T(α)
    h = Nu * fl.k / Dh
    m = sqrt(2h / (hs.k_solid * ww))
    mH = m * H
    η = tanh(mH) / mH
    G = h * (wc + 2η * H)                 # W/(m K) per unit cell, per unit length
    dTconv = heat_flux * pitch / G
    dTbase = heat_flux * hs.base_thickness / hs.k_solid

    ch = RectChannel(wc, H, hs.length)
    Re = fl.ρ * u * Dh / fl.μ
    dpc = pressure_drop(ch, fl, u)
    dpm = (hs.K_contraction + hs.K_expansion) * 0.5 * fl.ρ * u^2
    dTmax = dTc + dTconv + dTbase
    R = qtot > 0 ? dTmax / qtot : NaN
    return HeatSinkResult(N, u, Re, Re < 2300, h, η, Tm, dTc, dTconv, dTbase,
                          dTmax, R, dpc, dpc + dpm)
end

function flow_rate_for_pressure_drop(hs::HeatSink, Δp::Real; heat_flux::Real,
                                     T_in::Real = 20.0, nusselt::Symbol = :H1,
                                     fluid_model = water)
    Δp > 0 || throw(ArgumentError("Δp must be positive"))
    lo, hi = log(1e-12), log(1e-2)
    for _ in 1:100
        mid = 0.5 * (lo + hi)
        r = evaluate_heatsink(hs; flow_rate = exp(mid), heat_flux = heat_flux,
                              T_in = T_in, nusselt = nusselt, fluid_model = fluid_model)
        if r.dp_total > Δp
            hi = mid
        else
            lo = mid
        end
    end
    return exp(0.5 * (lo + hi))
end

function optimize_channel_width(hs::HeatSink; pressure_drop::Real, heat_flux::Real,
                                T_in::Real = 20.0,
                                wc_range::Tuple{<:Real,<:Real} = (20e-6, 500e-6),
                                n_scan::Integer = 48, laminar_only::Bool = true,
                                nusselt::Symbol = :H1, fluid_model = water)
    function eval_at(wc)
        h2 = modify(hs; channel_width = wc)
        Q = flow_rate_for_pressure_drop(h2, pressure_drop; heat_flux = heat_flux,
                                        T_in = T_in, nusselt = nusselt,
                                        fluid_model = fluid_model)
        r = evaluate_heatsink(h2; flow_rate = Q, heat_flux = heat_flux, T_in = T_in,
                              nusselt = nusselt, fluid_model = fluid_model)
        return r, Q
    end
    cost(wc) = begin
        r, _ = eval_at(wc)
        (laminar_only && !r.laminar) ? Inf : r.R_total
    end

    lw0, lw1 = log(wc_range[1]), log(wc_range[2])
    grid = range(lw0, lw1; length = n_scan)
    vals = [cost(exp(lw)) for lw in grid]
    ibest = argmin(vals)
    isfinite(vals[ibest]) ||
        error("no laminar design in the scanned range; widen `wc_range` or set laminar_only = false")

    a = grid[max(ibest - 1, 1)]
    b = grid[min(ibest + 1, n_scan)]
    φ = (sqrt(5) - 1) / 2
    c = b - φ * (b - a)
    d = a + φ * (b - a)
    fc = cost(exp(c))
    fd = cost(exp(d))
    for _ in 1:40
        if fc < fd
            b, d, fd = d, c, fc
            c = b - φ * (b - a)
            fc = cost(exp(c))
        else
            a, c, fc = c, d, fd
            d = a + φ * (b - a)
            fd = cost(exp(d))
        end
    end
    wopt = exp(0.5 * (a + b))
    r, Q = eval_at(wopt)
    return (channel_width = wopt, R_total = r.R_total, flow_rate = Q, result = r)
end

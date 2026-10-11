function _check_alpha(α)
    (0.0 <= α <= 1.0) ||
        throw(DomainError(α, "aspect ratio α* must satisfy 0 ≤ α* ≤ 1"))
    return float(α)
end

function fRe_fanning(α::Real)
    a = _check_alpha(α)
    return 24 * (1 - 1.3553a + 1.9467a^2 - 1.7012a^3 + 0.9564a^4 - 0.2537a^5)
end

"Darcy friction factor × Re, `f_D Re = 4 f_F Re`."
fRe_darcy(α::Real) = 4 * fRe_fanning(α)

function nusselt_T(α::Real)
    a = _check_alpha(α)
    return 7.541 * (1 - 2.610a + 4.970a^2 - 5.119a^3 + 2.702a^4 - 0.548a^5)
end

function nusselt_H1(α::Real)
    a = _check_alpha(α)
    return 8.235 * (1 - 2.0421a + 3.0853a^2 - 2.4765a^3 + 1.0578a^4 - 0.1861a^5)
end

function apparent_fRe(α::Real, Lstar::Real)
    (0.0 < α <= 1.0) || throw(DomainError(α, "need 0 < α* ≤ 1"))
    Lstar > 0 || throw(DomainError(Lstar, "need L/(Dh Re) > 0"))
    s = (1 + α) / (2 * sqrt(α))          # √A / Dh for a rectangle
    Lplus = Lstar / s^2                    # L / (√A Re_√A) = Lstar Dh² / A
    fReA = fRe_fanning(α) * s              # f Re_√A
    return sqrt((3.44 / sqrt(Lplus))^2 + fReA^2) / s
end


hydrodynamic_entrance_length(Re::Real) = 0.315 / (1 + 0.0175Re) + 0.011Re

thermal_entrance_length(Re::Real, Pr::Real) = 0.05 * Re * Pr

function pressure_drop(c::RectChannel, f::Fluid, u::Real)
    Dh = hydraulic_diameter(c)
    Re = f.ρ * u * Dh / f.μ
    Fapp = apparent_fRe(aspect_ratio(c), c.length / (Dh * Re))
    return 2 * Fapp * f.μ * u * c.length / Dh^2
end

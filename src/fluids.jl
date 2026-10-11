struct Fluid
    name::String
    ρ::Float64
    μ::Float64
    k::Float64
    cp::Float64
    function Fluid(name, ρ, μ, k, cp)
        all(x -> x > 0, (ρ, μ, k, cp)) ||
            throw(ArgumentError("all fluid properties must be positive"))
        return new(String(name), Float64(ρ), Float64(μ), Float64(k), Float64(cp))
    end
end

"Kinematic viscosity ν = μ/ρ [m²/s]."
kinematic_viscosity(f::Fluid) = f.μ / f.ρ
"Thermal diffusivity α = k/(ρ cp) [m²/s]."
thermal_diffusivity(f::Fluid) = f.k / (f.ρ * f.cp)
"Prandtl number Pr = cp μ / k."
prandtl(f::Fluid) = f.cp * f.μ / f.k

function water(T::Real)
    (0.0 <= T <= 100.0) ||
        throw(ArgumentError("water(T) is valid for 0 ≤ T ≤ 100 °C, got T = $T"))
    Tc = float(T)
    Tk = Tc + 273.15
    ρ = (999.83952 + 16.945176 * Tc - 7.9870401e-3 * Tc^2 - 46.170461e-6 * Tc^3 +
         105.56302e-9 * Tc^4 - 280.54253e-12 * Tc^5) / (1 + 16.879850e-3 * Tc)
    μ = 2.414e-5 * 10.0^(247.8 / (Tk - 140.0))
    k = -0.58166 + 6.3555e-3 * Tk - 7.9643e-6 * Tk^2
    cp = 4217.4 - 3.7206 * Tc + 0.1413 * Tc^2 - 2.654e-3 * Tc^3 + 2.093e-5 * Tc^4
    return Fluid("water@$(Tc)C", ρ, μ, k, cp)
end

"Water at 20 °C."
const WATER_20C = water(20.0)

function Base.show(io::IO, f::Fluid)
    print(io, "Fluid(", f.name, "; ρ=", round(f.ρ; digits = 2), ", μ=",
          round(f.μ * 1e3; digits = 4), " mPa·s, k=", round(f.k; digits = 4),
          ", cp=", round(f.cp; digits = 1), ", Pr=", round(prandtl(f); digits = 2), ")")
end

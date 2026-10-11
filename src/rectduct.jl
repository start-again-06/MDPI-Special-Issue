function _phi_mean_series(w, h, nterms)
    s = 0.0
    for m in 1:2:nterms, n in 1:2:nterms
        s += 64.0 / (π^6 * m^2 * n^2 * (m^2 / w^2 + n^2 / h^2))
    end
    return s
end

function exact_fRe(α::Real; nterms::Integer = 199)
    (0.0 < α <= 1.0) || throw(DomainError(α, "need 0 < α* ≤ 1"))
    w, h = 1.0, float(α)
    A = w * h
    P = 2 * (w + h)
    return 8 * A^2 / (P^2 * _phi_mean_series(w, h, nterms))
end

function exact_velocity(y::Real, z::Real, w::Real, h::Real; nterms::Integer = 99)
    φ = 0.0
    for m in 1:2:nterms, n in 1:2:nterms
        c = 16.0 / (π^4 * m * n * (m^2 / w^2 + n^2 / h^2))
        φ += c * sin(m * π * y / w) * sin(n * π * z / h)
    end
    return φ / _phi_mean_series(w, h, nterms)
end

function _lap1d(n::Int, d::Float64)
    main = fill(-2.0, n)
    main[1] = -3.0
    main[n] = -3.0
    return spdiagm(-1 => ones(n - 1), 0 => main, 1 => ones(n - 1)) / d^2
end

struct DuctCrossSection
    α::Float64
    ny::Int
    nz::Int
    dy::Float64
    dz::Float64
    lap::SparseMatrixCSC{Float64,Int}
    fact::Factorization{Float64}
    phi::Vector{Float64}
end

function solve_duct(α::Real; nz::Integer = 32)
    (0.0 < α <= 1.0) || throw(DomainError(α, "need 0 < α* ≤ 1"))
    nz >= 4 || throw(ArgumentError("need at least 4 cells across the short side"))
    a = float(α)
    nzi = Int(nz)
    ny = max(round(Int, nzi / a), nzi)
    dy = 1.0 / ny
    dz = a / nzi
    lap = kron(sparse(1.0I, nzi, nzi), _lap1d(ny, dy)) +
          kron(_lap1d(nzi, dz), sparse(1.0I, ny, ny))
    fact = lu(lap)
    phi = fact \ fill(-1.0, ny * nzi)
    return DuctCrossSection(a, ny, nzi, dy, dz, lap, fact, phi)
end

"Normalised velocity field `u/ū` as an `ny × nz` matrix."
normalized_velocity(s::DuctCrossSection) =
    reshape(s.phi ./ (sum(s.phi) / length(s.phi)), s.ny, s.nz)

function fRe_numeric(s::DuctCrossSection)
    A = s.α
    P = 2 * (1 + s.α)
    return 8 * A^2 / (P^2 * (sum(s.phi) / length(s.phi)))
end

function nusselt_H1_numeric(s::DuctCrossSection)
    n = length(s.phi)
    ψ = s.phi ./ (sum(s.phi) / n)
    Θ = s.fact \ (-ψ)
    A = s.α
    P = 2 * (1 + s.α)
    Dh = 4A / P
    return Dh^2 / (4 * (sum(ψ .* Θ) / n))
end

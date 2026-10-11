Base.@kwdef struct PlaneChannelConfig
    L::Float64 = 16.0
    H::Float64 = 1.0
    U::Float64 = 1.0
    Re::Float64 = 100.0
    Nx::Int = 160
    Ny::Int = 20
    gamma::Float64 = 0.3
    cfl::Float64 = 0.4
    tol::Float64 = 1e-8
    maxiter::Int = 200_000
    check_every::Int = 200
end

struct PlaneChannelSolution
    cfg::PlaneChannelConfig
    u::Matrix{Float64}
    v::Matrix{Float64}
    p::Matrix{Float64}
    dx::Float64
    dy::Float64
    nu::Float64
    dt::Float64
    iterations::Int
    residual::Float64
end

# Pressure Poisson operator: Neumann at inlet and walls, Dirichlet p = 0 at the
# outlet face (ghost value -p).  Column-major ordering: index = i + Nx (j-1).
function _pressure_matrix(Nx::Int, Ny::Int, dx::Float64, dy::Float64)
    mx = fill(-2.0, Nx)
    mx[1] = -1.0
    mx[Nx] = -3.0
    Dxx = spdiagm(-1 => ones(Nx - 1), 0 => mx, 1 => ones(Nx - 1)) / dx^2
    my = fill(-2.0, Ny)
    my[1] = -1.0
    my[Ny] = -1.0
    Dyy = spdiagm(-1 => ones(Ny - 1), 0 => my, 1 => ones(Ny - 1)) / dy^2
    return kron(sparse(1.0I, Ny, Ny), Dxx) + kron(Dyy, sparse(1.0I, Nx, Nx))
end

function _predictor!(us, vs, u, v, Nx, Ny, dx, dy, dt, nu, gam, U)
    # --- u-momentum at interior vertical faces i = 2:Nx (face Nx+1 = outlet)
    @inbounds for j in 1:Ny
        for i in 2:Nx
            uc = u[i, j]
            ujp = j < Ny ? u[i, j+1] : -u[i, Ny]
            ujm = j > 1 ? u[i, j-1] : -u[i, 1]
            ue = 0.5 * (uc + u[i+1, j])
            uw = 0.5 * (u[i-1, j] + uc)
            du2dx = (ue * ue - uw * uw) / dx +
                    gam / dx * (abs(ue) * 0.5 * (uc - u[i+1, j]) -
                                abs(uw) * 0.5 * (u[i-1, j] - uc))
            vn = 0.5 * (v[i-1, j+1] + v[i, j+1])
            vsv = 0.5 * (v[i-1, j] + v[i, j])
            un = 0.5 * (uc + ujp)
            usv = 0.5 * (ujm + uc)
            duvdy = (vn * un - vsv * usv) / dy +
                    gam / dy * (abs(vn) * 0.5 * (uc - ujp) -
                                abs(vsv) * 0.5 * (ujm - uc))
            lap = (u[i+1, j] - 2uc + u[i-1, j]) / dx^2 + (ujp - 2uc + ujm) / dy^2
            us[i, j] = uc + dt * (-du2dx - duvdy + nu * lap)
        end
    end
    # --- v-momentum at interior horizontal faces j = 2:Ny
    @inbounds for j in 2:Ny
        for i in 1:Nx
            vc = v[i, j]
            vip = i < Nx ? v[i+1, j] : v[Nx, j]
            vim = i > 1 ? v[i-1, j] : -v[1, j]
            vn = 0.5 * (vc + v[i, j+1])
            vsv = 0.5 * (v[i, j-1] + vc)
            dv2dy = (vn * vn - vsv * vsv) / dy +
                    gam / dy * (abs(vn) * 0.5 * (vc - v[i, j+1]) -
                                abs(vsv) * 0.5 * (v[i, j-1] - vc))
            ue_ = 0.5 * (u[i+1, j-1] + u[i+1, j])
            uw_ = 0.5 * (u[i, j-1] + u[i, j])
            ve = 0.5 * (vc + vip)
            vw = 0.5 * (vim + vc)
            duvdx = (ue_ * ve - uw_ * vw) / dx +
                    gam / dx * (abs(ue_) * 0.5 * (vc - vip) -
                                abs(uw_) * 0.5 * (vim - vc))
            lap = (vip - 2vc + vim) / dx^2 + (v[i, j+1] - 2vc + v[i, j-1]) / dy^2
            vs[i, j] = vc + dt * (-duvdx - dv2dy + nu * lap)
        end
    end
    # --- convective outflow for u* at the outlet face
    @inbounds for j in 1:Ny
        us[Nx+1, j] = u[Nx+1, j] - dt * U * (u[Nx+1, j] - u[Nx, j]) / dx
    end
    return nothing
end

function _divergence!(R, us, vs, Nx, Ny, dx, dy, dt)
    @inbounds for j in 1:Ny, i in 1:Nx
        R[i, j] = ((us[i+1, j] - us[i, j]) / dx + (vs[i, j+1] - vs[i, j]) / dy) / dt
    end
    return nothing
end

function _correct!(u, v, us, vs, p, Nx, Ny, dx, dy, dt)
    @inbounds for j in 1:Ny
        for i in 2:Nx+1
            pe = i <= Nx ? p[i, j] : -p[Nx, j]      # outlet ghost: p = 0 at face
            u[i, j] = us[i, j] - dt * (pe - p[i-1, j]) / dx
        end
    end
    @inbounds for j in 2:Ny, i in 1:Nx
        v[i, j] = vs[i, j] - dt * (p[i, j] - p[i, j-1]) / dy
    end
    return nothing
end

function solve_plane_channel(cfg::PlaneChannelConfig = PlaneChannelConfig();
                             verbose::Bool = false)
    (; L, H, U, Re, Nx, Ny, gamma, cfl, tol, maxiter, check_every) = cfg
    Nx >= 4 && Ny >= 4 || throw(ArgumentError("need Nx ≥ 4 and Ny ≥ 4"))
    nu = U * 2H / Re
    dx = L / Nx
    dy = H / Ny
    F = lu(_pressure_matrix(Nx, Ny, dx, dy))

    u = fill(U, Nx + 1, Ny)
    v = zeros(Nx, Ny + 1)
    us = copy(u)
    vs = copy(v)
    R = zeros(Nx, Ny)
    p = zeros(Nx, Ny)
    uold = similar(u)

    Umax = 1.6U      # fully developed centre-line velocity is 1.5 U
    dt = cfl * min(dx / Umax, 3dy / Umax,
                   0.5 / (nu * (1 / dx^2 + 1 / dy^2)),
                   2nu / Umax^2)

    res = Inf
    it = 0
    while it < maxiter
        it += 1
        copyto!(uold, u)
        u[1, :] .= U
        us[1, :] .= U
        vs[:, 1] .= 0.0
        vs[:, Ny+1] .= 0.0
        _predictor!(us, vs, u, v, Nx, Ny, dx, dy, dt, nu, gamma, U)
        _divergence!(R, us, vs, Nx, Ny, dx, dy, dt)
        p = reshape(F \ vec(R), Nx, Ny)
        _correct!(u, v, us, vs, p, Nx, Ny, dx, dy, dt)
        if it % check_every == 0
            res = maximum(abs.(u .- uold)) / dt
            verbose && println("iter ", it, "  residual ", res)
            res < tol && break
        end
    end
    return PlaneChannelSolution(cfg, u, v, p, dx, dy, nu, dt, it, res)
end

"""
    centerline_velocity(sol) -> Vector

Centre-line axial velocity `u(x, H/2)` at the `Nx+1` vertical faces.
"""
function centerline_velocity(sol::PlaneChannelSolution)
    Ny = sol.cfg.Ny
    if iseven(Ny)
        return 0.5 .* (sol.u[:, Ny÷2] .+ sol.u[:, Ny÷2+1])
    else
        return sol.u[:, (Ny + 1) ÷ 2]
    end
end

"Cross-stream-averaged pressure at the `Nx` cell-centre abscissae."
pressure_profile(sol::PlaneChannelSolution) = vec(sum(sol.p; dims = 2)) ./ sol.cfg.Ny

# plain least-squares line  y ≈ a + b x ;  returns (a, b)
function _linfit(x::AbstractVector, y::AbstractVector)
    xm = sum(x) / length(x)
    ym = sum(y) / length(y)
    b = sum((x .- xm) .* (y .- ym)) / sum((x .- xm) .^ 2)
    return ym - b * xm, b
end

function fully_developed_fRe(sol::PlaneChannelSolution; frac::Real = 0.6)
    cfg = sol.cfg
    x = ((1:cfg.Nx) .- 0.5) .* sol.dx
    mask = x .> frac * cfg.L
    _, slope = _linfit(x[mask], pressure_profile(sol)[mask])
    Dh = 2 * cfg.H
    return -slope * Dh / (2 * cfg.U^2) * cfg.Re
end

function entrance_length(sol::PlaneChannelSolution; fraction::Real = 0.99)
    ucl = centerline_velocity(sol)
    target = fraction * ucl[end-2]
    i = findfirst(>=(target), ucl)
    i === nothing && return Inf
    return (i - 1) * sol.dx / (2 * sol.cfg.H)
end

function incremental_pressure_drop(sol::PlaneChannelSolution)
    cfg = sol.cfg
    Nx, Ny = cfg.Nx, cfg.Ny
    uc = 0.5 .* (sol.u[1:Nx, :] .+ sol.u[2:Nx+1, :])      # cell-centred u
    M = vec(sum(uc .^ 2; dims = 2)) .* sol.dy              # momentum flux ∫u² dy
    M0 = cfg.U^2 * cfg.H
    tau = sol.nu .* uc[:, 1] ./ (sol.dy / 2)               # wall shear (ρ = 1)
    cum = cumsum(tau) .* sol.dx
    dpstar = ((M .- M0) ./ cfg.H .+ 2 .* cum ./ cfg.H) ./ (0.5 * cfg.U^2)
    x = ((1:Nx) .- 0.5) .* sol.dx
    K = dpstar .- 96 .* x ./ (2 * cfg.H * cfg.Re)
    return x, K
end

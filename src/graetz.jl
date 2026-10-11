struct ThermalEntrance
    α::Float64
    xstar::Vector{Float64}
    theta_b::Vector{Float64}
    xmid::Vector{Float64}
    Nu::Vector{Float64}
    Nu_mean::Vector{Float64}
end

function thermal_entrance(α::Real; nz::Integer = 24, x_end::Real = 0.3,
                          per_decade::Integer = 40, x0::Real = 1e-7,
                          x_switch::Real = 1e-2, dx_late::Real = 0.002)
    sec = solve_duct(α; nz = nz)
    n = sec.ny * sec.nz
    A = sec.α
    P = 2 * (1 + sec.α)
    Dh = 4A / P
    ψ = sec.phi ./ (sum(sec.phi) / n)
    Ld = sec.lap .* Dh^2                       # Laplacian in units of Dh
    Ψ = spdiagm(0 => ψ)

    θ = ones(n)
    xs = Float64[0.0]
    tb = Float64[sum(ψ .* θ) / n]
    x = 0.0

    # Rannacher start-up: four backward-Euler sub-steps from 0 to x0
    dx0 = x0 / 4
    Fb = lu(Ψ ./ dx0 - Ld)
    for _ in 1:4
        θ = Fb \ ((ψ ./ dx0) .* θ)
        x += dx0
        push!(xs, x)
        push!(tb, sum(ψ .* θ) / n)
    end

    # logarithmic stages: [xa, 10 xa] split into `per_decade` CN steps
    xa = float(x0)
    while 10 * xa <= x_switch * (1 + 1e-7)
        dx = 9 * xa / per_decade
        Fi = lu(Ψ ./ dx - 0.5 .* Ld)
        Me = Ψ ./ dx + 0.5 .* Ld
        for _ in 1:per_decade
            θ = Fi \ (Me * θ)
            x += dx
            push!(xs, x)
            push!(tb, sum(ψ .* θ) / n)
        end
        xa *= 10
    end

    # uniform late stage
    if x_end > x
        nlate = ceil(Int, (x_end - x) / dx_late)
        dx = (x_end - x) / nlate
        Fi = lu(Ψ ./ dx - 0.5 .* Ld)
        Me = Ψ ./ dx + 0.5 .* Ld
        for _ in 1:nlate
            θ = Fi \ (Me * θ)
            x += dx
            push!(xs, x)
            push!(tb, sum(ψ .* θ) / n)
        end
    end

    lt = log.(tb)
    dxs = diff(xs)
    Nu = -0.25 .* diff(lt) ./ dxs
    xmid = 0.5 .* (xs[1:end-1] .+ xs[2:end])
    Nu_mean = -0.25 .* lt[2:end] ./ xs[2:end]
    return ThermalEntrance(sec.α, xs, tb, xmid, Nu, Nu_mean)
end

function thermal_entrance_length_numeric(res::ThermalEntrance; tol::Real = 0.02)
    Ninf = res.Nu[end]
    for i in eachindex(res.Nu)
        if abs(res.Nu[i] - Ninf) <= tol * Ninf
            return res.xmid[i]
        end
    end
    return res.xmid[end]
end

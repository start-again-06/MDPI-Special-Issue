struct RectChannel
    width::Float64
    height::Float64
    length::Float64
    function RectChannel(width, height, length)
        all(x -> x > 0, (width, height, length)) ||
            throw(ArgumentError("channel dimensions must be positive"))
        return new(Float64(width), Float64(height), Float64(length))
    end
end

"Hydraulic diameter `Dh = 4A/P` of a `w × h` rectangle."
hydraulic_diameter(w::Real, h::Real) = 2 * w * h / (w + h)
hydraulic_diameter(c::RectChannel) = hydraulic_diameter(c.width, c.height)

"Aspect ratio `α* = min(w,h)/max(w,h)` ∈ (0, 1]."
aspect_ratio(w::Real, h::Real) = min(w, h) / max(w, h)
aspect_ratio(c::RectChannel) = aspect_ratio(c.width, c.height)

cross_section_area(c::RectChannel) = c.width * c.height
wetted_perimeter(c::RectChannel) = 2 * (c.width + c.height)

"Reynolds number `Re = ρ u Dh / μ` for mean velocity `u` [m/s]."
reynolds_number(c::RectChannel, f::Fluid, u::Real) =
    f.ρ * u * hydraulic_diameter(c) / f.μ

function classify_channel(Dh::Real)
    Dh > 0 || throw(ArgumentError("hydraulic diameter must be positive"))
    Dh <= 0.1e-6 && return :molecular_nanochannel
    Dh <= 1e-6 && return :transitional_channel
    Dh <= 10e-6 && return :transitional_microchannel
    Dh <= 200e-6 && return :microchannel
    Dh <= 3e-3 && return :minichannel
    return: conventional
end
classify_channel(c::RectChannel) = classify_channel(hydraulic_diameter(c))

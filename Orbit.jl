using Plots, LinearAlgebra

const Rₑ = 6_371.0
const hₚ = 200.0
const μₑ  = 398_600.4418
const ωₑ = 7.2921150e-5

function RotationMatrix(i_deg::Float64, Ω_deg::Float64, ω_deg::Float64)::Matrix{Float64}
    i, Ω, ω = deg2rad.([i_deg, Ω_deg, ω_deg])
    Ci, Si = cos(i), sin(i)
    CΩ, SΩ = cos(Ω), sin(Ω)
    Cω, Sω = cos(ω), sin(ω)
    return [ CΩ*Cω - SΩ*Sω*Ci    -CΩ*Sω - SΩ*Cω*Ci     SΩ*Si ;
             SΩ*Cω + CΩ*Sω*Ci    -SΩ*Sω + CΩ*Cω*Ci    -CΩ*Si ;
                 Sω*Si                Cω*Si              Ci   ]
end

function Perifocal(a::Float64, e::Float64, θs::AbstractVector{<:Real})
    p̲ = a * (1 .- e^2) ./ (1 .+ e .* cos.(θs))
    return hcat(p̲ .* cos.(θs), p̲ .* sin.(θs), zeros(length(θs)))'
end

function TrueAnomaly(e::Float64, N::Int)
    θs = Float64[]
    for M in range(0, 2π; length=N)
        E = M # Initial guess for eccentric anomaly
        for _ in 1:6 # Newton-Raphson method
            E -= (E - e*sin(E) - M) / (1 - e*cos(E))
        end
        θ = 2atan(√((1+e)/(1-e)) * tan(E/2))
        push!(θs, mod(θ, 2π))
    end
    return θs
end

function rECI(a, e, i, Ω, ω; N::Int=1000)
    θs = range(0, 2π; length=N)
    r_PQW = Perifocal(a, e, θs)
    return RotationMatrix(i, Ω, ω) * r_PQW
end

function Sphere(R::Float64 = Rₑ, n::Int = 50)
    θs = range(0, π; length=n)
    ϕs = range(0, 2π; length=2n)
    x̲ = [R * sin(θ) * cos(ϕ) for θ in θs, ϕ in ϕs]
    y̲ = [R * sin(θ) * sin(ϕ) for θ in θs, ϕ in ϕs]
    z̲ = [R * cos(θ)          for θ in θs, ϕ in ϕs]
    return x̲, y̲, z̲
end

Orbit = Dict(
    :a => 15_000.0,
    :e => 0.50,
    :i => 30.0,
    :Ω => 45.0,
    :ω => 30.0
)

r̲ₒ = rECI(Orbit[:a], Orbit[:e], Orbit[:i], Orbit[:Ω], Orbit[:ω])
x̲ₒ, y̲ₒ, z̲ₒ = eachrow(r̲ₒ)
x̲ₑ, y̲ₑ, z̲ₑ = Sphere()

Frames = 200
θ̲_sat = TrueAnomaly(Orbit[:e], Frames)
p_sat = Orbit[:a] * (1 .- Orbit[:e]^2) ./ (1 .+ Orbit[:e] .* cos.(θ̲_sat))
r_PQW_sat = hcat(p_sat .* cos.(θ̲_sat), p_sat .* sin.(θ̲_sat), zeros(Frames))'
r_sat = RotationMatrix(Orbit[:i], Orbit[:Ω], Orbit[:ω]) * r_PQW_sat

T_orb = 2π * sqrt(Orbit[:a]^3 / μₑ)
Δt = T_orb / Frames
lat = zeros(Frames)
lon = zeros(Frames)

λ₀ = -π - atan(r_sat[2,1], r_sat[1,1])

for k in 1:Frames
    r = r_sat[:, k]
    t = (k-1) * Δt
    ϕ = asin(r[3] / norm(r))
    λ = atan(r[2], r[1])
    λ_ECEF = λ + λ₀ - ωₑ * t
    λ_deg = mod(rad2deg(λ_ECEF) + 180, 360) - 180
    lat[k] = rad2deg(ϕ)
    lon[k] = λ_deg
end

anim = @animate for k in 1:Frames
    Azim = 360 * (k-1) / Frames
    Elev = 30.0

    # Earth orbit plot
    Plt3D = plot(legend=false, aspect_ratio=:equal,
                 camera=(Azim, Elev),
                 xlims=(-2Orbit[:a], 2Orbit[:a]),
                 ylims=(-2Orbit[:a], 2Orbit[:a]),
                 zlims=(-2Orbit[:a], 2Orbit[:a]),
                 framestyle=:none, showaxis=false)

    surface!(Plt3D, x̲ₑ, y̲ₑ, z̲ₑ, color=:lightblue, alpha=0.5)
    plot!(Plt3D, x̲ₒ, y̲ₒ, z̲ₒ, linewidth=2, color=:red)
    scatter!(Plt3D, [r_sat[1,k]], [r_sat[2,k]], [r_sat[3,k]], color=:gold, markersize=6)

    # ECI axis lines and labels
    L = 1.5Orbit[:a]
    quiver!(Plt3D, [0.0], [0.0], [0.0], quiver=([L], [0.0], [0.0]), color=:blue, linewidth=2)
    quiver!(Plt3D, [0.0], [0.0], [0.0], quiver=([0.0], [L], [0.0]), color=:blue, linewidth=2)
    quiver!(Plt3D, [0.0], [0.0], [0.0], quiver=([0.0], [0.0], [L]), color=:blue, linewidth=2)

    annotate!(Plt3D, [(1.3L, 0.0, 0.0, text("e₁", :blue, 12)),
                      (0.0, 1.3L, 0.0, text("e₂", :blue, 12)),
                      (0.0, 0.0, 1.3L, text("e₃", :blue, 12))])

    # Ground track plot
    Plt2D = plot(xlims=(-180, 180), ylims=(-90, 90), legend=false,
                 xlabel="Longitude", ylabel="Latitude",
                 framestyle=:box, aspect_ratio=:equal)

    for j in 2:k
        Δlon = abs(lon[j] - lon[j-1])
        if Δlon < 180
            plot!(Plt2D, lon[j-1:j], lat[j-1:j], color=:red, linewidth=1.5)
        end
    end
    scatter!(Plt2D, [lon[k]], [lat[k]], color=:gold, markersize=4)

    plot(Plt3D, Plt2D, layout = @layout [a{0.7h}; b{0.3h}])
end

gif(anim, joinpath(@__DIR__, "Orbit.gif"), fps=30)
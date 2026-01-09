using JuMP
using CairoMakie
using Random
using LinearAlgebra
using PoissonDiskSampling

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")
include("polygon.jl")

Random.seed!(7); # 3, 7 (5 fail)

f = CairoMakie.Figure(size = (512, 920))
ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())



nh = 64

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-8, 0, 0, 0, 0)
cdt_f = ( 8, 0, 0, 0, 0)

# matrices et vecteurs représentatifs des polygones
# collision using polyhedra
poly_C = []
poly_d = []

# collision using polytop
poly = []

N_faces = 5
poly_area = 1;
center = PoissonDiskSampling.generate(2.8, (-6, 6), (-4, 4))
N_poly = length(center)

function normalize(u)
    return u ./ norm(u)
end

function dir(a, b)
    return normalize(b .- a)
end

function normal(u)
    return [-u[2], u[1]]
end

function dot(u, v)
    return sum(u .* v)
end

function pointsToEdge(p1, p2, p3)
    norm = normal(dir(p1, p2));
    norm .*= dot(dir(p3, p1), norm) >= 0 ? 1 : -1
    d = dot(p1, norm)
    return norm, d
end

#=
for i in 1:N_poly
    local center = rand(Float64, 2)
    center .-= 0.5
    center .*= [6, 10]

    local points = rand(Float64, (N_faces, 2))
    points .-= 0.5;
    points .*= 3
    points = [Tuple(p .+ center) for p in eachrow(points)]

    
    local edges = [pointsToEdge(points[i], points[(i)%N_faces+1], points[(i+1)%N_faces+1]) for i in eachindex(points)]
    append!(poly_C, (n -> n[1]).(edges))
    append!(poly_d, (n -> n[2]).(edges))

    poly!(ax, points, strokecolor=:blue, strokewidth=1, color=:white)
end=#

for c in center
    points = convex_volume_polygon(N_faces, poly_area)
    points = [p .+ c for p in points]

    # collision polyhedra
    edges = [pointsToEdge(points[i], points[(i)%N_faces+1], points[(i+1)%N_faces+1]) for i in eachindex(points)]
    append!(poly_C, (n -> n[1]).(edges))
    append!(poly_d, (n -> n[2]).(edges))

    # collision polytop
    push!(poly, points)

    poly!(ax, [Tuple(p) for p in points], strokecolor=:blue, strokewidth=1, color=:white)
end

model = Model()
dynamic_model!(model, nh, cdt_0, cdt_f; m=105)

Npoly_rect_2017_penetration_collisions!(model, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)); kappa=20)
#Npoly_rect_2023_collisions!(model, nh, poly)

solve!(model, max_iter=2000)

#=
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "max_iter", 1000)
JuMP.optimize!(model)
=#

x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
u = JuMP.value.(model[:u]).data
r = JuMP.value.(model[:r]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

# rendering
times = LinRange(0, T, nh+1)

CairoMakie.Label(f[-1, :], "Solution pour traverser un champ de polygones", fontsize = 18)
CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts\nCondition initiale en ligne droite \nSeulement le basculement dans les virages est considéré", justification = :left, fontsize = 12, halign=:left)

plot_positions!(ax, pos, (1.128, 0.720), 1)
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])

ax_u = CairoMakie.Axis(f[2, :], xlabel="temps", ylabel="commande u (m/s)")
ax_r = CairoMakie.Axis(f[3, :], xlabel="temps", ylabel="commande r (rad/s)")
CairoMakie.linkxaxes!(ax_u, ax_r)
CairoMakie.lines!(ax_u, times, u)
CairoMakie.lines!(ax_r, times, r)


ul = JuMP.value.(model[:ul]).data
ur = JuMP.value.(model[:ur]).data
ax_ulr = CairoMakie.Axis(f[4, :], xlabel="temps", ylabel="tensions (V)")
CairoMakie.linkxaxes!(ax_ulr, ax_r)
CairoMakie.lines!(ax_ulr, times, ul, label="ul")
CairoMakie.lines!(ax_ulr, times, ur, label="ur")
axislegend(ax_ulr)


#=μ = 0.5
b = 1.128*0.4
Tl = JuMP.value.(model[:Tl]).data
Tr = JuMP.value.(model[:Tr]).data
τl = JuMP.value.(model[:τl]).data
τr = JuMP.value.(model[:τr]).data
nl = JuMP.value.(model[:nl]).data
nr = JuMP.value.(model[:nr]).data

limitr = μ .* (nr)# .- abs.(τl+τr) ./ b)
limitl = μ .* (nl)# .- abs.(τl+τr) ./ b)
ax_τlr = Axis(f[5, :], ylabel="f. tangentielles (N)")
linkxaxes!(ax_τlr, ax_r)
lines = [lines!(ax_τlr, times, norm.(Tl, τl), label="T left"),
lines!(ax_τlr, times, norm.(Tr, τr), label="T right"),
lines!(ax_τlr, times, limitl, label="T_lim left"),
lines!(ax_τlr, times, limitr, label="T_lim right")]

#lines!(ax_τlr, times, u.*r.*10) 
#=lines!(ax_τlr, times, τl .+ τr)
lines!(ax_τlr, times, JuMP.value.(model[:τplus]).data)
lines!(ax_τlr, times, JuMP.value.(model[:τmoins]).data)=#
#lines!(ax_τlr, times, nl)
#lines!(ax_τlr, times, nr)=#

#Legend(f[6, :], lines, ["T left", "T right", "T lim left", "T lim right"], orientation=:horizontal)


f
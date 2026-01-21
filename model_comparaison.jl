using JuMP
using CairoMakie
using Random
using LinearAlgebra
using PoissonDiskSampling

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")
include("polygon.jl")

Random.seed!(14); # 11, 3, 5, 7

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

start = initAstar(nh, poly, (-8.2, 8.2), (-5, 5), cdt_0, cdt_f; spacing=0.3, dmin=0.4)


# collision polyhedra
model1 = Model()
dynamic_model!(model1, nh, cdt_0, cdt_f; start=start, m=105)
Npoly_rect_2012_collisions!(model1, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
solve!(model1, max_iter=1000)


model2 = Model()
dynamic_model!(model2, nh, cdt_0, cdt_f; start=start, m=105)
Npoly_rect_2017_penetration_collisions!(model2, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)); kappa=1e+2)
solve!(model2, max_iter=1000)

# collision polytop
model3 = Model()
dynamic_model!(model3, nh, cdt_0, cdt_f; start=start, m=105)
Npoly_rect_2023_collisions!(model3, nh, poly; d_min=0.05)
solve!(model3, max_iter=1000)

# rendering
x1 = JuMP.value.(model1[:x]).data
y1 = JuMP.value.(model1[:y]).data
#u1 = JuMP.value.(model1[:u]).data
#r1 = JuMP.value.(model1[:r]).data
T1 = JuMP.value.(model1[:T])
N1 = JuMP.barrier_iterations(model1)

x2 = JuMP.value.(model2[:x]).data
y2 = JuMP.value.(model2[:y]).data
#u2 = JuMP.value.(model2[:u]).data
#r2 = JuMP.value.(model2[:r]).data
T2 = JuMP.value.(model2[:T])
N2 = JuMP.barrier_iterations(model2)

x3 = JuMP.value.(model3[:x]).data
y3 = JuMP.value.(model3[:y]).data
#u3 = JuMP.value.(model3[:u]).data
#r3 = JuMP.value.(model3[:r]).data
T3 = JuMP.value.(model3[:T])
N3 = JuMP.barrier_iterations(model3)

#ϕ = JuMP.value.(model[:ϕ]).data
pos1 = ((a, b) -> (a, b)).(x1, y1)
pos2 = ((a, b) -> (a, b)).(x2, y2)
pos3 = ((a, b) -> (a, b)).(x3, y3)

# rendering
f = CairoMakie.Figure(size = (792, 512))
infos = GridLayout()
plots = GridLayout()

#times = LinRange(0, T, nh+1)

#legend = axislegend(ax)#/, [m1, m2, m3], ["modèle 2012", "modèle 2017", "modèle 2024"], alignmode = Inside());
#infos[3, 1] = legend
label = CairoMakie.Label(f, "Nombre de points : $(nh+1) pts \nInitialisaion avec A*", justification = :left, fontsize = 12, halign=:left)
infos[3, 1] = label

for (i, (name, T, N)) in enumerate([("2012", T1, N1), ("2017", T2, N2), ("2023", T3, N3)])
    info = CairoMakie.Label(f, "Modèle : $name \nTemps de trajet : T = $(trunc(T, digits=3, base=10)) s\nIterations = $(N)", justification = :left, fontsize = 12, halign=:left)
    infos[2, i] = info
end
title = CairoMakie.Label(f, "Solution pour un champ de polygones", fontsize = 18)
infos[1, :] = title


ax = Axis(f, aspect = CairoMakie.DataAspect(), alignmode = Inside())
plots[1, 1] = ax

plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])
plot_trajectory!(ax, pos1, col=:red, label="2012")
plot_trajectory!(ax, pos2, col=:orange, label="2017")
plot_trajectory!(ax, pos3, col=:yellow, label="2023")
axislegend(ax)

for points in poly
    poly!(ax, [Tuple(p) for p in points], strokecolor=:blue, strokewidth=1, color=:white)
end
#=
ax_u = CairoMakie.Axis(f, xlabel="temps", ylabel="commande u", alignmode = Inside())
ax_r = CairoMakie.Axis(f, xlabel="temps", ylabel="commande r", alignmode = Inside())
CairoMakie.linkxaxes!(ax_u, ax_r)
plots[2, :] = ax_u
plots[3, :] = ax_r

m1 = CairoMakie.lines!(ax_u, times, u1, color=:red)
m2 = CairoMakie.lines!(ax_u, times, u2, color=:orange)
m3 = CairoMakie.lines!(ax_u, times, u3, color=:yellow)
CairoMakie.lines!(ax_r, times, r1, color=:red)
CairoMakie.lines!(ax_r, times, r2, color=:orange)
CairoMakie.lines!(ax_r, times, r3, color=:yellow)
=#

CairoMakie.colsize!(f.layout, 1, CairoMakie.Relative(1.0))


f.layout[1, 1] = infos
f.layout[2, 1] = plots

save("figure.svg", f)
f
using JuMP
using CairoMakie
using Random
using LinearAlgebra

include("../models.jl")
include("../renderer.jl")

nh = 32

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-1, 0,  π, 1.0, 0)
cdt_f = ( 7, 0,  0, 4, 0)

N_poly = 1
N_faces = 2
lp = 0.2 # largeur du mur
C = [-1 0; 0 -1; 0 1]
d = [lp/2; 1.0; lp/2]
poly = [[[-lp/2, 100], [-lp/2, -1], [lp/2, -1], [lp/2, 100]]]

#start = initAstar(nh, poly, (-3, 5), (-3, 1), cdt_0, cdt_f; spacing=0.2, dmin=0.6, smoothing=6)

model = Model()

#speed_model!(model, nh, cdt_0, cdt_f; rk=CrankNicolson())#; epsilon=1e-4)#, μ=0.5, m=105)
dynamic_model!(model, nh, cdt_0, cdt_f; start=nothing, epsilon=1e-4)#, μ=0.5, m=105)
#dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=start, μ=0.1)


#Npoly_rect_2012_collisions!(model, nh, (1, 2), C, d)
#Npoly_rect_2017_collisions!(model, nh, (1, 2), C, d)
Npoly_rect_2017_penetration_collisions!(model, nh, (1, 2), C, d; kappa=1e+5)
#Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.1)

limit_time!(model, 5.0)
solve!(model, max_iter=1000)


x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

# rendering
f = CairoMakie.Figure(size = (512, 200))
ax = CairoMakie.Axis(f[1, :], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

lines!(ax, [(-lp/2, 4), (-lp/2, -1), (lp/2, -1), (lp/2, 4)], color=:blue)

#CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nlargeur : l = $lp m\nposition : p = $mp m\nN = $(nh+1) pts\nCondition initiale en ligne droite", justification = :left, fontsize = 12, halign=:left)

plot_positions!(ax, pos, (1.128, 0.720), 1)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])
#plot_start!(ax, start)
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)


#CairoMakie.colsize!(f.layout, 1, CairoMakie.Aspect(1, 1.0))

save("figure.svg", f)
f
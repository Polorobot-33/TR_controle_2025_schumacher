using JuMP
using CairoMakie
using Random
using LinearAlgebra

include("models.jl")
include("renderer.jl")


f = CairoMakie.Figure(size = (512, 860))
ax = CairoMakie.Axis(f[1, 1], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

nh = 99

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-2.5, 2.0, 0, 0, 0)
cdt_f = (0, -0.4, 0, 0, 0)

N_poly = 3
N_faces = 2
lp = 2.1 # longueur de la place de parking
pp = 0.8 # profondeur de la place de parking
C = [0 1; 1 0; 0 1; 0 -1; -1 0; 0 1]
d = [0; -lp/2; -pp; pp+10; -lp/2; 0]
lines!(ax, [(-3, 0), (-lp/2, 0), (-lp/2, -pp), (lp/2, -pp), (lp/2, 0), (3, 0)], color=:blue)


model = Model()
speed_model!(model, nh, cdt_0, cdt_f)
Npoly_rect_2017_penetration_collisions!(model, nh, (3, 2), C, d; kappa=100)
solve!(model, max_iter=500)


x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
u = JuMP.value.(model[:u]).data
r = JuMP.value.(model[:r]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

# rendering
times = LinRange(0, T, nh+1)

CairoMakie.Label(f[-1, :], "Solution pour une place de parking", fontsize = 18)
CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nlongueur : L = $lp m\nprofondeur : P = $pp m\nN = $(nh+1) pts\nCondition initiale en ligne droite", justification = :left, fontsize = 12, halign=:left)

plot_positions!(ax, pos, (1.128, 0.720), 2)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)

ax_u = CairoMakie.Axis(f[2, :], xlabel="temps", ylabel="commande u")
ax_r = CairoMakie.Axis(f[3, :], xlabel="temps", ylabel="commande r")
CairoMakie.linkxaxes!(ax_u, ax_r)
CairoMakie.lines!(ax_u, times, u)
CairoMakie.lines!(ax_r, times, r)

#CairoMakie.rowsize!(f.layout, 1, CairoMakie.Aspect(1, 1.0))

save("figure.svg", f)
f
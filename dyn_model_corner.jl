using JuMP
using CairoMakie

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")

f = CairoMakie.Figure(size = (512, 892))
ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

nh = 99

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (0, -3, π/2, 2, 0)
cdt_f = (3,  0, 0,   2,   0)

model = Model()
dynamic_model!(model, nh, cdt_0, cdt_f)

# collisions
l1_c = 2.4 # largeur principale des couloirs
l2_c = 0.9 # largeur des couloirs latéraux
C = [0 -1; -1 0; 1 0; 0 -1; 0 1; 1 0; -1 0; 0 1]
d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
    

# collision polyhedra
Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d; kappa=20)

# collision polytop
# model, init = robot_rect_custom_polytop(nh, cdt_0, cdt_f, poly; d_min=1e-1)

solve!(model, max_iter=1000)
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

CairoMakie.Label(f[-1, :], "Solution pour le virage", fontsize = 18)
CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts\nCondition initiale en ligne droite", justification = :left, fontsize = 12, halign=:left)

plot_terrain!(ax, l1_c, l2_c)
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


CairoMakie.rowsize!(f.layout, 1, Aspect(1, 1.0))

f
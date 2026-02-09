using JuMP
using CairoMakie

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")

nh = 99

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (0, -2.5, π/2, 2.7, 0)
cdt_f = (2.5,  0, 0,   2.7, 0)

# collisions
l1_c = 2.4 # largeur principale des couloirs
l2_c = 0.9 # largeur des couloirs latéraux
C = [0 -1; -1 0; 1 0; 0 -1; 0 1; 1 0; -1 0; 0 1]
d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
poly = [[p .*[dx;dy] for p in [[l1_c/2;l2_c/2], [150;l2_c/2], [150;150], [l1_c/2;150]]] for (dx, dy) in [(1,1), (1,-1), (-1,-1), (-1,1)]]


# modèle
start = initAstar(nh, poly, (-3, 5), (-5, 3), cdt_0, cdt_f; spacing=0.18, dmin=0.3, smoothing=12)

model = Model()
#speed_model!(model, nh, cdt_0, cdt_f; start=start, rk=CrankNicolson())#; epsilon=1e-4)#, μ=0.5, m=105)
dynamic_model!(model, nh, cdt_0, cdt_f; start=start)#, μ=0.5, m=105)
#dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=start, μ=0.1)


#Npoly_rect_2012_collisions!(model, nh, (4, 2), C, d)
#Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d)
Npoly_rect_2017_penetration_collisions!(model, nh, (4, 2), C, d; kappa=1e+2)
#Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)


#=# Additional constraints
u = model[:u]
r = model[:r]
@constraints(model, begin 
    [i=0:nh], -3.0 <= 105 * u[i]*r[i] <= 3.0
    [i=0:nh], -2.7 <= u[i] <= 2.7
    [i=0:nh], -1.5 <= r[i] <= 1.5
end)=#



#limit_time!(model, 2.5)
solve!(model, max_iter=1000)

x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
#a = JuMP.value.(model[:a]).data
#r = JuMP.value.(model[:r]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)


# rendering
f = Figure(size = (600, 1024))
ax = CairoMakie.Axis(f[1, 1], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

times = LinRange(0, T, nh+1)

Label(f[-1, :], "Solution pour le virage, comparaison des différents modèles", fontsize = 18)
Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts\nCondition initiale avec A*", justification = :left, fontsize = 12, halign=:left)

plot_terrain!(ax, l1_c, l2_c)
plot_positions!(ax, pos, (1.128, 0.720), 1)
#pl_ref = plot_ref!(ax, (l1_c-1.228) / 2)
pl_st = plot_start!(ax, start)
pl_traj = plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])


# commandes

#=# ul, ur
ax_u = Axis(f[2, :], xlabel="temps", ylabel="tensions u (V)")
lines!(ax_u, times, var(model, :ul), label="ul")
lines!(ax_u, times, var(model, :ur), label="ur")
axislegend(ax)=#

# u et r
ax_u = Axis(f[2, :], xlabel="temps", ylabel="vitesse u (m/s)")
lines!(ax_u, times, var(model, :u))
ax_r = Axis(f[3, :], xlabel="temps", ylabel="rotation r (rad/s)")
lines!(ax_r, times, var(model, :r))
linkxaxes!(ax_u, ax_r)

# a et r
#=ax_a = Axis(f[2, :], xlabel="temps", ylabel="vitesse a (m/s²)")
lines!(ax_a, times, var(model, :a))
ax_r = Axis(f[3, :], xlabel="temps", ylabel="rotation r (rad/s²)")
lines!(ax_r, times, var(model, :r))
linkaxes!(ax_a, ax_r)=#

rowsize!(f.layout, 1, Aspect(1, 1.0))

#make_animation(x, y, ϕ, (ax) -> plot_terrain!(ax, l1_c, l2_c))
save("figure.svg", f)
f
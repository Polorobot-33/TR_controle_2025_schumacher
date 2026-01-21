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
#speed_model!(model, nh, cdt_0, cdt_f; start=start)#; epsilon=1e-4)#, μ=0.5, m=105)
#dynamic_model!(model, nh, cdt_0, cdt_f; start=start, epsilon=1e-4)#, μ=0.5, m=105)
dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=start, μ=0.1)


#Npoly_rect_2012_collisions!(model, nh, (4, 2), C, d)
Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d)
#Npoly_rect_2017_penetration_collisions!(model, nh, (4, 2), C, d; kappa=1e+2)
#Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)

#limit_time!(model, 3.5)
solve!(model, max_iter=1500)

x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
#a = JuMP.value.(model[:a]).data
#r = JuMP.value.(model[:r]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)


# rendering
f = Figure(size = (560, 920))#(512, 920))
ax = CairoMakie.Axis(f[1, 1], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

times = LinRange(0, T, nh+1)

Label(f[-1, :], "Solution pour le virage, comparaison des différents modèles", fontsize = 18)
Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts\nCondition initiale avec A*", justification = :left, fontsize = 12, halign=:left)

plot_terrain!(ax, l1_c, l2_c)
plot_positions!(ax, pos, (1.128, 0.720), 1)
pl_ref = plot_ref!(ax, (l1_c-1.228) / 2)
#pl_traj = plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])



#=ax_u = Axis(f[2, :], xlabel="temps", ylabel="vitesse u (m/s)")
lines!(ax_u, times, var(model, :u))=#
ax_a = Axis(f[2, :], xlabel="temps", ylabel="acc a (m/s²)")
lines!(ax_a, times, (var(model, :τl) .+ var(model, :τl)) ./ (0.05 * 105))

#Legend(f[1, :], [pl_ref, pl_traj], ["trajectoire de référence", "trajectoire optimisée"], halign=:right, valign=:top, orientation=:horizontal)

#= commandes
ax_a = Axis(f[2, :], xlabel="temps", ylabel="commande a (m/s²)")
ax_r = Axis(f[3, :], xlabel="temps", ylabel="commande r (rad/s²)")
linkxaxes!(ax_a, ax_r)
lines!(ax_a, times, a)
lines!(ax_r, times, r)=#

#ax_F = Axis(f[4, :], xlabel="temps", ylabel="F (N)")#; limits=(nothing, (-0.0001, 0.0001)))
#lines!(ax_F, times, var(model, :xplus) .- var(model, :xmoins), label="x+ - x-")
#lines!(ax_F, times, var(model, :xplus), label="x+")
#lines!(ax_F, times, var(model, :xmoins), label="x-")
#axislegend(ax_F)
#lines!(ax_F, times, -var(model, :sgn))

#lines!(ax_F, times, 0.1 .* var(model,:F))
#lines!(ax_F, times, 0.1 .* (transpose.(var(model,:v)) .* var(model,:e_t) .* var(model,:ϕdot)))
#lines!(ax_F, times, (transpose.(var(model,:v)) .* var(model,:e_n)))


#ax_v = Axis(f[2, :], xlabel="temps")
#ax_ϕdot = Axis(f[3, :], xlabel="temps")
#linkxaxes!(ax_v, ax_ϕdot)
#lines!(ax_v, times, norm.(var(model,:v)))
#lines!(ax_v, times, transpose.(var(model,:v)) .* var(model,:e_n))
#lines!(ax_v, times, transpose.(var(model,:v)) .* var(model,:e_t))
#lines!(ax_ϕdot, times, var(model,:ϕdot))
#=ax_a = Axis(f[2, :], xlabel="temps", ylabel="a (m/s²)")
ax_r = Axis(f[3, :], xlabel="temps", ylabel="r (rad/s²)")
linkxaxes!(ax_a, ax_r)
lines!(ax_a, times, var(model,:a))
lines!(ax_r, times, var(model,:r))=#
#lines!(ax_v, times, transpose.(var(model,:v)) .* var(model,:e_t))

#make_animation(x, y, ϕ, (ax) -> plot_terrain!(ax, l1_c, l2_c))
rowsize!(f.layout, 1, Aspect(1, 1.0))

save("figure.svg", f)
f
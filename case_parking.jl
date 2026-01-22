using JuMP
using CairoMakie
using Random
using LinearAlgebra

include("models.jl")
include("renderer.jl")

nh = 99

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-2.5, 2.0, 0, 0, 0)
cdt_f = (-0.4, -0.4, 0, 0, 0)

N_poly = 3
N_faces = 2
lp = 2.1 # longueur de la place de parking
pp = 0.8 # profondeur de la place de parking
rp = 3.0 # largeur de la route
C = [0 1; 1 0;   0 1; 0 -1;   -1 0; 0 1;   0 -1; 0 1]
d = [0; -lp/2;   -pp; pp+10;   -lp/2; 0;   -rp; rp+10]
poly=[[[-10;0], [-lp/2;0], [-lp/2;-10], [-10;-10]],
      [[-lp/2; -pp], [-lp/2;-pp-10], [lp/2; -pp-10], [lp/2;-pp]],
      [[10;0], [lp/2;0], [lp/2;-10], [10;-10]],
      [[-10;rp], [-10;rp+10], [10;rp+10], [10;rp]]]

start = initAstar(nh, poly, (-4, 2), (-1.5, 3), cdt_0, cdt_f; spacing=0.2, dmin=0.3, smoothing=12)

model = Model()

speed_model!(model, nh, cdt_0, cdt_f; start=start, rk=CrankNicolson())#; epsilon=1e-4)#, μ=0.5, m=105)
#dynamic_model!(model, nh, cdt_0, cdt_f; start=start, epsilon=1e-4)#, μ=0.5, m=105)
#dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=start, μ=0.1)


#Npoly_rect_2012_collisions!(model, nh, (4, 2), C, d)
#Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d)
#Npoly_rect_2017_penetration_collisions!(model, nh, (4, 2), C, d; kappa=1e+2)
Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.01)

limit_time!(model, 15.0)
solve!(model, max_iter=1000)


x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

# rendering
f = CairoMakie.Figure(size = (512, 800))
ax = CairoMakie.Axis(f[1, :], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

lines!(ax, [(-3, 0), (-lp/2, 0), (-lp/2, -pp), (lp/2, -pp), (lp/2, 0), (3, 0)], color=:blue)
lines!(ax, [(-3, rp), (3, rp)], color=:blue)

times = LinRange(0, T, nh+1)

CairoMakie.Label(f[-1, :], "Solution pour une place de parking", fontsize = 18)
CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nlongueur : L = $lp m\nprofondeur : P = $pp m\nN = $(nh+1) pts\nCondition initiale en ligne droite", justification = :left, fontsize = 12, halign=:left)

plot_positions!(ax, pos, (1.128, 0.720), 1)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])
plot_start!(ax, start)
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)





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

CairoMakie.rowsize!(f.layout, 1, CairoMakie.Aspect(1, 1.0))

save("figure.svg", f)
f
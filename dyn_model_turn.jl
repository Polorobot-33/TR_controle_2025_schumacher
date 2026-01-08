using JuMP
using CairoMakie

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")

f = CairoMakie.Figure(size = (512, 892))
ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

nh = 109

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (0,  1, 0, 3, 0)
cdt_f = (3, -3, 0, 3, 0)

model = Model()
dynamic_model_slide!(model, nh, cdt_0, cdt_f; m=10)

# collisions
l1_c = 0.2 # largeur du mur
l2_c = 3 # longueur du mur 
l3_c = 2 # espacement entre les virages
C = [0 1; 1 0; 0 -1;   0 1; -1 0; 0 -1]
d = [l1_c/2; l2_c; l1_c/2;  l1_c/2-l3_c; 0; l1_c/2+l3_c]

# collision polyhedra
Npoly_rect_2017_collisions!(model, nh, (2, 3), C, d; kappa=150)

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
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

# rendering
times = LinRange(0, T, nh+1)

Label(f[-1, :], "Solution pour un virage en épingle", fontsize = 18)
Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts\nCondition initiale en ligne droite", justification = :left, fontsize = 12, halign=:left)

lines!(ax, [(-1, l1_c/2), (l2_c, l1_c/2), (l2_c, -l1_c/2), (-1, -l1_c/2)]; color=:blue)
lines!(ax, [(l2_c+1, l1_c/2-l3_c), (0, l1_c/2-l3_c), (0, -l1_c/2-l3_c), (l2_c+1, -l1_c/2-l3_c)]; color=:blue)
plot_positions!(ax, pos, (1.128, 0.720), 1)
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])

ax_a = Axis(f[2, :], xlabel="temps", ylabel="commande a (m/s²)")
ax_r = Axis(f[3, :], xlabel="temps", ylabel="commande r (rad/s²)")
linkxaxes!(ax_a, ax_r)
lines!(ax_a, times, JuMP.value.(model[:a]).data)
lines!(ax_r, times, JuMP.value.(model[:r]).data)

ax_sgn = Axis(f[4, :], xlabel="temps", ylabel="sgn(n(ϕ)^T v)")
lines!(ax_sgn, times, JuMP.value.(model[:sign]).data)
#lines!(ax_sgn, times, (JuMP.value.(model[:xplus]).data .- JuMP.value.(model[:xmoins]).data))
#lines!(ax_sgn, times, (transpose.(JuMP.value.(model[:e_n]).data) .* JuMP.value.(model[:v]).data))
#=
ul = JuMP.value.(model[:ul]).data
ur = JuMP.value.(model[:ur]).data
ax_ulr = CairoMakie.Axis(f[4, :], xlabel="temps", ylabel="tensions (V)")
linkxaxes!(ax_ulr, ax_r)
lines!(ax_ulr, times, ul, label="ul")
lines!(ax_ulr, times, ur, label="ur")
axislegend(ax_ulr)=#

rowsize!(f.layout, 1, Aspect(1, 1.0))


# animation
framerate = 25
timestamps = 1:(nh+1)

fig = Figure(size=(512, 512))
axanim = Axis(fig[1, 1], aspect = CairoMakie.DataAspect())
lines!(axanim, [(-1, l1_c/2), (l2_c, l1_c/2), (l2_c, -l1_c/2), (-1, -l1_c/2)]; color=:blue)
lines!(axanim, [(l2_c+1, l1_c/2-l3_c), (0, l1_c/2-l3_c), (0, -l1_c/2-l3_c), (l2_c+1, -l1_c/2-l3_c)]; color=:blue)
plot_trajectory!(axanim, [(px, py) for (px, py, _) in pos], col=:red)
time = Observable(1)
robot_anim = @lift(robot_points((x[$time], y[$time], ϕ[$time]), (1.128, 0.720)))
lines!(axanim, robot_anim, color=:green)

record(fig, "time_animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end

f
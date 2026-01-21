using JuMP
using CairoMakie

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")

nh = 99
f_c = 0.1

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (0,  1, 0, 1, 0)
cdt_f = (3, -3.5, 0, 1, 0)


# collisions
l1_c = 0.5 # largeur du mur
l2_c = 3.0 # longueur du mur 
l3_c = 2.5 # espacement entre les virages
C = [0 1; 1 0; 0 -1;   0 1; -1 0; 0 -1]
d = [l1_c/2; l2_c; l1_c/2;  l1_c/2-l3_c; 0; l1_c/2+l3_c]
poly = [[[-10; l1_c/2], [l2_c; l1_c/2], [l2_c; -l1_c/2], [-10; -l1_c/2]],
        [[20; l1_c/2-l3_c], [0; l1_c/2-l3_c], [0; -l1_c/2-l3_c], [20; -l1_c/2-l3_c]]]


# modèle
start = initAstar(nh, poly, (-2, 6), (-4, 2), cdt_0, cdt_f; spacing=0.2, dmin=0.6, smoothing=8)

model = Model()
#dynamic_model!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=start, f=f_c)
#dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=start, μ=f_c)

#Npoly_rect_2012_collisions!(model, nh, (2, 3), C, d)#; epsilon=1e-7)
#Npoly_rect_2017_collisions!(model, nh, (2, 3), C, d)
#Npoly_rect_2017_penetration_collisions!(model, nh, (2, 3), C, d; kappa=1e+2)
#Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)

#limit_time!(model, 20.0)
solve!(model, max_iter=1000)


x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)


# rendering
f = CairoMakie.Figure(size = (512, 512))#(512, 892))
ax = CairoMakie.Axis(f[1, 1], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

times = LinRange(0, T, nh+1)

Label(f[-1, :], "Solution pour un virage en épingle", fontsize = 18)
Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts\nf = $(f_c)", justification = :left, fontsize = 12, halign=:left)


function plot_turn!(ax)
        lines!(ax, [(-1, l1_c/2), (l2_c, l1_c/2), (l2_c, -l1_c/2), (-1, -l1_c/2)]; color=:blue)
        lines!(ax, [(l2_c+1, l1_c/2-l3_c), (0, l1_c/2-l3_c), (0, -l1_c/2-l3_c), (l2_c+1, -l1_c/2-l3_c)]; color=:blue)
end
plot_turn!(ax)
plot_positions!(ax, pos, (1.128, 0.720), 1)
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])

#=ax_u = Axis(f[2, :], xlabel="temps", ylabel="vitesse u (m/s)")
ax_r = Axis(f[3, :], xlabel="temps", ylabel="vitesse r (rad/s)")
linkxaxes!(ax_u, ax_r)
lines!(ax_u, times, JuMP.value.(model[:u]).data)
lines!(ax_r, times, JuMP.value.(model[:r]).data)=#

#=
ax_a = Axis(f[2, :], xlabel="temps", ylabel="commande a (m/s²)")
ax_r = Axis(f[3, :], xlabel="temps", ylabel="commande r (rad/s²)")
linkxaxes!(ax_a, ax_r)
lines!(ax_a, times, JuMP.value.(model[:a]).data)
lines!(ax_r, times, JuMP.value.(model[:r]).data)

ax_F = Axis(f[4, :], xlabel="temps", ylabel="F (N)")
lines!(ax_F, times, JuMP.value.(model[:F]).data)
#lines!(ax_sgn, times, (JuMP.value.(model[:xplus]).data .- JuMP.value.(model[:xmoins]).data))
#lines!(ax_sgn, times, (transpose.(JuMP.value.(model[:e_n]).data) .* JuMP.value.(model[:v]).data))

ul = JuMP.value.(model[:ul]).data
ur = JuMP.value.(model[:ur]).data
ax_ulr = CairoMakie.Axis(f[4, :], xlabel="temps", ylabel="tensions (V)")
linkxaxes!(ax_ulr, ax_r)
lines!(ax_ulr, times, ul, label="ul")
lines!(ax_ulr, times, ur, label="ur")
axislegend(ax_ulr)

rowsize!(f.layout, 1, Aspect(1, 1.0))=#

#=
ax_F = Axis(f[2, :], xlabel="temps", ylabel="F (N)")#; limits=(nothing, (-0.0001, 0.0001)))
lines!(ax_F, times, var(model,:F))
lines!(ax_F, times, (transpose.(var(model,:v)) .* var(model,:e_t) .* var(model,:ϕdot)))
=#

save("figure.svg", f)

f
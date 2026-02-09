using JuMP
using CairoMakie
using Random
using LinearAlgebra
using PoissonDiskSampling

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")
include("polygon.jl")

Random.seed!(7); # 14, 11, 3, 5, 7

f = CairoMakie.Figure(size = (512, 920))#(512, 920))
ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside(), xlabel="position x (m)", ylabel="position y (m)")



nh = 127

# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-8, 0, 0, 0, 0)
cdt_f = ( 8, 0, 0, 0, 0)

# matrices et vecteurs représentatifs des polygones
poly_C = []
poly_d = []
poly = []

N_faces = 5
poly_area = 1;
center = PoissonDiskSampling.generate(2.8, (-6, 6), (-4, 4))
N_poly = length(center)



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

model = Model()
speed_model!(model, nh, cdt_0, cdt_f; start=start, rk=CrankNicolson())#; epsilon=1e-4)#, μ=0.5, m=105)
#dynamic_model!(model, nh, cdt_0, cdt_f; start=start)#, μ=0.5, m=105)
#dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=start)
#dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=start, μ=0.1)

#Npoly_rect_2012_collisions!(model, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
#Npoly_rect_2017_collisions!(model, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
#Npoly_rect_2017_penetration_collisions!(model, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)); kappa=20)
#Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)

solve!(model, max_iter=1000)

x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
T = JuMP.value.(model[:T])
pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

# rendering
times = LinRange(0, T, nh+1)

CairoMakie.Label(f[-1, :], "Solution pour traverser un champ de polygones", fontsize = 18)
CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $(nh+1) pts", justification = :left, fontsize = 12, halign=:left)

plot_positions!(ax, pos, (1.128, 0.720), 1)
plot_start!(ax, start)
plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
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

CairoMakie.rowsize!(f.layout, 1, CairoMakie.Aspect(1, 1.0))

save("figure.svg", f)
f
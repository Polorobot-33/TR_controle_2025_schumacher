using JuMP
using CairoMakie

include("models.jl")
include("initial_conditions.jl")
include("renderer.jl")

nh = 64
u_max = 0.9:0.2:2.7

f = CairoMakie.Figure(size = (512, 512))

CairoMakie.Label(f[-1, :], "Solutions optimales en fonction de U max", fontsize = 18)

ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect())
plot_terrain!(ax, 2.4, 0.9)

l = []
for u in u_max
    local model, init = robot_rect_model(nh, u_max_p = u, r_max_p = 1)

    # Solve problem with Ipopt.
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "max_iter", 1000)
    JuMP.optimize!(model)

    # Get solution
    local x = JuMP.value.(model[:x]).data
    local y = JuMP.value.(model[:y]).data

    line = plot_trajectory!(ax, ((a, b) -> (a, b)).(x, y))
    push!(l, line)
end

plot_endpoints!(ax, (0, -3, π/2), (3, 0, 0))
CairoMakie.Label(f[0, :], "N = $nh pts\nu min = 0.9 m/2\nu max = 2.7 m/s\nr max = 1 rad/s\nCondition initale en ligne droite", justification = :left, fontsize = 12, halign=:left)
CairoMakie.Legend(f[1, 2], l, ["$(u) m/s" for u in u_max])

f
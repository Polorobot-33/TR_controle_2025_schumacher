
using JuMP

include("models.jl")
include("initial_conditions.jl")
include("renderer.jl")




# Load problem
nh = 100
model, init = robot_rect_model(nh)
# Solve problem with Ipopt.
JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "max_iter", 1000)
JuMP.optimize!(model)
# Get solution


x = JuMP.value.(model[:x]).data
y = JuMP.value.(model[:y]).data
ϕ = JuMP.value.(model[:ϕ]).data
u = JuMP.value.(model[:u]).data
r = JuMP.value.(model[:r]).data
T = JuMP.value.(model[:T])

sol = x, y, ϕ, u, r, T
show_results(sol, [init,]; title = "Solution Optimale avec contraintes géométriquess", metadata = ["JuMP with Ipopt", "u max = 2.7 m/s"])
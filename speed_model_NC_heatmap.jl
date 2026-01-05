using JuMP
using CairoMakie
using LinearAlgebra

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")

# nombre de step de la simulation
nh = 100

# conditions initales (x, y, ϕ, u, r)
cdt_0 = (0, 0, 0, 0, 0)

function getTime(x, y) 
    norm((x, y)) < 0.1 && return 0


    cdt_f = (x, y, 0, 0, 0)
    model = Model()
    speed_model!(model, nh, cdt_0, cdt_f)
    solve!(model)

    return JuMP.value.(model[:T])
end


f = CairoMakie.Figure(size = (512, 512))
ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())
hm = heatmap!(ax, -0.2:0.1:0.2, -0.2:0.1:0.2, getTime, colormap = :deep, interpolate=true)

Label(f[0, :], "Temps de trajet en partant de (0, 0),\n tourné vers la droite,\n jusqu'à la position finale,\n orientée vers la droite", fontsize = 18)

Colorbar(f[1, end+1], hm)

f

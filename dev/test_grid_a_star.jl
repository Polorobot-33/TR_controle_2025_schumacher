using CairoMakie
using PoissonDiskSampling
using Graphs

include("../initial_conditions.jl")
include("../polygon.jl")


# parameters
rangex = (-5, 5)
rangey = (-5, 5)
spacing = 0.2



f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect())


poly = []
center = PoissonDiskSampling.generate(2.8, rangex, rangey)
for c in center
    points = convex_volume_polygon(5, 1.0)
    points = [p .+ c for p in points]

    # collision polytop
    push!(poly, points)
end

px = [x for x in (rangex[1]:spacing:rangex[2])]
py = [y for y in (rangey[1]:spacing:rangey[2])]

@time gridp = near_poly_grid(poly, px, py, spacing)

mesh = Graphs.SimpleGraphs.grid([length(px), length(py)])
mask = Int32[]
for i in 1:length(px), j in 1:length(py)
    gridp[i, j] && push!(mask, i + (j-1)*length(px))
end
vmap = rem_vertices!(mesh, mask)

path_raw = [e.src for e in Graphs.a_star(mesh, 50, 1800)]
push!(path_raw, 1800)
path = vmap[path_raw]

heatmap!(ax, px, py, gridp, colormap=:greens)
for vert in poly 
    poly!(ax, [Tuple(p) for p in vert], strokecolor=:red, strokewidth=1, color=:transparent)
end

Nx = length(px)
Ny = length(py)
pos(i) = ((((i-1) % Nx)) / (Nx-1) * (10) - 5, (((i-1) ÷ Nx)) / (Ny-1) * (10) - 5)
lines!(ax, pos.(path))


f
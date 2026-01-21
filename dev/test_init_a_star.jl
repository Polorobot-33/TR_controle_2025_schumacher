using CairoMakie
using PoissonDiskSampling
using Graphs

include("../initial_conditions.jl")
include("../polygon.jl")

# parameters
rangex = (-5, 5)
rangey = (-5, 5)
spacing = 0.2

Random.seed!(15)

f = Figure()
ax = Axis(f[1, 1], aspect=DataAspect(), xlabel="position x (m)", ylabel="position y (m)")


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

#=heatmap!(ax, px, py, gridp, colormap=:greens)
for vert in poly 
    poly!(ax, [Tuple(p) for p in vert], strokecolor=:red, strokewidth=1, color=:transparent)
end=#

Nx = length(px)
Ny = length(py)
get_pos(i) = ((((i-1) % Nx)) / (Nx-1) * (10) - 5, (((i-1) ÷ Nx)) / (Ny-1) * (10) - 5)
#lines!(ax, get_pos.(path))


# extract the initial conditions from the path
nh = 32
ind2x(i) = map((i-1) % Nx, 0, Nx-1, rangex[1], rangex[2])
ind2y(i) = map((i-1) ÷ Nx, 0, Ny-1, rangey[1], rangey[2])
interp = [map(i, 0, nh, 1, length(path)) for i in 0:nh]
x_i = [lerp(ind2x(path[floor(Int, t)]), ind2x(path[ceil(Int, t)]), t%1) for t in interp]
y_i = [lerp(ind2y(path[floor(Int, t)]), ind2y(path[ceil(Int, t)]), t%1) for t in interp]
#=x_i[begin] = x_0
x_i[end]   = x_f
y_i[begin] = y_0
y_i[end]   = y_f=#
x_i, y_i = smooth(x_i, y_i; s=2)

lines!(ax, x_i, y_i, color=:green)
scatter!(ax, x_i, y_i, color=:red)

scatter!(ax, [get_pos(path[begin]), get_pos(path[end])], marker=:utriangle, color=:black, markersize=15)

f
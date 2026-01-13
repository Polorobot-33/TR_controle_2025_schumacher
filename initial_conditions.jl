using DataInterpolations
using LinearAlgebra
using Graphs

include("polygon.jl")

function init_straight(N, pos_0, pos_f, u)
    distance = norm(pos_f .- pos_0)
    T = distance / u

    interp = LinearInterpolation([pos_0, pos_f], [0, N])
    pos = reduce(hcat, interp(0:N))

    phi = atan(pos_f[2] .- pos_0[2], pos_f[1] .- pos_0[1])

    return pos[1,:], pos[2,:], fill(phi, N+1), fill(u, N+1), fill(0, N+1), T
end

function init_arc(N::Int, y_0::Real, x_f::Real, u::Real, r::Real, radius::Real)
    angle = LinRange(π, π/2, 64)
	p_x = [0.; radius .* (cos.(angle) .+ 1); x_f]
	p_y = [y_0; radius .* (sin.(angle) .- 1); 0.]

    d0 = abs(y_0) - radius
    df = abs(x_f) - radius
    d_tot = d0 + df + 0.5*π*radius
    p_d = [0; ((π .- angle) .* radius .+ d0) .* (N/d_tot); N]
    
    interp_x = LinearInterpolation(p_x, p_d)
    interp_y = LinearInterpolation(p_y, p_d)
    x = interp_x(0:N)
    y = interp_y(0:N)

    u_turn = min(u, radius * r)
    T = (d0 + df) / u + 0.5*π*radius / u_turn

    phi = [[atan(y[i+1]-y[i], x[i+1]-x[i]) for i in 1:N]; 0]
    r = [0; (phi[2:N+1] .- phi[1:N]) ./ (T / N)]
    u = [abs(r) > 0 ? u_turn : u for r in r]

    return x, y, phi, u, r, T
end

function map(x, frommin, frommax, tomin, tomax; clamp=false)
    if frommin == frommax
        throw("invalid map range")
    end
    y = (x - frommin) / (frommax - frommin) * (tomax - tomin) + tomin
    return clamp ? Base.clamp(y, tomin, tomax) : y
end

function lerp(a, b, x)
    return b * x + (1 - x) * a
end

function smooth(x, y; s=2)
    N = length(x)
    x_o = fill(0.0, N)
    y_o = fill(0.0, N)
    for i in 1:N
        count = min(N, i+s) - max(1, i-s) + 1
        x_o[i] = sum(x[max(1,i-s) : min(N,i+s)]) / count
        y_o[i] = sum(y[max(1,i-s) : min(N,i+s)]) / count
    end
    return x_o, y_o
end

function initAstar(nh, poly, rangex, rangey, cdt_0, cdt_f; spacing=0.1, dmin=0.1)
    x_0, y_0, ϕ_0, u_0, r_0  = cdt_0
    x_f, y_f, ϕ_f, u_f, r_f  = cdt_f
    from = (x_0, y_0)
    to = (x_f, y_f)
    
    # initialisation of the grid
    px = [x for x in (rangex[1]:spacing:rangex[2])]
    py = [y for y in (rangey[1]:spacing:rangey[2])]
    Nx = length(px)
    Ny = length(py)

    grid = near_poly_grid(poly, px, py, dmin)

    mesh = Graphs.SimpleGraphs.grid([length(px), length(py)])
    mask = Int32[]
    for i in 1:length(px), j in 1:length(py)
        grid[i, j] && push!(mask, i + (j-1)*length(px))
    end
    vmap = rem_vertices!(mesh, mask)

    # find start and index positions in the graph
    start_x = round(map(from[1], rangex[1], rangex[2], 0, Nx-0.001; clamp=true))
    start_y = round(map(from[2], rangey[1], rangey[2], 0, Ny-0.001; clamp=true))
    ind_start = findfirst(e -> e == start_x + start_y * Nx + 1, vmap)
    
    end_x = round(map(to[1], rangex[1], rangex[2], 0, Nx-0.001))
    end_y = round(map(to[2], rangey[1], rangey[2], 0, Ny-0.001))
    ind_end = findfirst(e -> e == end_x + end_y * Nx + 1, vmap)

    if isnothing(ind_start) || isnothing(ind_end) 
        throw("start or end position is impossible")
    end

    # calculate the path
    path_raw = [e.src for e in Graphs.a_star(mesh, ind_start, ind_end)]
    if length(path_raw) == 0 throw("couldn't find a path") end
    push!(path_raw, ind_end)
    path = vmap[path_raw]

    # extract the initial conditions from the path
    ind2x(i) = map((i-1) % Nx, 0, Nx-1, rangex[1], rangex[2])
    ind2y(i) = map((i-1) ÷ Nx, 0, Ny-1, rangey[1], rangey[2])
    interp = [map(i, 0, nh, 1, length(path)) for i in 0:nh]
    x_i = [lerp(ind2x(path[floor(Int, t)]), ind2x(path[ceil(Int, t)]), t%1) for t in interp]
    y_i = [lerp(ind2y(path[floor(Int, t)]), ind2y(path[ceil(Int, t)]), t%1) for t in interp]
    x_i, y_i = smooth(x_i, y_i; s=12)

    ϕ_i = [atan(y_i[i+1] - y_i[i-1], x_i[i+1] - x_i[i-1]) for i in 2:nh]
    pushfirst!(ϕ_i, ϕ_0)
    push!(ϕ_i, ϕ_f)

    u_avg = max((u_0 + u_f) / 2, 1.0)
    u = fill(u_avg, nh+1)

    totaldist = sum([norm([x_i[i+1]-x_i[i], y_i[i+1]-y_i[i]]) for i in 1:nh])
    T = totaldist / u_avg

    δt = T / nh
    r = [(ϕ_i[i+1] - ϕ_i[i-1]) / (2*δt) for i in 2:nh]
    pushfirst!(r, r_0)
    push!(r, r_f)

    return x_i, y_i, ϕ_i, u, r, T
end
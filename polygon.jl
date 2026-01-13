using SimpleRandom
using Random
using LinearAlgebra

function convex_volume_polygon(n, v)
    # generate random coordinates
    x = rand(n)
    y = rand(n)

    # sort them
    sort!(x)
    sort!(y)

    # extract two chains, one ascending, the other one descending
    chain_choice_x = collect(random_subset(n-2, (n-1) ÷ 2)) .+ 1
    chain_choice_y = collect(random_subset(n-2, (n-1) ÷ 2)) .+ 1
    sort!(chain_choice_x)
    sort!(chain_choice_y)

    forward_x = [x[begin]; [x[i] for i in chain_choice_x]; x[end]]
    forward_y = [y[begin]; [y[i] for i in chain_choice_y]; y[end]]
    backward_x = [x[begin]; [x_ for (i, x_) in enumerate(x) if i>1 && i<n && !(i in chain_choice_x)]; x[end]]
    backward_y = [y[begin]; [y_ for (i, y_) in enumerate(y) if i>1 && i<n && !(i in chain_choice_y)]; y[end]]

    # consider the vectors from one element to the next in each chain
    vec_x = [forward_x[2:end] .- forward_x[1:end-1]; backward_x[1:end-1] .- backward_x[2:end]]
    vec_y = [forward_y[2:end] .- forward_y[1:end-1]; backward_y[1:end-1] .- backward_y[2:end]]
    shuf_y = shuffle(vec_y)

    # randomly pair them
    vec = ((vx, vy) -> [vx, vy]).(vec_x, shuf_y)

    # sort them by angle
    sort!(vec, by = (u->atan(u[2], u[1])))

    # construct the polygon and calculate it's barycentre
    center = [0., 0.]
    points = [[0., 0.]]
    for v in vec 
        push!(points, points[end] .+ v)
        center .+= points[end]
    end
    pop!(points)
    center ./= n

    # center it
    points = [p .- center for p in points]

    # scale i to achieve the desired area
    area = polygon_area(points)
    scale = sqrt(v / area);
    points = [p .* scale for p in points]

    return points
end

function polygon_area(verts)
    area = triangle_area(verts[1:3]...)
    for i in 3:length(verts)-1
        area += triangle_area(verts[1], verts[i:i+1]...)
    end
    return area
end

function triangle_area(a, b, c)
    return abs(cross2(c .- a, b .- a)) * 0.5
end

function cross2(u, v)
    return u[1]*v[2] - u[2]*v[1]
end

function normalize(u)
    return u ./ norm(u)
end

function dir(a, b)
    return normalize(b .- a)
end

function normal(u)
    return [-u[2], u[1]]
end

function dot(u, v)
    return sum(u .* v)
end

function pointsToEdge(p1, p2, p3)
    norm = normal(dir(p1, p2));
    norm .*= dot(dir(p3, p1), norm) >= 0 ? 1 : -1
    d = dot(p1, norm)
    return norm, d
end

function vertices2Edges(poly)
    poly_C = []
    poly_d = []
    N_faces = 0
    for points in poly
        N_faces = length(points)

        edges = [pointsToEdge(points[i], points[(i)%N_faces+1], points[(i+1)%N_faces+1]) for i in eachindex(points)]
        append!(poly_C, (n -> n[1]).(edges))
        append!(poly_d, (n -> n[2]).(edges))
    end

    return poly_C, poly_d, N_faces
end

"""
    returns true if the point p lies inside of the convex polygon defined by A and b
"""
function in_poly(A, b, p)
    return sum(A*p .> b) == 0
end

"""
    returns true if the point p lies inside of the set of N_faces convex polygons defined by A and b
"""
function in_poly(A, b, N_faces, p)
    N = length(A)
    for i in 1:N_faces:N
        A_p = transpose(reduce(hcat, A[i:i+N_faces-1]))
        in_poly(A_p, b[i:i+N_faces-1], p) && return true
    end
    return false
end

function near_single_poly(vert::Vector{Any}, p::Vector{Float64}, dmin::Float64)
    i = 1
    prev_move = 0
    N = length(vert)

    count = 0
    while(count <= N)
        count += 1

        v1 = vert[(i-1) % N + 1]
        v2 = vert[(i) % N + 1]
        v3 = vert[(i+1) % N + 1]
        d = norm(v2 .- v1)

        tan = dir(v1, v2)
        n = normal(tan)
        n .*= dot(dir(v3, v1), n) >= 0 ? 1 : -1

        d_proj = dot(p .- v1, tan)

        side = (d_proj > d) - (d_proj < 0)
        proj = v1 .+ clamp(d_proj, 0, d) .* dir(v1, v2)


        if dot(p, n) <= dot(v1, n)
            if prev_move != 0 return (norm(proj .- p) < dmin) end
            i += 1
            continue
        end

        (side * prev_move < 0 || side == 0) && return (norm(proj .- p) < dmin)
        #print(side)
        
        prev_move = side
        i = (i + side + length(vert) - 1) % length(vert) + 1
    end

    # the point is already inside the polygon
    #println(p)
    return true
end

function near_poly(poly, p, dmin)
    #println(p)
    for vert in poly
        near_single_poly(vert, p, dmin) && return true
    end
    return false
end

function near_single_poly_help(N, vert, p, dmin::Float64, n, t, d)
    i = 1
    prev_move = 0
    
    count = 0
    while(count <= N)
        count += 1

        v = vert[i]
        d_proj = dot(p .- v, t[i])

        side = (d_proj > d[i]) - (d_proj < 0)
        proj = v .+ clamp(d_proj, 0, d[i]) .* t[i]


        if dot(p, n[i]) <= dot(v, n[i])
            if prev_move != 0 return (norm(proj .- p) < dmin) end
            i = (i + 1 + N - 1) % N + 1
            continue
        end

        (side * prev_move < 0 || side == 0) && return (norm(proj .- p) < dmin)
        
        prev_move = side
        i = (i + side + N - 1) % N + 1
    end

    # the point is already inside the polygon
    return true
end

function near_poly_grid(poly, px, py, dmin)
    Nx = length(px)
    Ny = length(py)
    grid = fill(false, (Nx, Ny))

    for vert in poly
        N = length(vert)

        n = fill([0., 0.], N)
        t = fill([0., 0.], N)
        d = fill(0., N)
        bounds = [vert[begin][1], vert[begin][1], vert[begin][2], vert[begin][2]]
        for (i, v1) in enumerate(vert)
            v2 = vert[(i) % N + 1]
            v3 = vert[(i+1) % N + 1]

            d[i] = norm(v2 .- v1)
            t[i] = dir(v1, v2)
            n[i] = normal(t[i])
            n[i] .*= dot(dir(v3, v1), n[i]) >= 0 ? 1 : -1

            bounds[1] = min(bounds[1], v1[1])
            bounds[2] = max(bounds[2], v1[1])
            bounds[3] = min(bounds[3], v1[2])
            bounds[4] = max(bounds[4], v1[2])
        end

        beginx, beginy, endx, endy = 0, 0, 0, 0
        for (i, x) in enumerate(px) if (x >= bounds[1] - dmin) beginx=i; break end end
        for (i, x) in Iterators.reverse(enumerate(px)) if (x <= bounds[2] + dmin) endx=i; break end end
        for (i, y) in enumerate(py) if (y >= bounds[3] - dmin) beginy=i; break end end
        for (i, y) in Iterators.reverse(enumerate(py)) if (y <= bounds[4] + dmin) endy=i; break end end

        for ix in max(beginx,1):min(endx, Nx), iy in max(beginy,1):min(endy, Ny)
            near_single_poly_help(N, vert, [px[ix], py[iy]], dmin, n, t, d) && (grid[ix, iy] = true)
        end
    end

    return grid
end



;
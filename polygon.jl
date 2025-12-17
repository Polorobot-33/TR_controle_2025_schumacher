using SimpleRandom
using Random

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

;
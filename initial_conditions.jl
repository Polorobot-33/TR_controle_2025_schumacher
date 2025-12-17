using DataInterpolations
using LinearAlgebra

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
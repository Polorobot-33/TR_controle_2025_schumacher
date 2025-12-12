using DataInterpolations
using LinearAlgebra

function init_straight(N, pos_0, pos_f, u)
    distance = norm(pos_f .- pos_0)
    T = distance / u

    interp = LinearInterpolation([pos_0, pos_f], [0, N])
    pos = interp(0:N)

    phi = atan(pos_f[2] - pos_0[2], pos_f[2] - pos_0[1])


    return pos[0], pos[1], fill(phi, N+1), fill(u, N+1), fill(0, N+1), T
end

function init_arc(N::Int, y_0::Real, x_f::Real, u::Real, radius::Real)
    angle = LinRange(π, π/2, 64)
	p_x = [0; radius .* (cos.(angle) .+ 1); x_f]
	p_y = [y_0; radius .* (sin.(angle) .- 1); 0]

    d0 = abs(y_0) - radius
    df = abs(x_f) - radius
    d_tot = d0 + df + 0.5*π*radius
    p_d = [0; (π .- angle) .* radius .+ d0; d_tot]
    p_d .*= N / d_tot
    
    interp_x = LinearInterpolation(p_x, p_d)
    interp_y = LinearInterpolation(p_y, p_d)
    x = interp_x(0:N)
    y = interp_y(0:N)

    T = (d0 + df + 0.5*π*radius) / u

    phi = [π/2; [atan(y[i+1]-y[i-1], x[i+1]-x[i-1]) for i in 2:N]; 0]
    r = [(phi[2:N+1] - phi[1:N]) * T / (N+1); 0]

    return x, y, phi, fill(u, N+1), r, T
end
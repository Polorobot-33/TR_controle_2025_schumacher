using ADTypes
using FiniteDiff
using Ipopt
using LinearAlgebra
using Optimization
using OptimizationMOI
using Zygote
using CairoMakie
using Random



function dyn(x::Vector, u::Real, r::Real)
	_, _, ϕ = x
	return [u * cos(ϕ), u * sin(ϕ), r]
end

function dynamics_cons(x::Vector, u::Vector, r::Vector, T::Real, N::Int)
	δt = T / N
	arr = ((xnp1, xn, xdot) -> xnp1 .- (xn .+ δt*T.*xdot)).(x[2:end], x[1:end-1], dyn.(x[1:end-1], u[1:end-1], r[1:end-1]))
	return reduce(vcat, arr)#collect(Iterators.flatten(arr))
end

function rotate(p, ϕ)
	return (cos(ϕ) * p[1] - sin(ϕ) * p[2], sin(ϕ) * p[1] + cos(ϕ) * p[2])
end

"""
	returns the coordinates of the ith corner of a centered rectangle rotated by ϕ
"""
function corner(ϕ, width, length, i)
	cx = (i ÷ 2) % 2 == 0 ? width/2 : -width/2
	cy = ((i+1) ÷ 2) % 2 == 0 ? length/2 : -length/2
	
	return rotate([cx, cy], ϕ)
end

function smooth_minmax(a, b, alpha)
    return (a * exp(a*alpha) + b * exp(b*alpha)) / (exp(a*alpha) + exp(b*alpha))
end

function corner_cons(ϕ, x, y, l1, l2, rw, rh, alpha, i)
	c = corner(ϕ, rw, rh, i) .+ [x, y]
	return [abs(c[1]) / (l1/2), abs(c[2]) / (l2/2)]
	#return 1 - smooth_minmax(abs(c[1]) / (l1/2), abs(c[2]) / (l2/2), -alpha)
end

function robot_corner_cons(ϕ::Real, x::Real, y::Real, l1::Real, l2::Real, rw::Real, rh::Real, alpha::Float64, i::Int)
	c = corner(0, l1, l2, i)

	# calculate the coordinates of the corner c in the base of the robot
	p = rotate(c .- (x, y), -ϕ)
	
	#return smooth_minmax(abs(p[1]) / (rw/2), abs(p[2]) / (rh/2), alpha) - 1
	return [abs(p[1]) / (rw/2), abs(p[2]) / (rh/2)]
end

function single_geom_cons(x::Vector, dimensions::Tuple, alpha::Float64)
	x_, y_, ϕ_ = x
	rw, ry, l1, l2 = dimensions
	


	return [[corner_cons(ϕ_, x_, y_, l1, l2, rw, ry, alpha, i) for i in 1:4]; [robot_corner_cons(ϕ_, x_, y_, l1, l2, rw, ry, alpha, t) for t in 1:4]]
end

function geom_cons(x::Vector, dimensions::Tuple, alpha::Float64)
	cons = (s -> single_geom_cons(s, dimensions, alpha)).(x)
	return reduce(vcat, cons)#collect(Iterators.flatten(cons))
end

function complement_cons(comp::Vector, x::Vector, epsilon::Float64)
	x_plus = comp[1:2:end]
	x_minus = comp[2:2:end]

	return [x_plus .- x_minus .- x; x_plus; x_minus; epsilon .- (x_plus .* x_minus)]
end

function problem_cons(res, x, p::Tuple)
	# unpack the parameters
	N, dimensions, U_max, cdt_0, cdt_f, epsilon = p
	
	# unpack the optimization vector
	x_ = [[x[i], x[i+1], x[i+2]] for i in 1:5:5*N]
	u_ = x[4:5:5*N]
	r  = x[5:5:5*N]
	comp = x[5*N+1:end-1]
	T  = x[end]

	# calculate geometric collision variables
	geom = geom_variables(x_, dimensions)

	# calculate constraints
	speed_cons = [(speed -> U_max[1] - abs(speed)).(u_) ; (turn -> U_max[2] - abs(turn)).(r)]
	end_cons = [[x_[begin]; u_[begin]; r[begin]] .- cdt_0; [x_[end]; u_[end]; r[end]] .- cdt_f]
	
	res .= [dynamics_cons(x_, u_, r, T, N); end_cons; complement_cons(comp, epsilon); geom_cons(x_, dimensions, alpha); speed_cons; x[end]]
end

function lerp(a, b, x)
	return b .* x .+ (1-x) .* a
end

function initialConditionsStraight(N, U_max, x_0, x_f)
	p0 = [x_0...][1:2]
	pf = [x_f...][1:2]
	ϕ = atan(pf[2] - p0[2], pf[1] - p0[1])

	x = [[lerp(p0, pf, (i-1)/(N-1)); ϕ] for i in 1:N]
	u = fill(U_max[1] / 2, N)
	r = zeros(N)
	
	distance = norm(pf .- p0)
	T = distance / (U_max[1] / 2)

	return x, u, r, T
end

function lerp_path(x, y, f)
    for (i, px) in enumerate(x)
        px > f && (return lerp(y[i-1], y[i], (f - x[i-1]) / (px - x[i-1])))
    end
    return x[end]
end

"""
    initialize a path that starts at p0, passes through the center to reach pf 
"""
function initial_conditions_rect(N, U_max, cdt_0, cdt_f)
	p0 = [cdt_0...][1:2]
	pf = [cdt_f...][1:2]

    lA = norm(p0)
    lB = norm(pf)
    d_tot = lA + lB
    T = d_tot / (U_max[1] / 2)
    ϕA = atan(-p0[2], -p0[1])
    ϕB = atan(pf[2], pf[1])

	x = [[lerp_path([0, lA/d_tot, d_tot], [p0, [0, 0], pf], (i-1)/(N-1)); i/(N) < (lA/d_tot) ? ϕA : ϕB] for i in 1:N]
    u = [[norm(x[i+1]-x[1]) / (T/(N-1)) for i in 1:N-1]; cdt_f[end]]
    r = zeros(N)

    return x, u, r, T
end

function objective(u::Vector, p::Tuple)
	return u[end]
end

function solve_problem(x0::Vector{Float64}, U_max::Tuple{Float64, Float64}, cdt_0::NTuple{5, Float64}, cdt_f::NTuple{5, Float64}, dimensions::NTuple{4, Float64})
	# paramètres : N, U_max, x_0, x_f, dimensions géométriques
	N = (length(x0) - 1) ÷ (5 + 16)
	params = (N, dimensions, U_max, cdt_0, cdt_f, 3.)
	
	# initialisation des contraintes
	# dynamiques : 3 * (N-1)
	# conditions initiales et finales : 10
	# constraintes de complémentarité : 32 * N (8 égalités, 24 inégalités)
	# vitesses maximales : 2 * N
	# géométriques : 8 * N
	# temps : 1
	lcons = zeros(3 * (N-1) + 10 + 32*N + 2*N + 8*N + 1)
	ucons = [zeros(3 * (N-1) + 10 + 8*N); fill(Inf,24*N + 2*N + 8*N + 1)]

	#ADTypes.AutoZygote()
	optim_f = OptimizationFunction(objective, ADTypes.AutoZygote(), cons = problem_cons)
	prob = Optimization.OptimizationProblem(optim_f, x0, params, lcons=lcons, ucons=ucons)
	sol = solve(prob, Ipopt.Optimizer())

	return sol.u
end

function plot_terrain!(ax, width, height)
  points = [(-3, height/2),  (-width/2, height/2),  (-width/2, 3),
  			(3, height/2),   (width/2, height/2),   (width/2, 3),
  			(-3, -height/2), (-width/2, -height/2), (-width/2, -3),
  			(3, -height/2),  (width/2, -height/2),  (width/2, -3)]
  linesegments!(ax, points[[1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11, 12]], color=:blue)
end

function plot_trajectory!(ax, data)
	points = [(x, y) for (x, y, _) in data]
	lines!(ax, points, color=:red)
end

"""
	plots a reference trajectory
"""
function plot_circle!(ax, robot_length, width)

	angle = linRange(-π/2, 0, 36)
	r = (width - robot_length) / 2
	x = [0; r .* (cos.(angle .- π/2) .+ 1); 3]
	y = [-3; r .* (sin.(angle .- π/2) .- 1); 0]

	lines!(ax, x, y, color=:yellow)
end

function plot_endpoints!(ax, x0, xf)
	scatter!(ax, [x0[1:2], xf[1:2]], marker=:utriangle, color=:black, markersize=15)
end

function plot_robot!(ax, pos, size)
	x, y, ϕ = pos
	w, l = size
	
	points = [(x, y) .+ Tuple(corner(ϕ, w, l, i)) for i in 1:5]
	lines!(ax, points, color=:green)
end

function plot_positions!(ax, pos, size, spacing)
    N = length(pos)
    for i in 1:spacing:N
		plot_robot!(ax, pos[i], size)
	end
end

function test(N)
	U_max = (2.4, 1.5)
	cdt_0 = (0., -2., π/2, 2., 0.)
	cdt_f(2., 0., 0., 2., 0.)

	x, u, r, T = initial_conditions_rect(N, U_max, cdt_0, cdt_f)
	x0 = Float64[]
	# append variables
	for elem in zip(x, u, r)
		append!(x0, [elem[1]...])
		push!(x0, elem[2], elem[3])
	end
	# append complementary constraints variables
	append!(x0, zeros(8*2*N))
	push!(x0, T)

    return solve_problem(x0, U_max, cdt_0, cdt_f, (1.128, 0.720, 2.4, 0.9))
end

#= function test_sub(sol, subsampling)
	N = (length(sol) - 1) ÷ 5
	
	# initialisation de x0, en interpolant linéairement la solution précédente
	x0 = Float64[]
	for i in 0:N-2
		for j in 1:subsampling
			interpolated = lerp(sol[5*i+1 : 5*i+3], sol[5*(i+1)+1 : 5*(i+1)+3], (j-1)/(subsampling))
			append!(x0, [interpolated; sol[5*i+4]; j == 1 ? sol[5*i+5]*(subsampling+1) : 0])
		end
	end
	append!(x0, sol[5*N-4:5*N])
	push!(x0, sol[end])

	return solve_problem(x0, (2.4, 1.5), (0., -2., π/2, 2., 0.), (2., 0., 0., 2., 0.), (1.128, 0.720, 2.4, 0.9))
end =#

function show_results(in)
	N = (length(in) - 1) ÷ (5 + 16)
	sol = reshape(in[1:5*N], (5, N))
	pos = (s -> Tuple(s)).(eachcol(sol[1:3, :]))
	u = sol[4, :]
	r = sol[5, :]
	T = in[end]

    f = Figure(size = (512, 820))
    times = LinRange(0, T, N)

    Label(f[-1, :], "Solution optimale avec contraintes géométriques", fontsize = 18)
    Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $N pts", justification = :left, fontsize = 12, halign=:left)

    ax = Axis(f[1, 1], aspect = DataAspect(), alignmode=Inside())
    plot_terrain!(ax, 2.4, 0.9)
	plot_circle!(ax, 1.172, 2.4)

    plot_positions!(ax, pos, (1.128, 0.720), 1)
    plot_trajectory!(ax, pos)
    plot_endpoints!(ax, (0., -2., π/2), (2., 0., 0.))

    ax_u = Axis(f[2, :], xlabel="temps", ylabel="commande u")
    ax_r = Axis(f[3, :], xlabel="temps", ylabel="commande r")
    linkxaxes!(ax_u, ax_r)

    lines!(ax_u, times, u)
    lines!(ax_r, times, r)

    rowsize!(f.layout, 1, Aspect(1, 1.0))

	return f
end
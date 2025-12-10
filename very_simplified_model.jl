using ADTypes
using FiniteDiff
using Ipopt
using LinearAlgebra
using Optimization
using OptimizationMOI
using Zygote
using CairoMakie
using Random


"""
	expression de la dérivée dans le modèle de dynamique : x' = f(x, u)
"""
function dyn(x::Vector, u::Real, r::Real)
	_, _, ϕ = x
	return [u * cos(ϕ), u * sin(ϕ), r]
end

"""
	contraintes associées à la dynamique (schéma d'euler explicite)
	Remarque : le schéma d'euler explicite est beaucoup plus rapide que le schéma d'euler implicite
"""
function dynamics_cons(x::Vector, u::Vector, r::Vector, T::Real, N::Int)
	δt = T / (N-1)
	arr = ((xnp1, xn, xdot) -> xnp1 .- (xn .+ δt*T.*xdot)).(x[2:end], x[1:end-1], dyn.(x[1:end-1], u[1:end-1], r[1:end-1]))
	return reduce(vcat, arr)#collect(Iterators.flatten(arr))
end

function rotate(p, ϕ)
	return (cos(ϕ) * p[1] - sin(ϕ) * p[2], sin(ϕ) * p[1] + cos(ϕ) * p[2])
end

"""
	renvoie les coordonnées du sommet i d'un rectangle de dimensions (width, length) tourné de ϕ
"""
function corner(ϕ, width, length, i)
	cx = (i ÷ 2) % 2 == 0 ? width/2 : -width/2
	cy = ((i+1) ÷ 2) % 2 == 0 ? length/2 : -length/2
	return rotate([cx, cy], ϕ)
end


function corner(width, length, i)
    cx = (i ÷ 2) % 2 == 0 ? width/2 : -width/2
	cy = ((i+1) ÷ 2) % 2 == 0 ? length/2 : -length/2
	return [cx, cy]
end


"""
	apprixomation régularisante de la fonction max (ou min en donnant un alpha < 0)
"""
function smooth_minmax(a, b, alpha)
    return (a * exp(a*alpha) + b * exp(b*alpha)) / (exp(a*alpha) + exp(b*alpha))
end

"""
	contrainte géométrique pour exprimer l'appartenance d'un point au domaine admissible (croix)
	Cette grandeur doit être ≥ 0
"""
function single_geom_cons(x::Vector, dimensions::Tuple, alpha::Float64)
	x_, y_, _ = x
	l1, l2 = dimensions
	
    return 1 - smooth_minmax(abs(x_) / (l1/2), abs(y_) / (l2/2), -alpha)
end

"""
	obtention d'un vecteur qui concentre toutes les contraintes géometriques des points de la trajectoire
"""
function geom_cons(x::Vector, dimensions::Tuple, alpha::Float64)
    return (s -> single_geom_cons(s, dimensions, alpha)).(x)
end

"""
	fonction contraintes utilisée dans le problème d'optimisation
"""
function problem_cons(res, x, p::Tuple)
	# unpack the parameters
	N, dimensions, U_max, cdt_0, cdt_f, alpha = p
	
	# unpack the optimization vector
	x_ = [[x[i], x[i+1], x[i+2]] for i in 1:5:5*N]
	u_ = x[4:5:5*N]
	r  = x[5:5:5*N]
	T  = x[end]

	# calculate constraints
	speed_cons = [(speed -> U_max[1] - abs(speed)).(u_) ; (turn -> U_max[2] - abs(turn)).(r)]
	end_cons = [[x_[begin]; u_[begin]; r[begin]] .- cdt_0; [x_[end]; u_[end]; r[end]] .- cdt_f]
	
	res .= [dynamics_cons(x_, u_, r, T, N); end_cons; geom_cons(x_, dimensions, alpha); speed_cons; x[end]]
end


function lerp(a, b, x)
	return b .* x .+ (1-x) .* a
end

"""
	génère un vecteur initial x0 suivant une trajectoire en ligne droite entre les points de départ et d'arrivée
"""
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

"""
	interpolation linéaire entre plusieurs vecteurs formant un chemin
"""
function lerp_path(x, y, f)
    for (i, px) in enumerate(x)
        px > f && (return lerp(y[i-1], y[i], (f - x[i-1]) / (px - x[i-1])))
    end
    return x[end]
end

"""
    génération d'un vectuer initial x0 partant de p0 et allant jusqu'à pf en passant par l'origine (0, 0)
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

"""
	fonction à optimiser (ici T)
"""
function objective(u::Vector, p::Tuple)
	return u[end]
end

function solve_problem(x0::Vector{Float64}, U_max::Tuple{Float64, Float64}, cdt_0::NTuple{5, Float64}, cdt_f::NTuple{5, Float64}, dimensions::NTuple{4, Float64})
	# paramètres : N, U_max, x_0, x_f, dimensions géométriques
	N = (length(x0) - 1) ÷ 5
	params = (N, dimensions, U_max, cdt_0, cdt_f, 3.)
	
	# initialisation des contraintes
	# dynamiques : 3 * (N-1)
	# conditions initiales et finales : 10
	# vitesses maximales : 2 * N
	# géométriques : N
	# temps : 1
	lcons = zeros(3 * (N-1) + 10 + 2*N + N + 1)
	ucons = [zeros(3 * (N-1) + 10); fill(Inf,2*N + N + 1)]

	#ADTypes.AutoZygote()
	optim_f = OptimizationFunction(objective, ADTypes.AutoZygote(), cons = problem_cons)
	prob = Optimization.OptimizationProblem(optim_f, x0, params, lcons=lcons, ucons=ucons)
	sol = solve(prob, Ipopt.Optimizer())

	return sol.u
end

"""
	affiche les obstacles
"""
function plot_terrain!(ax, width, height)
  	points = [(-2, height/2),  (-width/2, height/2),  (-width/2, 2),
  			(4, height/2),   (width/2, height/2),   (width/2, 2),
  			(-2, -height/2), (-width/2, -height/2), (-width/2, -4),
  			(4, -height/2),  (width/2, -height/2),  (width/2, -4)]
  	linesegments!(ax, points[[1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11, 12]], color=:blue)
end

"""
	affiche la trajectoire
"""
function plot_trajectory!(ax, data)
	points = [(x, y) for (x, y, _) in data]
	lines!(ax, points, color=:red)
end

"""
	affiche une trajectoire de référence
"""
function plot_circle!(ax, robot_length, width)

	angle = LinRange(π, π/2, 36)
	r = (width - robot_length) / 2
	x = [0; r .* (cos.(angle) .+ 1); 3]
	y = [-3; r .* (sin.(angle) .- 1); 0]

	lines!(ax, x, y, color=:magenta)
end

"""
	affichage du départ et de l'arrivée
"""
function plot_endpoints!(ax, x0, xf)
	scatter!(ax, [x0[1:2], xf[1:2]], marker=:utriangle, color=:black, markersize=15)
end

"""
	déssine un rectangle représentant le robot
"""
function plot_robot!(ax, pos, size)
	x, y, ϕ = pos
	w, l = size
	
	points = [(x, y) .+ Tuple(corner(ϕ, w, l, i)) for i in 1:5]
	lines!(ax, points, color=:green)
end

"""
	déssine les positions successives du robot
"""
function plot_positions!(ax, pos, size, spacing)
    N = length(pos)
    for i in 1:spacing:N
		plot_robot!(ax, pos[i], size)
	end
end

"""
	créé un problème à résoudre avec des valeurs prédéfinies, et le résout pour un nombre N de points dans la trajectoire
"""
function test(N)
	U_max = (2.4, 1.5)
	cdt_0 = (0., -3., π/2, 2.0, 0.)
	cdt_f = (3., 0., 0., 2.0, 0.)

	x, u, r, T = initial_conditions_rect(N, U_max, cdt_0, cdt_f)
	x0 = Float64[]
	# append variables
	for elem in zip(x, u, r)
		append!(x0, [elem[1]...])
		push!(x0, elem[2], elem[3])
	end
	push!(x0, T)

    return solve_problem(x0, U_max, cdt_0, cdt_f, (1.128, 0.720, 2.4, 0.9))
end

"""
	affiche le résultat de la fonction test
"""
function show_results(in)
	N = (length(in) - 1) ÷ 5
	sol = reshape(in[1:5*N], (5, N))
	pos = (s -> Tuple(s)).(eachcol(sol[1:3, :]))
	u = sol[4, :]
	r = sol[5, :]
	T = in[end]

    f = Figure(size = (512, 860))
    times = LinRange(0, T, N)

    Label(f[-1, :], "Solution optimale avec contraintes géométriques,\n modélisation par un point", fontsize = 18)
    Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $N pts", justification = :left, fontsize = 12, halign=:left)

    ax = Axis(f[1, 1], aspect = DataAspect(), alignmode=Inside())
    plot_terrain!(ax, 2.4, 0.9)
	plot_circle!(ax, 1.172, 2.4)

    plot_positions!(ax, pos, (1.128, 0.720), 1)
    plot_trajectory!(ax, pos)
    plot_endpoints!(ax, (0., -3., π/2), (3., 0., 0.))

    ax_u = Axis(f[2, :], xlabel="temps", ylabel="commande u")
    ax_r = Axis(f[3, :], xlabel="temps", ylabel="commande r")
    linkxaxes!(ax_u, ax_r)

    lines!(ax_u, times, u)
    lines!(ax_r, times, r)

    rowsize!(f.layout, 1, Aspect(1, 1.0))

	
	return f
end
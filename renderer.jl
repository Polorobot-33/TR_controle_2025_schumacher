using CairoMakie

"""
	affiche les obstacles
"""
function plot_terrain!(ax, width, height)
  	points = [(-2, height/2),  (-width/2, height/2),  (-width/2, 2),
  			(4, height/2),   (width/2, height/2),   (width/2, 2),
  			(-2, -height/2), (-width/2, -height/2), (-width/2, -4),
  			(4, -height/2),  (width/2, -height/2),  (width/2, -4)]
  	CairoMakie.linesegments!(ax, points[[1, 2, 2, 3, 4, 5, 5, 6, 7, 8, 8, 9, 10, 11, 11, 12]], color=:blue)
end

"""
	affiche la trajectoire
"""
function plot_trajectory!(ax, points; col=:red)
	CairoMakie.lines!(ax, points, color=col)
end

"""
	affiche une trajectoire de référence
"""
function plot_ref!(ax, r)

	angle = LinRange(π, π/2, 36)
	x = [0; r .* (cos.(angle) .+ 1); 3]
	y = [-3; r .* (sin.(angle) .- 1); 0]

	CairoMakie.lines!(ax, x, y, color=:magenta)
end

"""
	affichage du départ et de l'arrivée
"""
function plot_endpoints!(ax, x0, xf)
	CairoMakie.scatter!(ax, [x0[1:2], xf[1:2]], marker=:utriangle, color=:black, markersize=15)
end

"""
	déssine un rectangle représentant le robot
"""
function plot_robot!(ax, pos, size)
	x, y, ϕ = pos
	w, l = size
	
	points = [(x, y) .+ Tuple(corner(ϕ, w, l, i)) for i in 1:5]
	CairoMakie.lines!(ax, points, color=:green)
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


function show_results(sol, ref; title="Solution optimale", metadata="")
     x, y, ϕ, u, r, T = sol

    N = length(x)
	pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

    f = CairoMakie.Figure(size = (512, 900))
    times = LinRange(0, T, N)

    meta = "\n" * join(metadata, "\n")
    CairoMakie.Label(f[-1, :], title, fontsize = 18)
    CairoMakie.Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\nN = $N pts" * meta, justification = :left, fontsize = 12, halign=:left)

    ax = CairoMakie.Axis(f[1, 1], aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())
    plot_terrain!(ax, 2.4, 0.9)

    plot_positions!(ax, pos, (1.128, 0.720), 1)
    (traj -> plot_trajectory!(ax, ((a, b) -> (a, b)).(traj[1], traj[2]); col=:magenta)).(ref)
    plot_trajectory!(ax, [(px, py) for (px, py, _) in pos])
    plot_endpoints!(ax, (0., -3., π/2), (3., 0., 0.))

    ax_u = CairoMakie.Axis(f[2, :], xlabel="temps", ylabel="commande u")
    ax_r = CairoMakie.Axis(f[3, :], xlabel="temps", ylabel="commande r")
    CairoMakie.linkxaxes!(ax_u, ax_r)

    CairoMakie.lines!(ax_u, times, u)
    CairoMakie.lines!(ax_r, times, r)

    CairoMakie.rowsize!(f.layout, 1, CairoMakie.Aspect(1, 1.0))

	return f
end
using JuMP
using CairoMakie

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")

#= CHANGES
#  
#  format calc_time (2 justificatifs) + " s"
#
=#





collocation_p = [("implicit_Euler", ImplicitEuler()),
                ("Runge_Kutta", RungeKutta()),
                ("Crank_Nicolson", CrankNicolson()), 
                ("RadauIIA-2", RadauIIA(2))]
model_p = [("cinematic", speed_model!),
           ("dynamic",   dynamic_model!)]
nh_p = [50, 100]

nb_tests = length(collocation_p) * length(model_p) * length(nh_p) * 4 * 4
test_count = 0

function save_result(model, nh, case, model_name, collocation, collision_name, terrain, start)
    println(test_count)
    println("completion : $(trunc(test_count / nb_tests * 100, digits=1, base=10))")
    
    f = Figure(size = (600, 1024))
    ax = CairoMakie.Axis(f[1, 1], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

    T = JuMP.value(model[:T])
    nb_iter = JuMP.barrier_iterations(model)
    calc_time = JuMP.solve_time(model)
    times = LinRange(0, T, nh+1)

    Label(f[-1, :], "Solution pour \"$(case)\", comparaison des différents modèles", fontsize = 18)
    Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\n
                    modèle : $(model_name)\n
                    collisions : $(collision_name)\n                
                    collocation : $(collocation)\n
                    nb iterations : $(nb_iter)\n
                    tps de calcul : $(trunc(calc_time, digits=1, base=10))\n
                    N = $(nh+1) pts\n
                    Condition initiale avec A*", justification = :left, fontsize = 12, halign=:left)

    x = var(model, :x)
    y = var(model, :y)
    ϕ = var(model, :ϕ)
    pos = ((a, b, c) -> (a, b, c)).(x, y, ϕ)

    terrain(ax)
    plot_positions!(ax, pos, (1.128, 0.720), 1)
    plot_start!(ax, start)
    plot_trajectory!(ax, [(px, py) for (px, py, _) in pos], col=:red)
    plot_endpoints!(ax, cdt_0[1:3], cdt_f[1:3])


    # commandes
    if :ul ∈ keys(object_dictionary(model)) && :ur ∈ keys(object_dictionary(model))
        ax_u = Axis(f[2, :], xlabel="temps", ylabel="tensions u (V)")
        lines!(ax_u, times, var(model, :ul), label="ul")
        lines!(ax_u, times, var(model, :ur), label="ur")
        axislegend(ax)
    elseif :u ∈ keys(object_dictionary(model)) && :r ∈ keys(object_dictionary(model))
        ax_u = Axis(f[2, :], xlabel="temps", ylabel="vitesse u (m/s)")
        lines!(ax_u, times, var(model, :u))
        ax_r = Axis(f[3, :], xlabel="temps", ylabel="rotation r (rad/s)")
        lines!(ax_r, times, var(model, :r))
        linkxaxes!(ax_u, ax_r)
    end

    rowsize!(f.layout, 1, Aspect(1, 1.0))

    save("comp_figures/$(case)_$(model_name)_$(collision_name)_$(collocation)_$(nh).svg", f)
end






#======== corner ========#


# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (0, -2.5, π/2, 2.7, 0)
cdt_f = (2.5,  0, 0,   2.7, 0)

# collisions
l1_c = 2.4 # largeur principale des couloirs
l2_c = 0.9 # largeur des couloirs latéraux
C = [0 -1; -1 0; 1 0; 0 -1; 0 1; 1 0; -1 0; 0 1]
d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
poly = [[p .*[dx;dy] for p in [[l1_c/2;l2_c/2], [150;l2_c/2], [150;150], [l1_c/2;150]]] for (dx, dy) in [(1,1), (1,-1), (-1,-1), (-1,1)]]


for nh in nh_p
    start = initAstar(nh, poly, (-3, 5), (-5, 3), cdt_0, cdt_f; spacing=0.18, dmin=0.3, smoothing=12)

    for (rk_name, rk) in collocation_p
        for (model_name, m!) in model_p
            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2012_collisions!(model, nh, (4, 2), C, d)
            solve!(model, max_iter=1000)
            save_result(model, nh, "corner", model_name, rk_name, "2012", ax -> plot_terrain!(ax, l1_c, l2_c), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d)
            solve!(model, max_iter=1000)
            save_result(model, nh, "corner", model_name, rk_name, "2017", ax -> plot_terrain!(ax, l1_c, l2_c), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_penetration_collisions!(model, nh, (4, 2), C, d; kappa=1e+3)
            solve!(model, max_iter=1000)
            save_result(model, nh, "corner", model_name, rk_name, "2017 penetration", ax -> plot_terrain!(ax, l1_c, l2_c), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.02)
            solve!(model, max_iter=1000)
            save_result(model, nh, "corner", model_name, rk_name, "2023", ax -> plot_terrain!(ax, l1_c, l2_c), start)
        end
    end
end
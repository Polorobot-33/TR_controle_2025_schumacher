using JuMP
using CairoMakie
using DataFrames
using CSV
using Random
using LinearAlgebra
using PoissonDiskSampling

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")
include("polygon.jl")

#= CHANGES
#  
#  
#
=#


results = DataFrame(case_name=String[], model_name=[], collision=[], collocation=String[], nh=Int[], Time=Float64[], nb_iter=Int[], calc_time=Float64[])



collocation_p = [("implicit_Euler", ImplicitEuler()),
                ("Runge_Kutta", RungeKutta()),
                ("Crank_Nicolson", CrankNicolson()), 
                ("RadauIIA-2", RadauIIA(2))]
model_p = [("cinematic", speed_model!),
           ("dynamic",   dynamic_model!)]
nh_p = [50, 100]

nb_tests = length(collocation_p) * length(model_p) * length(nh_p) * (3+5) * 4
test_count = 0

function save_result(model, nh, case, model_name, collocation, collision_name, terrain, start)
    global test_count += 1
    println("completion : $(trunc(test_count / nb_tests * 100, digits=1, base=10)) %")
    
    f = Figure(size = (600, 1024))
    ax = CairoMakie.Axis(f[1, 1], xlabel="position x (m)", ylabel="position y (m)", aspect = CairoMakie.DataAspect(), alignmode=CairoMakie.Inside())

    T = JuMP.value(model[:T])
    nb_iter = JuMP.barrier_iterations(model)
    calc_time = JuMP.solve_time(model)
    times = LinRange(0, T, nh+1)

    push!(results, (case, model_name, collision_name, collocation, nh+1, T, nb_iter, calc_time))


    Label(f[-1, :], "Solution pour \"$(case)\", comparaison des différents modèles", fontsize = 18, lineheight=0.7)
    Label(f[0, :], "temps de trajet : T = $(trunc(T, digits=3, base=10)) s\n
                    modèle : $(model_name)\n
                    collisions : $(collision_name)\n                
                    collocation : $(collocation)\n
                    nb iterations : $(nb_iter)\n
                    tps de calcul : $(trunc(calc_time, digits=1, base=10)) s\n
                    N = $(nh+1) pts\n
                    Condition initiale avec A*", justification = :left, fontsize = 12, halign=:left, lineheight=0.7)

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
    start = initAstar(nh, poly, (-3, 5), (-5, 3), cdt_0, cdt_f; spacing=0.18, dmin=0.3, smoothing=nh÷8)

    for (rk_name, rk) in collocation_p
        for (model_name, m!) in model_p
            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2012_collisions!(model, nh, (4, 2), C, d)
            solve!(model; max_iter=1000, verbose=0)
            save_result(model, nh, "corner", model_name, rk_name, "2012", ax -> plot_terrain!(ax, l1_c, l2_c), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d)
            solve!(model; max_iter=1000, verbose=0)
            save_result(model, nh, "corner", model_name, rk_name, "2017", ax -> plot_terrain!(ax, l1_c, l2_c), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_penetration_collisions!(model, nh, (4, 2), C, d; kappa=1e+3)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "corner", model_name, rk_name, "2017 penetration", ax -> plot_terrain!(ax, l1_c, l2_c), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "corner", model_name, rk_name, "2023", ax -> plot_terrain!(ax, l1_c, l2_c), start)
        end
    end
end







#======== S-turn ========#


# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (0,  1,   0, 1, 0)
cdt_f = (3, -3.5, 0, 1, 0)

# collisions
l1_c = 0.5 # largeur du mur
l2_c = 3.0 # longueur du mur 
l3_c = 2.5 # espacement entre les virages
C = [0 1; 1 0; 0 -1;   0 1; -1 0; 0 -1]
d = [l1_c/2; l2_c; l1_c/2;  l1_c/2-l3_c; 0; l1_c/2+l3_c]
poly = [[[-10; l1_c/2], [l2_c; l1_c/2], [l2_c; -l1_c/2], [-10; -l1_c/2]],
        [[20; l1_c/2-l3_c], [0; l1_c/2-l3_c], [0; -l1_c/2-l3_c], [20; -l1_c/2-l3_c]]]

function plot_s_turn!(ax)
        lines!(ax, [(-1, l1_c/2), (l2_c, l1_c/2), (l2_c, -l1_c/2), (-1, -l1_c/2)]; color=:blue)
        lines!(ax, [(l2_c+1, l1_c/2-l3_c), (0, l1_c/2-l3_c), (0, -l1_c/2-l3_c), (l2_c+1, -l1_c/2-l3_c)]; color=:blue)
end

for nh in nh_p
    start = initAstar(nh, poly, (-2, 6), (-4, 2), cdt_0, cdt_f; spacing=0.2, dmin=0.6, smoothing=nh÷8)

    for (rk_name, rk) in collocation_p
        for (model_name, m!) in model_p
            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2012_collisions!(model, nh, (2, 3), C, d)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "S-turn", model_name, rk_name, "2012", ax -> plot_s_turn!(ax), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_collisions!(model, nh, (2, 3), C, d)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "S-turn", model_name, rk_name, "2017", ax -> plot_s_turn!(ax), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_penetration_collisions!(model, nh, (2, 3), C, d; kappa=1e+3)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "S-turn", model_name, rk_name, "2017 penetration", ax -> plot_s_turn!(ax), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "S-turn", model_name, rk_name, "2023", ax -> plot_s_turn!(ax), start)
        end
    end
end







#======== parking ========#


# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-2.5,  1.2,  0, 0, 0)
cdt_f = (-0.5, -0.35, 0, 0, 0)

# collisions
lp = 2.3 # longueur de la place de parking
pp = 0.8 # profondeur de la place de parking
rp = 2.6 # largeur de la route
C = [0 1; 1 0;   0 1; 0 -1;   -1 0; 0 1;   0 -1; 0 1]
d = [0; -lp/2;   -pp; pp+10;   -lp/2; 0;   -rp; rp+10]
poly=[[[-10;0], [-lp/2;0], [-lp/2;-10], [-10;-10]],
      [[-lp/2; -pp], [-lp/2;-pp-10], [lp/2; -pp-10], [lp/2;-pp]],
      [[10;0], [lp/2;0], [lp/2;-10], [10;-10]],
      [[-10;rp], [-10;rp+10], [10;rp+10], [10;rp]]]

function plot_parking!(ax)    
    lines!(ax, [(-3, 0), (-lp/2, 0), (-lp/2, -pp), (lp/2, -pp), (lp/2, 0), (3, 0)], color=:blue)
    lines!(ax, [(-3, rp), (3, rp)], color=:blue)
end


for nh in nh_p
    #start = init_straight(nh, [cdt_0[1:2]...], [cdt_f[1:2]...], 2.7)
    start = initAstar(nh, poly, (-4, 2), (-1.5, 3), cdt_0, cdt_f; spacing=0.2, dmin=0.3, smoothing=nh÷8)

    for (rk_name, rk) in collocation_p
        for (model_name, m!) in model_p
            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2012_collisions!(model, nh, (4, 2), C, d)
            limit_time!(model, 4)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "parking", model_name, rk_name, "2012", ax -> plot_parking!(ax), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_collisions!(model, nh, (4, 2), C, d)
            limit_time!(model, 4)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "parking", model_name, rk_name, "2017", ax -> plot_parking!(ax), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2017_penetration_collisions!(model, nh, (4, 2), C, d; kappa=1e+3)
            limit_time!(model, 4)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "parking", model_name, rk_name, "2017 penetration", ax -> plot_parking!(ax), start)

            model = Model()
            m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
            Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)
            limit_time!(model, 4)
            solve!(model, max_iter=1000, verbose=0)
            save_result(model, nh, "parking", model_name, rk_name, "2023", ax -> plot_parking!(ax), start)
        end
    end
end








#======== poly-field ========#


# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-8, 0, 0, 0, 0)
cdt_f = ( 8, 0, 0, 0, 0)

# collisions
N_faces = 5
poly_area = 1;
seeds = [3, 5, 7, 11, 14]

function plot_poly_field!(ax, poly)
    for points in poly
        poly!(ax, [Tuple(p) for p in points], strokecolor=:blue, strokewidth=1, color=:white)
    end
end

for seed in seeds
    Random.seed!(seed)

    # matrices et vecteurs représentatifs des polygones
    poly_C = []
    poly_d = []
    poly = []
    center = PoissonDiskSampling.generate(2.8, (-6, 6), (-4, 4))
    N_poly = length(center)

    for c in center
        points = convex_volume_polygon(N_faces, poly_area)
        points = [p .+ c for p in points]

        # collision polyhedra
        edges = [pointsToEdge(points[i], points[(i)%N_faces+1], points[(i+1)%N_faces+1]) for i in eachindex(points)]
        append!(poly_C, (n -> n[1]).(edges))
        append!(poly_d, (n -> n[2]).(edges))

        # collision polytop
        push!(poly, points)
    end
    C = stack(poly_C, dims=1)
    d = reshape(poly_d, (N_faces*N_poly, 1))

    for nh in nh_p
        start = initAstar(nh, poly, (-8.2, 8.2), (-5, 5), cdt_0, cdt_f; spacing=0.3, dmin=0.4)

        for (rk_name, rk) in collocation_p
            for (model_name, m!) in model_p
                model = Model()
                m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
                Npoly_rect_2012_collisions!(model, nh, (N_poly, N_faces), C, d)
                solve!(model, max_iter=1000, verbose=0)
                save_result(model, nh, "poly-field-$(seed)", model_name, rk_name, "2012", ax -> plot_poly_field!(ax, poly), start)

                model = Model()
                m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
                Npoly_rect_2017_collisions!(model, nh, (N_poly, N_faces), C, d)
                solve!(model, max_iter=1000, verbose=0)
                save_result(model, nh, "poly-field-$(seed)", model_name, rk_name, "2017", ax -> plot_poly_field!(ax, poly), start)

                model = Model()
                m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
                Npoly_rect_2017_penetration_collisions!(model, nh, (N_poly, N_faces), C, d; kappa=1e+3)
                solve!(model, max_iter=1000, verbose=0)
                save_result(model, nh, "poly-field-$(seed)", model_name, rk_name, "2017 penetration", ax -> plot_poly_field!(ax, poly), start)

                model = Model()
                m!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
                Npoly_rect_2023_collisions!(model, nh, poly; d_min=0.05)
                solve!(model, max_iter=1000, verbose=0)
                save_result(model, nh, "poly-field-$(seed)", model_name, rk_name, "2023", ax -> plot_poly_field!(ax, poly), start)
            end
        end
    end
end







# save results in .csv format
CSV.write("comp_figures/results.csv", results)
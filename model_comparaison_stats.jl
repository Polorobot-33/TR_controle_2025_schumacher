using JuMP
using Random
using LinearAlgebra
using PoissonDiskSampling

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")
include("polygon.jl")

# nombre de tests
N = 2

mod2012  = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)
mod2017a = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)
mod2017b = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)
mod2023  = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)


nh = 63
# conditions initales et finales (x, y, ϕ, u, r)
cdt_0 = (-7, 0, 0, 0, 0)
cdt_f = ( 7, 0, 0, 0, 0)

N_faces = 5
poly_area = 1;


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

for i in 1:N
    # matrices et vecteurs représentatifs des polygones
    # collision using polyhedra
    poly_C = []
    poly_d = []

    # collision using polytop
    poly = []

    center = PoissonDiskSampling.generate(2.5, (-4, 4), (-3, 3))
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

    start = initAstar(nh, poly, (-7.5, 7.5), (-5, 5), cdt_0, cdt_f; spacing=0.3, dmin=0.4)

    # collision polyhedra
    model1 = Model()
    dynamic_model!(model1, nh, cdt_0, cdt_f; start=start, m=105)
    Npoly_rect_2012_collisions!(model1, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
    solve!(model1, max_iter=1000)


    model2a = Model()
    dynamic_model!(model2a, nh, cdt_0, cdt_f; start=start, m=105)
    Npoly_rect_2017_collisions!(model2a, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
    solve!(model2a, max_iter=1000)

    model2b = Model()
    dynamic_model!(model2b, nh, cdt_0, cdt_f; start=start, m=105)
    Npoly_rect_2017_penetration_collisions!(model2b, nh, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)); kappa=1e+2)
    solve!(model2b, max_iter=1000)

    # collision polytop
    model3 = Model()
    dynamic_model!(model3, nh, cdt_0, cdt_f; start=start, m=105)
    Npoly_rect_2023_collisions!(model3, nh, poly; d_min=0.05)
    solve!(model3, max_iter=1000)

    # extract stats
    mod2012["tps_trajet"] += JuMP.value.(model1[:T])
    mod2012["tps_calcul"] += JuMP.solve_time(model1)
    mod2012["nb_iter"]    += JuMP.barrier_iterations(model1)

    mod2017a["tps_trajet"] += JuMP.value.(model2a[:T])
    mod2017a["tps_calcul"] += JuMP.solve_time(model2a)
    mod2017a["nb_iter"]    += JuMP.barrier_iterations(model2a)

    mod2017b["tps_trajet"] += JuMP.value.(model2b[:T])
    mod2017b["tps_calcul"] += JuMP.solve_time(model2b)
    mod2017b["nb_iter"]    += JuMP.barrier_iterations(model2b)

    mod2023["tps_trajet"] += JuMP.value.(model3[:T])
    mod2023["tps_calcul"] += JuMP.solve_time(model3)
    mod2023["nb_iter"]    += JuMP.barrier_iterations(model3)
end

println("\t\t2012\t2017a\t2017b\t2023");
for (m1, m2a, m2b, m3) in zip(mod2012, mod2017a, mod2017b, mod2023)
    field = m1[1]
    
    print(field)
    print(" :\t")

    v1 = m1[2]
    v2 = m2a[2]
    v3 = m2b[2]
    v4 = m3[2]

    printstyled(round(v1/N; digits=3), color=(v1<v2 && v1<v3 && v1<v4) ? :green : :red)
    print("\t")
    printstyled(round(v2/N; digits=3), color=(v2<v1 && v2<v3 && v2<v4) ? :green : :red)
    print("\t")
    printstyled(round(v3/N; digits=3), color=(v3<v1 && v3<v2 && v3<v4) ? :green : :red)
    print("\t")
    printstyled(round(v4/N; digits=3), color=(v4<v1 && v4<v2 && v4<v3) ? :green : :red)
    print("\n");
end
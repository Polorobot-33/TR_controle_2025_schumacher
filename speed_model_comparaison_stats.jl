using JuMP
using Random
using LinearAlgebra
using PoissonDiskSampling

include("models.jl")
include("renderer.jl")
include("initial_conditions.jl")
include("polygon.jl")

# nombre de tests
N = 20

mod2012 = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)
mod2017 = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)
mod2023 = Dict("tps_trajet" => 0., "tps_calcul" => 0., "nb_iter" => 0.)


nh = 99
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

    # collision polyhedra
    model1, init = robot_rect_custom_model(nh, cdt_0, cdt_f, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
    JuMP.set_optimizer(model1, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model1, "max_iter", 500)
    JuMP.optimize!(model1)

    model2, _ = robot_rect_custom_polyhedra(nh, cdt_0, cdt_f, (N_poly, N_faces), stack(poly_C, dims=1), reshape(poly_d, (N_faces*N_poly, 1)))
    JuMP.set_optimizer(model2, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model2, "max_iter", 500)
    JuMP.optimize!(model2)

    # collision polytop
    model3, init = robot_rect_custom_polytop(nh, cdt_0, cdt_f, poly; d_min=1e-1)
    JuMP.set_optimizer(model3, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model3, "max_iter", 500)
    JuMP.optimize!(model3)

    # extract stats
    mod2012["tps_trajet"] += JuMP.value.(model1[:T])
    mod2012["tps_calcul"] += JuMP.solve_time(model1)
    mod2012["nb_iter"]    += JuMP.barrier_iterations(model1)

    mod2017["tps_trajet"] += JuMP.value.(model2[:T])
    mod2017["tps_calcul"] += JuMP.solve_time(model2)
    mod2017["nb_iter"]    += JuMP.barrier_iterations(model2)

    mod2023["tps_trajet"] += JuMP.value.(model3[:T])
    mod2023["tps_calcul"] += JuMP.solve_time(model3)
    mod2023["nb_iter"]    += JuMP.barrier_iterations(model3)
end

println("\t\t2012\t2017\t2023");
for (m1, m2, m3) in zip(mod2012, mod2017, mod2023)
    field = m1[1]
    
    print(field)
    print(" :\t")

    v1 = m1[2]
    v2 = m2[2]
    v3 = m3[2]

    printstyled(round(v1/N; digits=3), color=(v1<v2 && v1<v3) ? :green : :red)
    print("\t")
    printstyled(round(v2/N; digits=3), color=(v2<v1 && v2<v3) ? :green : :red)
    print("\t")
    printstyled(round(v3/N; digits=3), color=(v3<v1 && v3<v2) ? :green : :red)
    print("\n");
end
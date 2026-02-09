using JuMP
using Ipopt
using LinearAlgebra

include("collocation.jl")
include("initial_conditions.jl")

function var(model, var)
    return JuMP.value.(model[var]).data
end

function limit_time!(model, T_max)
    T = model[:T]
    @constraint(model, T <= T_max)
end



function robot_point_model(nh)
    l1_c = 2.4 # largeur principale des couloirs
    l2_c = 0.9 # largeur des couloirs latéraux
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot
    d_c  = rw_c * 0.4 # espacement des roues * 0.5
    u_max = 2.7 # vitesse maximale d'une roue
    r_max = 1.5 # vitesse maximale de rotation

    # conditions initiales
    x_0 = 0
    y_0 = -3
    ϕ_0 = π/2
    u_0 = u_max
    r_0 = 0

    # conditions finales
    x_f = 3
    y_f = 0
    ϕ_f = 0
    u_f = u_max
    r_f = 0

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_arc(nh, -3, 3, u_max/2, 1)
    #x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [0, -3], [3, 0], u_max/2)
    step = 1 / nh

    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [0; 0; 0; 0]
    C = [0 -1;    -1 0;     1 0;     0 -1;    0 1;     1 0;    -1 0;     0 1]
    d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
    S(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    R(x, y) = [x ; y]
    epsilon = 1e-6

    model = Model()

    @variables(model, begin
        x[i=0:nh],                    (start = x_i[i+1])
        y[i=0:nh],                    (start = y_i[i+1])
        ϕ[i=0:nh],                    (start = ϕ_i[i+1])
        0.0    <= u[i=0:nh] <= u_max, (start = u_i[i+1])
        -r_max <= r[i=0:nh] <= r_max, (start = r_i[i+1])
        0.0 <= T,                     (start = T_i)
        ω[t=0:nh, i=0:3, j=0:5] >= 0, (start = 0)
    end)

    @expressions(model, begin
        δt, T * step
    end)

    @objective(model, Min, T)

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * 0.5 * cos((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_y[i=1:nh], y[i] == y[i-1] + δt * 0.5 * sin((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * 0.5 * (r[i-1] + r[i])
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:3], transpose([A * transpose(S(ϕ[t])); C[[2*i+1, 2*i+2], :]]) * ω[t, i, :] == 0
        con_c2[t=1:nh, i=0:3], transpose([b .+ A * transpose(S(ϕ[t])) * R(x[t], y[t]); d[[2*i+1, 2*i+2]]]) * ω[t, i, :] <= -epsilon
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)

    return model, (x_i, y_i)
end

function robot_rect_model(nh; u_max_p=2.7, r_max_p=1.5, epsilon_p=1e-6)
    l1_c = 2.4 # largeur principale des couloirs
    l2_c = 0.9 # largeur des couloirs latéraux
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot
    d_c  = rw_c * 0.4 # espacement des roues * 0.5
    u_max = u_max_p # vitesse maximale d'une roue
    r_max = r_max_p # vitesse maximale de rotation

    # conditions initiales
    x_0 = 0
    y_0 = -3
    ϕ_0 = π/2
    u_0 = u_max
    r_0 = 0

    # conditions finales
    x_f = 3
    y_f = 0
    ϕ_f = 0
    u_f = u_max
    r_f = 0

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_arc(nh, -3, 3, u_max/2, r_max*0.8, 0.5)
    #x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [0, -3], [3, 0], u_max/2)
    step = 1 / nh

    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = [0 -1; -1 0; 1 0; 0 -1; 0 1; 1 0; -1 0; 0 1]
    d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
    S(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    R(x, y) = [x ; y]
    epsilon = epsilon_p

    model = Model()

    @variables(model, begin
        x[i=0:nh],                    (start = x_i[i+1])
        y[i=0:nh],                    (start = y_i[i+1])
        ϕ[i=0:nh],                    (start = ϕ_i[i+1])
        0.0    <= u[i=0:nh] <= u_max, (start = u_i[i+1])
        -r_max <= r[i=0:nh] <= r_max, (start = r_i[i+1])
        0.0 <= T,                     (start = T_i)
        ω[t=0:nh, i=0:3, j=0:5] >= 0, (start = 0)
    end)

    @expressions(model, begin
        δt, T * step
    end)

    @objective(model, Min, T)

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * 0.5 * cos((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_y[i=1:nh], y[i] == y[i-1] + δt * 0.5 * sin((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * 0.5 * (r[i-1] + r[i])
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:3], transpose([A * transpose(S(ϕ[t])); C[[2*i+1, 2*i+2], :]]) * ω[t, i, :] == 0
        con_c2[t=1:nh, i=0:3], transpose([b .+ A * transpose(S(ϕ[t])) * R(x[t], y[t]); d[[2*i+1, 2*i+2]]]) * ω[t, i, :] <= -epsilon
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)

    return model, (x_i, y_i)
end

# no collision, rectangular hitbox
function robot_rect(nh, cdt_0, cdt_f; u_max_p=2.7, r_max_p=1.5)
    rw_c = 0.72 # largeur du robot
    d_c  = rw_c * 0.4 # espacement des roues * 0.5
    u_max = u_max_p # vitesse maximale d'une roue
    r_max = r_max_p # vitesse maximale de rotation

    # conditions initiales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0

    # conditions finales
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [x_0, y_0], [x_f, y_f], u_max/2)
    step = 1 / nh

    model = Model()

    @variables(model, begin
        x[i=0:nh],                    (start = x_i[i+1])
        y[i=0:nh],                    (start = y_i[i+1])
        ϕ[i=0:nh],                    (start = ϕ_i[i+1])
        -u_max <= u[i=0:nh] <= u_max, (start = u_i[i+1])
        -r_max <= r[i=0:nh] <= r_max, (start = r_i[i+1])
        0.0 <= T,                     (start = T_i)
    end)

    @expressions(model, begin
        δt, T * step
    end)

    @objective(model, Min, T)

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * 0.5 * cos((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_y[i=1:nh], y[i] == y[i-1] + δt * 0.5 * sin((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * 0.5 * (r[i-1] + r[i])
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)

    return model, (x_i, y_i)
end

function robot_rect_custom_model(nh, cdt_0, cdt_f, dims, C_p, d_p; u_max_p=2.7, r_max_p=1.5, epsilon_p=1e-6)
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot
    d_c  = rw_c * 0.4 # espacement des roues * 0.5
    u_max = u_max_p # vitesse maximale d'une roue
    r_max = r_max_p # vitesse maximale de rotation

    # conditions initiales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0

    # conditions finales
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [x_0, y_0], [x_f, y_f], u_max/2)
    step = 1 / nh

    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = C_p
    d = d_p
    N_poly, N_faces = dims
    S(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    R(x, y) = [x ; y]
    epsilon = epsilon_p

    model = Model()

    @variables(model, begin
        x[i=0:nh],                    (start = x_i[i+1])
        y[i=0:nh],                    (start = y_i[i+1])
        ϕ[i=0:nh],                    (start = ϕ_i[i+1])
        -u_max <= u[i=0:nh] <= u_max, (start = u_i[i+1])
        -r_max <= r[i=0:nh] <= r_max, (start = r_i[i+1])
        0.0 <= T,                     (start = T_i)
        ω[t=0:nh, i=0:(N_poly-1), j=0:(4+N_faces-1)] >= 0, (start = 0)
    end)

    @expressions(model, begin
        δt, T * step
    end)

    @objective(model, Min, T)

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * 0.5 * cos((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_y[i=1:nh], y[i] == y[i-1] + δt * 0.5 * sin((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * 0.5 * (r[i-1] + r[i])
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose([A * transpose(S(ϕ[t])); C[N_faces*i+1 : N_faces*(i+1), :]]) * ω[t, i, :] == 0
        con_c2[t=1:nh, i=0:(N_poly-1)], transpose([b .+ A * transpose(S(ϕ[t])) * R(x[t], y[t]); d[N_faces*i+1 : N_faces*(i+1)]]) * ω[t, i, :] <= -epsilon
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)

    return model, (x_i, y_i)
end

function robot_rect_custom_polyhedra(nh, cdt_0, cdt_f, dims, C_p, d_p; u_max_p=2.7, r_max_p=1.5, kappa_p=10)
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot
    d_c  = rw_c * 0.4 # espacement des roues * 0.5
    u_max = u_max_p # vitesse maximale d'une roue
    r_max = r_max_p # vitesse maximale de rotation

    # conditions initiales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0

    # conditions finales
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [x_0, y_0], [x_f, y_f], u_max/2)
    step = 1 / nh

    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = C_p
    d = d_p
    N_poly, N_faces = dims
    R(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    tr(x, y) = [x ; y]

    model = Model()

    @variables(model, begin
        x[i=0:nh],                    (start = x_i[i+1])
        y[i=0:nh],                    (start = y_i[i+1])
        ϕ[i=0:nh],                    (start = ϕ_i[i+1])
        -u_max <= u[i=0:nh] <= u_max, (start = u_i[i+1])
        -r_max <= r[i=0:nh] <= r_max, (start = r_i[i+1])
        0.0 <= T,                     (start = T_i)
        s[t=0:nh, i=0:(N_poly-1)]                  >= 0, (start = 0)
        μ[t=0:nh, i=0:(N_poly-1), j=0:3]           >= 0, (start = 0)
        λ[t=0:nh, i=0:(N_poly-1), j=0:(N_faces-1)] >= 0, (start = 0)
    end)

    @expressions(model, begin
        δt, T * step
    end)

    @objective(model, Min, T + kappa_p * sum(s))

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * 0.5 * cos((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_y[i=1:nh], y[i] == y[i-1] + δt * 0.5 * sin((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * 0.5 * (r[i-1] + r[i])
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose(transpose(C[N_faces*i+1 : N_faces*(i+1), :]) * λ[t, i, :]) * (transpose(C[N_faces*i+1 : N_faces*(i+1), :]) * λ[t, i, :]) == 1
        con_c2[t=1:nh, i=0:(N_poly-1)], transpose(A) * μ[t, i, :] + transpose(C[N_faces*i+1 : N_faces*(i+1), :] * R(ϕ[t])) * λ[t, i, :] == 0
        con_c3[t=1:nh, i=0:(N_poly-1)], -transpose(b) * μ[t, i, :] + transpose(C[N_faces*i+1 : N_faces*(i+1), :] * tr(x[t], y[t]) - d[N_faces*i+1 : N_faces*(i+1)]) * λ[t, i, :] >= -s[t, i]
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)

    return model, (x_i, y_i)
end

function robot_rect_custom_polytop(nh, cdt_0, cdt_f, poly; u_max_p=2.7, r_max_p=1.5, d_min=1e-4)
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot
    d_c  = rw_c * 0.4 # espacement des roues ÷ 2
    u_max = u_max_p # vitesse maximale d'une roue
    r_max = r_max_p # vitesse maximale de rotation

    # conditions initiales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0

    # conditions finales
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [x_0, y_0], [x_f, y_f], u_max/2)
    step = 1 / nh

    # contraintes de Collisions
    N_poly = length(poly)
    Ve = [(-rh_c/2) (rh_c/2) (rh_c/2) (-rh_c/2); (rw_c/2) (rw_c/2) (-rw_c/2) (-rw_c/2)]
    Vo = [reduce(hcat, vertices) for vertices in poly]

    R(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    tr(x, y) = [x ; y]

    model = Model()

    @variables(model, begin
        x[i=0:nh],                    (start = x_i[i+1])
        y[i=0:nh],                    (start = y_i[i+1])
        ϕ[i=0:nh],                    (start = ϕ_i[i+1])
        -u_max <= u[i=0:nh] <= u_max, (start = u_i[i+1])
        -r_max <= r[i=0:nh] <= r_max, (start = r_i[i+1])
        0.0 <= T,                     (start = T_i)
        # collisions
        ξ[t=0:nh, i=0:(N_poly-1), j=0:1], (start = 0.0)
        μ_r[t=0:nh, i=0:(N_poly-1)], (start = 0.0)
        μ_o[t=0:nh, i=0:(N_poly-1)], (start = 0.0)
    end)

    @expressions(model, begin
        δt, T * step
    end)

    @objective(model, Min, T)

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * 0.5 * cos((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_y[i=1:nh], y[i] == y[i-1] + δt * 0.5 * sin((ϕ[i] + ϕ[i-1]) * 0.5) * (u[i] + u[i-1])
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * 0.5 * (r[i-1] + r[i])
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose(ξ[t, i, :]) * ξ[t, i, :] / 4 + μ_r[t, i] + μ_o[t, i] + d_min^2 <= 0
        con_c2[t=1:nh, i=0:(N_poly-1)], -transpose(R(ϕ[t]) * Ve .+ tr(x[t], y[t])) * ξ[t, i, :] .- μ_r[t, i] <= 0
        con_c3[t=1:nh, i=0:(N_poly-1)], transpose(Vo[i+1]) * ξ[t, i, :] .- μ_o[t, i] <= 0
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)

    return model, (x_i, y_i)
end


# model = Model();
# speed_model!(model, ...)
# collisions!(model, ...)
# solve()

function base_model!(model, nh, cdt_0, cdt_f; start=nothing, u_max_p=2.75, rk::RKScheme=CrankNicolson())
    nc = length(rk.c)
    
    u_max = u_max_p # vitesse maximale d'une roue
    x_0, y_0, ϕ_0, _= cdt_0 # conditions initiales
    x_f, y_f, ϕ_f, _ = cdt_f # conditions finales

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = (isnothing(start) ? init_straight(nh, [x_0, y_0], [x_f, y_f], u_max/2) : start)

    step = 1 / nh

    @variables(model, begin
        T >= 0.1,  (start = T_i)

        qx[i=0:nh, j=0:nc], (start = x_i[i+1])
        qy[i=0:nh, j=0:nc], (start = y_i[i+1])
        qϕ[i=0:nh, j=0:nc], (start = ϕ_i[i+1])
    end)

    @expressions(model, begin 
        δt, T * step
        
        x[i=0:nh], qx[i, 0]
        y[i=0:nh], qy[i, 0]
        ϕ[i=0:nh], qϕ[i, 0]
    end)

    @objective(model, Min, T)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
    end)

    return x_i, y_i, ϕ_i, u_i, r_i, T_i
end

function speed_model!(model, nh, cdt_0, cdt_f; start=nothing, u_max=2.7, r_max=5.0, rw_c=0.72, rk::RKScheme=CrankNicolson())
    _, _, _, u_i, r_i, _ = base_model!(model, nh, cdt_0, cdt_f; start=start, u_max_p=u_max, rk=rk)
    nc = length(rk.c)

    d_c = rw_c*0.4;

    @variables(model, begin
        u[i=0:nh], (start = u_i[i+1])
        r[i=0:nh], (start = r_i[i+1])
    end)
    
    qx = model[:qx]
    qy = model[:qy]
    qϕ = model[:qϕ]
    δt = model[:δt]

    @expressions(model, begin
        dx[i=0:nh-1, j=0:nc-1], u[i] * cos(qϕ[i, j])
        dy[i=0:nh-1, j=0:nc-1], u[i] * sin(qϕ[i, j])
        dϕ[i=0:nh-1, j=0:nc-1], r[i]
    end)

    # Dynamics
    @constraints(model, begin
        con_dx[i=0:nh-1, j=1:nc], qx[i, j] == qx[i, 0] + δt * sum(rk.a[j, k+1] * dx[i, k] for k in 0:nc-1)
        con_dy[i=0:nh-1, j=1:nc], qy[i, j] == qy[i, 0] + δt * sum(rk.a[j, k+1] * dy[i, k] for k in 0:nc-1)
        con_dϕ[i=0:nh-1, j=1:nc], qϕ[i, j] == qϕ[i, 0] + δt * sum(rk.a[j, k+1] * dϕ[i, k] for k in 0:nc-1)

        [i=0:nh-1], qx[i+1, 0] == qx[i, 0] + δt * sum(rk.b[j+1] * dx[i, j] for j in 0:nc-1)
        [i=0:nh-1], qy[i+1, 0] == qy[i, 0] + δt * sum(rk.b[j+1] * dy[i, j] for j in 0:nc-1)
        [i=0:nh-1], qϕ[i+1, 0] == qϕ[i, 0] + δt * sum(rk.b[j+1] * dϕ[i, j] for j in 0:nc-1)
    end)


    _, _, _, u_0, r_0 = cdt_0
    _, _, _, u_f, r_f = cdt_f

    # Commands
    @constraints(model, begin
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique

        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)
end

function dynamic_model!(model, nh, cdt_0, cdt_f; start=nothing, U_max=48, rw=0.72, rh=1.128, m=105, rk=CrankNicolson())
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.05;
    d = rw*0.4; # espacement des roues

    R = 0.25;# 1.43; # résistance du moteur
    K = 4#6.53 # rapport du réducteur
    k = 0.107 * K #0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues

    _, _, _, u_i, r_i, _ = base_model!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
    nc = length(rk.c)

    @variables(model, begin
        qu[i=0:nh, j=0:nc], (start = u_i[i+1])
        qr[i=0:nh, j=0:nc], (start = r_i[i+1])
    end)
    
    qx = model[:qx]
    qy = model[:qy]
    qϕ = model[:qϕ]
    δt = model[:δt]

    _, _, _, u_0, r_0 = cdt_0
    _, _, _, u_f, r_f = cdt_f


    @variables(model, begin
        -U_max <= ul[i=0:nh] <= U_max, (start = k * (u_i[i+1] - d*r_i[i+1]) / ρ)
        -U_max <= ur[i=0:nh] <= U_max, (start = k * (u_i[i+1] + d*r_i[i+1]) / ρ)
    end)

    @expressions(model, begin
        τl[i=0:nh], k * (ul[i] - k/ρ*(qu[i, 0] - d*qr[i, 0])) / R
        τr[i=0:nh], k * (ur[i] - k/ρ*(qu[i, 0] + d*qr[i, 0])) / R

        dx[i=0:nh-1, j=0:nc-1], qu[i, j] * cos(qϕ[i, j])
        dy[i=0:nh-1, j=0:nc-1], qu[i, j] * sin(qϕ[i, j])
        dϕ[i=0:nh-1, j=0:nc-1], qr[i, j]
        du[i=0:nh-1, j=0:nc-1], (τl[i] + τr[i]) / (ρ * Meq)
        dr[i=0:nh-1, j=0:nc-1], (τr[i] - τl[i])*d / (ρ * Ieq)
    end)

    # Dynamics
    @constraints(model, begin
        con_dx[i=0:nh-1, j=1:nc], qx[i, j] == qx[i, 0] + δt * sum(rk.a[j, k+1] * dx[i, k] for k in 0:nc-1)
        con_dy[i=0:nh-1, j=1:nc], qy[i, j] == qy[i, 0] + δt * sum(rk.a[j, k+1] * dy[i, k] for k in 0:nc-1)
        con_dϕ[i=0:nh-1, j=1:nc], qϕ[i, j] == qϕ[i, 0] + δt * sum(rk.a[j, k+1] * dϕ[i, k] for k in 0:nc-1)
        con_du[i=0:nh-1, j=1:nc], qu[i, j] == qu[i, 0] + δt * sum(rk.a[j, k+1] * du[i, k] for k in 0:nc-1)
        con_dr[i=0:nh-1, j=1:nc], qr[i, j] == qr[i, 0] + δt * sum(rk.a[j, k+1] * dr[i, k] for k in 0:nc-1)

        [i=0:nh-1], qx[i+1, 0] == qx[i, 0] + δt * sum(rk.b[j+1] * dx[i, j] for j in 0:nc-1)
        [i=0:nh-1], qy[i+1, 0] == qy[i, 0] + δt * sum(rk.b[j+1] * dy[i, j] for j in 0:nc-1)
        [i=0:nh-1], qϕ[i+1, 0] == qϕ[i, 0] + δt * sum(rk.b[j+1] * dϕ[i, j] for j in 0:nc-1)
        [i=0:nh-1], qu[i+1, 0] == qu[i, 0] + δt * sum(rk.b[j+1] * du[i, j] for j in 0:nc-1)
        [i=0:nh-1], qr[i+1, 0] == qr[i, 0] + δt * sum(rk.b[j+1] * dr[i, j] for j in 0:nc-1)
    end)


    @expressions(model, begin
        u[i=0:nh], qu[i, 0]
        r[i=0:nh], qr[i, 0]
    end)

    # Commands
    @constraints(model, begin
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0

        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
    end)
end

function dynamic_model_Tlim!(model, nh, cdt_0, cdt_f; start=nothing, U_max=48, rw=0.72, rh=1.128, m=105, epsilon = 1e-6, f=0.5, rk=CrankNicolson())
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.10;
    d = rw*0.4; # espacement des roues
    l = 0.1; # hauteur du centre de gravit au-dessus de l'axe des roues
    b = rh*0.4 # distance des roulettes


    R = 0.25;#1.43; # résistance du moteur
    k = 0.05#0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues
    g = 9.81

    _, _, _, u0, r0, _ = base_model!(model, nh, cdt_0, cdt_f; start=start, u_max_p=2.7, rk=rk)


    u = model[:u]
    r = model[:r]
    δt = model[:δt]

    @variables(model, begin
        -U_max <= ul[i=0:nh] <= U_max, (start = 0)#k * (u0[i+1] - d*r0[i+1]) / ρ)
        -U_max <= ur[i=0:nh] <= U_max, (start = 0)#k * (u0[i+1] + d*r0[i+1]) / ρ)

        # force tangente au sol
        tl[i=0:nh], (start = 0)
        tr[i=0:nh], (start = 0)
    end)

    @expressions(model, begin
        τl[i=0:nh], k * (ul[i] - k/ρ*(u[i] - d*r[i])) / R
        τr[i=0:nh], k * (ur[i] - k/ρ*(u[i] + d*r[i])) / R

        # forces normale sur chaque pneu
        nl[i=0:nh], (-(l+ρ)/d * Meq*r[i]*u[i] + Meq*g) / 2
        nr[i=0:nh], ( (l+ρ)/d * Meq*r[i]*u[i] + Meq*g) / 2
    end)    

    # Dynamics
    @constraints(model, begin
        con_u[i=1:nh], u[i] == u[i-1] + δt * 0.5 * (τl[i]+τl[i-1] + τr[i]+τr[i-1]) / (ρ * Meq)
        con_r[i=1:nh], r[i] == r[i-1] + δt * 0.5 * (τr[i]+τr[i-1] - τl[i]-τl[i-1])*d / (ρ * Ieq)
    end)

    
    # friction limit
    @constraints(model, begin # assume that nl and nr are always positive (the robot doesn't tip over)
        con_nl[i=0:nh], 0 <= nl[i]
        con_nr[i=0:nh], 0 <= nr[i]

        con_a[i=0:nh], -1.4 <= (τl[i] + τr[i]) / (ρ * Meq) <= 1.4

        con_T[i=0:nh], tl[i] + tr[i] == Meq * u[i] * r[i]
        # con_fl_min[i=0:nh], -f * nl[i] <= tl[i]
        # con_fl_max[i=0:nh],  f * nl[i] >= tl[i]
        # con_fr_min[i=0:nh], -f * nr[i] <= tr[i]
        # con_fr_max[i=0:nh],  f * nr[i] >= tr[i]

        con_fl[i=0:nh], tl[i]^2 + (τl[i]/d)^2 <= (f * nl[i])^2
        con_fr[i=0:nh], tr[i]^2 + (τr[i]/d)^2 <= (f * nr[i])^2
    end)
end

function dynamic_model_Tlim_full!(model, nh, cdt_0, cdt_f; start=nothing, U_max=48, rw=0.72, rh=1.128, m=105, epsilon = 1e-6)
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.10;
    d = rw*0.4; # espacement des roues
    l = 0.1; # hauteur du centre de gravit au-dessus de l'axe des roues
    b = rh*0.4 # distance des roulettes
    f = 10 # coefficient de frottements secs


    R = 0.25;#1.43; # résistance du moteur
    k = 0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues
    g = 9.81

    _, _, _, u0, r0, _ = base_model!(model, nh, cdt_0, cdt_f; start=start, u_max_p=2.7)


    u = model[:u]
    r = model[:r]
    δt = model[:δt]

    @variables(model, begin
        -U_max <= ul[i=0:nh] <= U_max, (start = 0)#k * (u0[i+1] - d*r0[i+1]) / ρ)
        -U_max <= ur[i=0:nh] <= U_max, (start = 0)#k * (u0[i+1] + d*r0[i+1]) / ρ)
        
        # forces (pour les limites de frottements)
        τplus[i=0:nh] >= 0, (start = 0)
        τmoins[i=0:nh] >= 0, (start = 0)
        Tl[i=0:nh], (start = 0) # forces normales à la direction des roues
        Tr[i=0:nh], (start = 0)
    end)

    @expressions(model, begin
        τl[i=0:nh], k * (ul[i] - k/ρ*(u[i] - d*r[i])) / R
        τr[i=0:nh], k * (ur[i] - k/ρ*(u[i] + d*r[i])) / R

        # forces
        τ[i=0:nh], τl[i] + τr[i] 
        nl[i=0:nh], (-(l+ρ) * Meq*r[i]*u[i] + d*Meq*g) / (2*d)
        nr[i=0:nh], ( (l+ρ) * Meq*r[i]*u[i] + d*Meq*g) / (2*d)
    end)

    # Dynamics
    @constraints(model, begin
        con_u[i=1:nh], u[i] == u[i-1] + δt * 0.5 * (τl[i]+τl[i-1] + τr[i]+τr[i-1]) / (ρ * Meq)
        con_r[i=1:nh], r[i] == r[i-1] + δt * 0.5 * (τr[i]+τr[i-1] - τl[i]-τl[i-1])*d / (ρ * Ieq)
    end)

    
    # friction limit
    @constraints(model, begin
        con_comp_τA[i=0:nh], τplus[i] - τmoins[i] == τl[i] + τr[i]
        con_comp_τB[i=0:nh], τplus[i] * τmoins[i] <= epsilon
        con_comp_T[i=0:nh],  Tl[i] + Tr[i] == r[i]*u[i]*Meq
        
        con_comp_N1[i=0:nh], Tl[i] #=+ (τl[i]/ρ)^2=# <= (f * (nl[i]#= - (2*τplus[i]-τ[i]) / b=#))
        con_comp_N3[i=0:nh], Tr[i] #=+ (τr[i]/ρ)^2=# <= (f * (nr[i]#= - (2*τplus[i]-τ[i]) / b=#))
        
        #=con_comp_N1[i=0:nh], Tl[i]^2 #=+ (τl[i]/ρ)^2=# <= (f * (nl[i]#= - (2*τplus[i]-τ[i]) / b=#))^2
        #con_comp_N2[i=0:nh], Tl[i]^2 #=+ (τl[i]/ρ)^2=# <= (f * (nl[i]#= - (2*τmoins[i]+τ[i]) / b=#))^2
        con_comp_N3[i=0:nh], Tr[i]^2 #=+ (τr[i]/ρ)^2=# <= (f * (nr[i]#= - (2*τplus[i]-τ[i]) / b=#))^2
        #con_comp_N4[i=0:nh], Tr[i]^2 #=+ (τr[i]/ρ)^2=# <= (f * (nr[i]#= - (2*τmoins[i]+τ[i]) / b=#))^2=#
    end)
end

function dynamic_model_slide!(model, nh, cdt_0, cdt_f; start=nothing, U_max=48, rw=0.72, rh=1.128, m=105, μ=0.5, epsilon = 1e-6)
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.10;
    d = rw*0.4; # espacement des roues
    # μ coefficient de friction

    R = 0.25;# 1.43; # résistance du moteur
    k = 0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues

    # conditions initiales et finales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f
    # etat initial
    x_i, y_i, ϕ_i, u_i, _, T_i = (isnothing(start) ? init_straight(nh, [x_0, y_0], [x_f, y_f], (u_0 + u_f) / 2) : start)
    step = 1 / nh

    @variables(model, begin
        x[i=0:nh], (start = x_i[i+1])
        y[i=0:nh], (start = y_i[i+1])
        ϕ[i=0:nh], (start = ϕ_i[i+1])
        ϕdot[i=0:nh], (start = 0)
        xdot[i=0:nh], (start = u_i[i+1]*cos(ϕ_i[i+1]))
        ydot[i=0:nh], (start = u_i[i+1]*sin(ϕ_i[i+1]))
        -1.4 <= a[i=0:nh] <= 1.4 , (start = 0)
        -2.5 <= r[i=0:nh] <= 2.5, (start = 0)
        
        # friction force
        xplus[i=0:nh] >= 0, (start = 0)
        xmoins[i=0:nh] >= 0, (start = 0)
        -1 <= sgn[i=0:nh] <= 1, (start = 0)

        0.0 <= T,  (start = T_i)
    end)

    @expressions(model, begin
        δt, T * step

        e_t[i=0:nh], [cos(ϕ[i]); sin(ϕ[i])]
        e_n[i=0:nh], [-sin(ϕ[i]); cos(ϕ[i])]
        v[i=0:nh], [xdot[i]; ydot[i]]

        F[i=0:nh], -μ * 1 * 9.81 * sgn[i] # Meq / Meq = 1
    end)

    @objective(model, Min, T)# + (1e+10)*sum(abs.((transpose.(v) .* e_t) .* ϕdot .- F)))

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * (xdot[i-1] + xdot[i]) / 2
        con_y[i=1:nh], y[i] == y[i-1] + δt * (ydot[i-1] + ydot[i]) / 2
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * (ϕdot[i-1] + ϕdot[i]) / 2

        con_v[i=1:nh], v[i] == v[i-1] + δt * (#=a[i-1].*e_t[i-1] .+ F[i-1].*e_n[i-1] .+ =#a[i].*e_t[i] .+ F[i].*e_n[i]) / 2
        con_ϕdot[i=1:nh], ϕdot[i] == ϕdot[i-1] + δt * (#=r[i-1] * (transpose(e_t[i-1])*v[i-1]) .+ =#r[i] * (transpose(e_t[i])*v[i])) / 2
    
        # frottements
        con_xsgn[i=0:nh], xplus[i] - xmoins[i] == transpose(e_n[i]) * v[i]
        con_orth[i=0:nh], xplus[i] * xmoins[i] <= epsilon
        con_sgn1[i=0:nh], (1 - sgn[i]) * xplus[i] <= epsilon
        con_sgn2[i=0:nh], (1 + sgn[i]) * xmoins[i] <= epsilon
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        xdot_ic, xdot[0] == u_0 * cos(ϕ_0)
        ydot_ic, ydot[0] == u_0 * sin(ϕ_0)
        ϕdot_ic, ϕdot[0] == 0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        xdot_fc, xdot[nh] == u_f * cos(ϕ_f)
        ydot_fc, ydot[nh] == u_f * sin(ϕ_f)
        ϕdot_fc, ϕdot[nh] == 0
    end)
end

function dynamic_model_slide2!(model, nh, cdt_0, cdt_f; start=nothing, U_max=48, rw=0.72, rh=1.128, m=105, μ=0.5, epsilon = 1e-6)
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.10;
    d = rw*0.4; # espacement des roues
    # μ coefficient de friction

    R = 0.25;# 1.43; # résistance du moteur
    k = 0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues

    # conditions initiales et finales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f
    # etat initial
    x_i, y_i, ϕ_i, u_i, _, T_i = (isnothing(start) ? init_straight(nh, [x_0, y_0], [x_f, y_f], (u_0 + u_f) / 2) : start)
    step = 1 / nh

    @variables(model, begin
        x[i=0:nh], (start = x_i[i+1])
        y[i=0:nh], (start = y_i[i+1])
        ϕ[i=0:nh], (start = ϕ_i[i+1])
        ϕdot[i=0:nh], (start = 0)
        xdot[i=0:nh], (start = u_i[i+1]*cos(ϕ_i[i+1]))
        ydot[i=0:nh], (start = u_i[i+1]*sin(ϕ_i[i+1]))
        -U_max <= ul[i=0:nh] <= U_max, (start = 0)#k * (u0[i+1] - d*r0[i+1]) / ρ)
        -U_max <= ur[i=0:nh] <= U_max, (start = 0)#k * (u0[i+1] + d*r0[i+1]) / ρ)
        
        # friction force
        xplus[i=0:nh] >= 0, (start = 0)
        xmoins[i=0:nh] >= 0, (start = 0)
        -1 <= sgn[i=0:nh] <= 1, (start = 0)

        0.0 <= T,  (start = T_i)
    end)

    @expressions(model, begin
        δt, T * step

        e_t[i=0:nh], [cos(ϕ[i]); sin(ϕ[i])]
        e_n[i=0:nh], [-sin(ϕ[i]); cos(ϕ[i])]
        v[i=0:nh], [xdot[i]; ydot[i]]

        F[i=0:nh], -μ * 1 * 9.81 * sgn[i] # Meq / Meq = 1
        τl[i=0:nh], k * (ul[i] - k/ρ*(transpose(v[i])*e_t[i] - d*ϕdot[i])) / R
        τr[i=0:nh], k * (ur[i] - k/ρ*(transpose(v[i])*e_t[i] + d*ϕdot[i])) / R
        a[i=0:nh], (τr[i] + τl[i]) / (ρ*Meq)
        r[i=0:nh], (τr[i] - τl[i])*d / (ρ*Ieq)
    end)

    @objective(model, Min, T)# + (1e+8)*sum(abs.((transpose.(v) .* e_t) .* ϕdot .- F)))

    # Dynamics
    @constraints(model, begin
        con_x[i=1:nh], x[i] == x[i-1] + δt * (xdot[i-1] + xdot[i]) / 2
        con_y[i=1:nh], y[i] == y[i-1] + δt * (ydot[i-1] + ydot[i]) / 2
        con_ϕ[i=1:nh], ϕ[i] == ϕ[i-1] + δt * (ϕdot[i-1] + ϕdot[i]) / 2

        con_v[i=1:nh], v[i] == v[i-1] + δt * (a[i-1].*e_t[i-1] .+ F[i-1].*e_n[i-1] .+ a[i].*e_t[i] .+ F[i].*e_n[i]) / 2
        con_ϕdot[i=1:nh], ϕdot[i] == ϕdot[i-1] + δt * (r[i-1] .+ r[i]) / 2
    
        # frottements
        con_xsgn[i=0:nh], xplus[i] - xmoins[i] == transpose(e_n[i]) * v[i]
        con_orth[i=0:nh], xplus[i] * xmoins[i] <= epsilon
        con_sgn1[i=0:nh], (1 - sgn[i]) * xplus[i] <= epsilon
        con_sgn2[i=0:nh], (1 + sgn[i]) * xmoins[i] <= epsilon
    end)

    # Boundary constraints
    @constraints(model, begin
        x_ic, x[0] == x_0
        y_ic, y[0] == y_0
        ϕ_ic, ϕ[0] == ϕ_0
        xdot_ic, xdot[0] == u_0 * cos(ϕ_0)
        ydot_ic, ydot[0] == u_0 * sin(ϕ_0)
        ϕdot_ic, ϕdot[0] == 0

        x_fc, x[nh] == x_f
        y_fc, y[nh] == y_f
        ϕ_fc, ϕ[nh] == ϕ_f
        xdot_fc, xdot[nh] == u_f * cos(ϕ_f)
        ydot_fc, ydot[nh] == u_f * sin(ϕ_f)
        ϕdot_fc, ϕdot[nh] == 0
    end)
end

function dynamic_model_slide31!(model, nh, cdt_0, cdt_f; start=nothing, U_max=48, rw=0.72, rh=1.128, m=105, μ=0.5, rk=CrankNicolson())
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.05;
    d = rw*0.4; # espacement des roues

    R = 0.25;# 1.43; # résistance du moteur
    K = 6.53 # rapport du réducteur
    k = 0.107 * K #0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues
    Iω  = 0.1;
    g = 9.81

    kappa = 1e+2
    epsilon = 1e-5

    _, _, ϕ_i, u_i, r_i, _ = base_model!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
    nc = length(rk.c)

    # state variables
    @variables(model, begin
        vx[i=0:nh, j=0:nc], (start = u_i[i+1]*cos(ϕ_i[i+1]))
        vy[i=0:nh, j=0:nc], (start = u_i[i+1]*sin(ϕ_i[i+1]))
        vϕ[i=0:nh, j=0:nc], (start = r_i[i+1])

        ωl[i=0:nh, j=0:nc], (start = (u_i[i+1] - d*r_i[i+1]) / ρ)
        ωr[i=0:nh, j=0:nc], (start = (u_i[i+1] + d*r_i[i+1]) / ρ)

        s_l[i=0:nh],  (start = 0) # RSG
        ru_l[i=0:nh], (start = 0)
        rv_l[i=0:nh], (start = 0)

        s_r[i=0:nh],  (start = 0) # RSG
        ru_r[i=0:nh], (start = 0)
        rv_r[i=0:nh], (start = 0)
    end)
    
    qx = model[:qx]
    qy = model[:qy]
    qϕ = model[:qϕ]
    δt = model[:δt]

    _, _, _, u_0, r_0 = cdt_0
    _, _, _, u_f, r_f = cdt_f


    @variables(model, begin
        -U_max <= ul[i=0:nh] <= U_max, (start = k * (u_i[i+1] - d*r_i[i+1]) / ρ)
        -U_max <= ur[i=0:nh] <= U_max, (start = k * (u_i[i+1] + d*r_i[i+1]) / ρ)
    end)

    @expressions(model, begin
        vu[i=0:nh, j=0:nc-1], cos(qϕ[i, j])*vx[i, j] + sin(qϕ[i, j])*vy[i, j]

        τl[i=0:nh, j=0:nc-1], k * (ul[i] - k/ρ*(vu[i, j] - d*vϕ[i, j])) / R
        τr[i=0:nh, j=0:nc-1], k * (ur[i] - k/ρ*(vu[i, j] + d*vϕ[i, j])) / R

        dqx[i=0:nh-1, j=0:nc-1], vx[i, j]
        dqy[i=0:nh-1, j=0:nc-1], vy[i, j]
        dqϕ[i=0:nh-1, j=0:nc-1], vϕ[i, j]
        
        dvx[i=0:nh-1, j=0:nc-1], ((ru_r[i]+ru_l[i]) * cos(qϕ[i, nc]) + (rv_r[i]+rv_l[i]) * sin(qϕ[i, nc])) / Meq
        dvy[i=0:nh-1, j=0:nc-1], ((ru_r[i]+ru_l[i]) * cos(qϕ[i, nc]) + (rv_r[i]+rv_l[i]) * sin(qϕ[i, nc])) / Meq
        dvϕ[i=0:nh-1, j=0:nc-1], (ru_r[i] - ru_l[i]) / Ieq * d

        dωr[i=0:nh-1, j=0:nc-1], (τr[i, 0] - ρ * ru_r[i]) / Iω
        dωl[i=0:nh-1, j=0:nc-1], (τl[i, 0] - ρ * ru_l[i]) / Iω
    end)

    # Dynamics
    @constraints(model, begin
        con_dqx[i=0:nh-1, j=1:nc], qx[i, j] == qx[i, 0] + δt * sum(rk.a[j, k+1] * dqx[i, k] for k in 0:nc-1)
        con_dqy[i=0:nh-1, j=1:nc], qy[i, j] == qy[i, 0] + δt * sum(rk.a[j, k+1] * dqy[i, k] for k in 0:nc-1)
        con_dqϕ[i=0:nh-1, j=1:nc], qϕ[i, j] == qϕ[i, 0] + δt * sum(rk.a[j, k+1] * dqϕ[i, k] for k in 0:nc-1)
        con_dvx[i=0:nh-1, j=1:nc], vx[i, j] == vx[i, 0] + δt * sum(rk.a[j, k+1] * dvx[i, k] for k in 0:nc-1)
        con_dvy[i=0:nh-1, j=1:nc], vy[i, j] == vy[i, 0] + δt * sum(rk.a[j, k+1] * dvy[i, k] for k in 0:nc-1)
        con_dωr[i=0:nh-1, j=1:nc], ωr[i, j] == ωr[i, 0] + δt * sum(rk.a[j, k+1] * dωr[i, k] for k in 0:nc-1)
        con_dωl[i=0:nh-1, j=1:nc], ωl[i, j] == ωl[i, 0] + δt * sum(rk.a[j, k+1] * dωl[i, k] for k in 0:nc-1)

        [i=0:nh-1], qx[i+1, 0] == qx[i, 0] + δt * sum(rk.b[j+1] * dqx[i, j] for j in 0:nc-1)
        [i=0:nh-1], qy[i+1, 0] == qy[i, 0] + δt * sum(rk.b[j+1] * dqy[i, j] for j in 0:nc-1)
        [i=0:nh-1], qϕ[i+1, 0] == qϕ[i, 0] + δt * sum(rk.b[j+1] * dqϕ[i, j] for j in 0:nc-1)
        [i=0:nh-1], vx[i+1, 0] == vx[i, 0] + δt * sum(rk.b[j+1] * dvx[i, j] for j in 0:nc-1)
        [i=0:nh-1], vy[i+1, 0] == vy[i, 0] + δt * sum(rk.b[j+1] * dvy[i, j] for j in 0:nc-1)
        [i=0:nh-1], ωr[i+1, 0] == ωr[i, 0] + δt * sum(rk.b[j+1] * dωr[i, j] for j in 0:nc-1)
        [i=0:nh-1], ωl[i+1, 0] == ωl[i, 0] + δt * sum(rk.b[j+1] * dωl[i, j] for j in 0:nc-1)
    end)

    # Friction
    h(vx, vy) = (vx^2 + vy^2)^0.5
    function ∇h(g::AbstractVector{T}, vx::T, vy::T) where {T}
        if vx == 0 && vy == 0
            g[1] = 1
            g[2] = 0
            return
        end
        g[1] = vx / (vx^2 + vy^2)^0.5
        g[2] = vy / (vx^2 + vy^2)^0.5
        return
    end
    function ∇²h(H::AbstractMatrix{T}, vx::T, vy::T) where {T}
        H[1, 1] = 0
        H[2, 1] = 0
        H[2, 2] = 0
        vx == 0 && vy == 0 && return
        H[1, 1] = (1 - vx^2 / (vx^2+vy^2)) / (vx^2 + vy^2)^0.5
        H[2, 1] = -vx*vy / (vx^2 + vy^2)^(3/2)
        H[2, 2] = (1 - vy^2 / (vx^2+vy^2)) / (vx^2 + vy^2)^0.5
        return
    end
    @operator(model, reg_norm, 2, h, ∇h, ∇²h)

    @expressions(model, begin
        vIu_l[i=0:nh], vu[i, 0] - d*vϕ[i, 0] - ωl[i, 0]*ρ
        vIv_l[i=0:nh], -sin(qϕ[i, 0])*vx[i, 0] + cos(qϕ[i, 0])*vy[i, 0]

        vIu_r[i=0:nh], vu[i, 0] + d*vϕ[i, 0] - ωr[i, 0]*ρ
        vIv_r[i=0:nh], -sin(qϕ[i, 0])*vx[i, 0] + cos(qϕ[i, 0])*vy[i, 0]
    end)

    @constraints(model, begin
        #[i=0:nh], -epsilon <= vIu_l[i]*ru_l[i] + vIv_l[i]*rv_l[i] + μ * s_l[i] * g*m/2 <= epsilon
        #[i=0:nh], (ru_l[i]^2 + rv_l[i]^2) <= (g*m/2 * μ)^2
        #[i=0:nh], s_l[i] == reg_norm(vIu_l[i], vIv_l[i])
        #[i=0:nh], s_l[i]^2 == vIu_l[i]^2 + vIv_l[i]^2
        #[i=0:nh], s_l[i] >= 0

        #[i=0:nh], -epsilon <= vIu_r[i]*ru_r[i] + vIv_r[i]*rv_r[i] + μ * s_r[i] * g*m/2 <= epsilon
        #[i=0:nh], (ru_r[i]^2 + rv_r[i]^2) <= (g*m/2 * μ)^2
        #[i=0:nh], s_r[i] == reg_norm(vIu_r[i], vIv_r[i])
        #[i=0:nh], s_r[i]^2 == vIu_r[i]^2 + vIv_r[i]^2
        #[i=0:nh], s_r[i] >= 0
    end)


    @expressions(model, begin
        u[i=0:nh],  vx[i, 0]*cos(qϕ[i, 0]) + vy[i, 0]*sin(qϕ[i, 0])
        v[i=0:nh], -vx[i, 0]*sin(qϕ[i, 0]) + vy[i, 0]*cos(qϕ[i, 0])
        r[i=0:nh],  vϕ[i, 0]
        τ_l[i=0:nh], τl[i, 0]
        τ_r[i=0:nh], τr[i, 0]
    end)

    # Commands
    @constraints(model, begin
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0
        v_ic, v[0] == 0

        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
        v_fc, v[nh] == 0
    end)
end

function dynamic_model_slide32!(model, nh, cdt_0, cdt_f; start=nothing, F_max=50, R_max=10, rw=0.72, rh=1.128, m=105, μ=0.5, rk=CrankNicolson())
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.05;
    d = rw*0.4; # espacement des roues

    R = 0.25;# 1.43; # résistance du moteur
    K = 6.53 # rapport du réducteur
    k = 0.107 * K #0.107; # constante de couple (0.56)

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues
    g = 9.81

    kappa = 1e+2
    epsilon = 1e-6

    _, _, ϕ_i, u_i, r_i, _ = base_model!(model, nh, cdt_0, cdt_f; start=start, rk=rk)
    nc = length(rk.c)

    # state variables
    @variables(model, begin
        vx[i=0:nh, j=0:nc], (start = u_i[i+1]*cos(ϕ_i[i+1]))
        vy[i=0:nh, j=0:nc], (start = u_i[i+1]*sin(ϕ_i[i+1]))
        vϕ[i=0:nh, j=0:nc], (start = r_i[i+1])

        s[i=0:nh]
        ru[i=0:nh]
        rv[i=0:nh]
    end)
    
    qx = model[:qx]
    qy = model[:qy]
    qϕ = model[:qϕ]
    δt = model[:δt]

    _, _, _, u_0, r_0 = cdt_0
    _, _, _, u_f, r_f = cdt_f


    @variables(model, begin
        -F_max <= F[i=0:nh] <= F_max, (start = 0)
        -R_max <= R[i=0:nh] <= R_max, (start = 0)
    end)

    @expressions(model, begin
        vu[i=0:nh],  cos(qϕ[i, 0])*vx[i, 0] + sin(qϕ[i, 0])*vy[i, 0]
        vv[i=0:nh], -sin(qϕ[i, 0])*vx[i, 0] + cos(qϕ[i, 0])*vy[i, 0]

        dqx[i=0:nh-1, j=0:nc-1], vx[i, j]
        dqy[i=0:nh-1, j=0:nc-1], vy[i, j]
        dqϕ[i=0:nh-1, j=0:nc-1], vϕ[i, j]
        
        dvx[i=0:nh-1, j=0:nc-1], ( (ru[i] + F[i]) * cos(qϕ[i, nc]) + rv[i] * sin(qϕ[i, nc])) / Meq
        dvy[i=0:nh-1, j=0:nc-1], (-(ru[i] + F[i]) * sin(qϕ[i, nc]) + rv[i] * cos(qϕ[i, nc])) / Meq
        dvϕ[i=0:nh-1, j=0:nc-1], R[i] / Ieq
    end)

    # Dynamics
    @constraints(model, begin
        con_dqx[i=0:nh-1, j=1:nc], qx[i, j] == qx[i, 0] + δt * sum(rk.a[j, k+1] * dqx[i, k] for k in 0:nc-1)
        con_dqy[i=0:nh-1, j=1:nc], qy[i, j] == qy[i, 0] + δt * sum(rk.a[j, k+1] * dqy[i, k] for k in 0:nc-1)
        con_dqϕ[i=0:nh-1, j=1:nc], qϕ[i, j] == qϕ[i, 0] + δt * sum(rk.a[j, k+1] * dqϕ[i, k] for k in 0:nc-1)
        con_dvx[i=0:nh-1, j=1:nc], vx[i, j] == vx[i, 0] + δt * sum(rk.a[j, k+1] * dvx[i, k] for k in 0:nc-1)
        con_dvy[i=0:nh-1, j=1:nc], vy[i, j] == vy[i, 0] + δt * sum(rk.a[j, k+1] * dvy[i, k] for k in 0:nc-1)

        [i=0:nh-1], qx[i+1, 0] == qx[i, 0] + δt * sum(rk.b[j+1] * dqx[i, j] for j in 0:nc-1)
        [i=0:nh-1], qy[i+1, 0] == qy[i, 0] + δt * sum(rk.b[j+1] * dqy[i, j] for j in 0:nc-1)
        [i=0:nh-1], qϕ[i+1, 0] == qϕ[i, 0] + δt * sum(rk.b[j+1] * dqϕ[i, j] for j in 0:nc-1)
        [i=0:nh-1], vx[i+1, 0] == vx[i, 0] + δt * sum(rk.b[j+1] * dvx[i, j] for j in 0:nc-1)
        [i=0:nh-1], vy[i+1, 0] == vy[i, 0] + δt * sum(rk.b[j+1] * dvy[i, j] for j in 0:nc-1)
    end)

    # Friction
    h(vx, vy) = (vx^2 + vy^2)^0.5
    function ∇h(g::AbstractVector{T}, vx::T, vy::T) where {T}
        if vx == 0 && vy == 0
            g[1] = 1
            g[2] = 0
            return
        end
        g[1] = vx / (vx^2 + vy^2)^0.5
        g[2] = vy / (vx^2 + vy^2)^0.5
        return
    end
    function ∇²h(H::AbstractMatrix{T}, vx::T, vy::T) where {T}
        H[1, 1] = 0
        H[2, 1] = 0
        H[2, 2] = 0
        vx == 0 && vy == 0 && return
        H[1, 1] = (1 - vx^2 / (vx^2+vy^2)) / (vx^2 + vy^2)^0.5
        H[2, 1] = -vx*vy / (vx^2 + vy^2)^(3/2)
        H[2, 2] = (1 - vy^2 / (vx^2+vy^2)) / (vx^2 + vy^2)^0.5
        return
    end
    @operator(model, reg_norm, 2, h, ∇h, ∇²h)

    @constraints(model, begin
        [i=0:nh], -epsilon <= vu[i]*ru[i] + vv[i]*rv[i] + μ * s[i] * g*m <= epsilon
        [i=0:nh], (ru[i]^2 + rv[i]^2) <= (g*m/2 * μ)^2
        #[i=0:nh], s[i] == reg_norm(vx[i, 0], vy[i, 0])
        [i=0:nh], s[i]^2 == vx[i,0]^2 + vy[i,0]^2
        [i=0:nh], s[i] >= 0
    end)


    @expressions(model, begin
        u[i=0:nh],  vx[i, 0]*cos(qϕ[i, 0]) + vy[i, 0]*sin(qϕ[i, 0])
        v[i=0:nh], -vx[i, 0]*sin(qϕ[i, 0]) + vy[i, 0]*cos(qϕ[i, 0])
        r[i=0:nh],  vϕ[i, 0]
    end)

    # Commands
    @constraints(model, begin
        u_ic, u[0] == u_0
        r_ic, r[0] == r_0
        v_ic, v[0] == 0

        u_fc, u[nh] == u_f
        r_fc, r[nh] == r_f
        v_fc, v[nh] == 0
    end)
end

function corner_point_2012_collisions!(model, nh; l1_c = 2.4, l2_c = 0.9, epsilon=1e-6)
    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [0; 0; 0; 0]
    C = [0 -1;    -1 0;     1 0;     0 -1;    0 1;     1 0;    -1 0;     0 1]
    d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
    S(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    R(x, y) = [x ; y]

    @variables(model, begin
        ω[t=0:nh, i=0:3, j=0:5] >= 0, (start = 0)
    end)

    x = model[:x]
    y = model[:y]
    ϕ = model[:ϕ]

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:3], transpose([A * transpose(S(ϕ[t])); C[[2*i+1, 2*i+2], :]]) * ω[t, i, :] == 0
        con_c2[t=1:nh, i=0:3], transpose([b .+ A * transpose(S(ϕ[t])) * R(x[t], y[t]); d[[2*i+1, 2*i+2]]]) * ω[t, i, :] <= -epsilon
    end)
end

function corner_rect_2012_collisions!(model, nh; l1_c = 2.4, l2_c = 0.9, rw_c = 0.72, epsilon_p=1e-6)
    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = [0 -1; -1 0; 1 0; 0 -1; 0 1; 1 0; -1 0; 0 1]
    d = [-l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2; -l2_c/2; -l1_c/2; -l1_c/2; -l2_c/2]
    S(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    R(x, y) = [x ; y]
    epsilon = epsilon_p

    @variables(model, begin ω[t=0:nh, i=0:3, j=0:5] >= 0, (start = 0)  end)

    x = model[:x]
    y = model[:y]
    ϕ = model[:ϕ]

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:3], transpose([A * transpose(S(ϕ[t])); C[[2*i+1, 2*i+2], :]]) * ω[t, i, :] == 0
        con_c2[t=1:nh, i=0:3], transpose([b .+ A * transpose(S(ϕ[t])) * R(x[t], y[t]); d[[2*i+1, 2*i+2]]]) * ω[t, i, :] <= -epsilon
    end)
end

function Npoly_rect_2012_collisions!(model, nh, dims, C_p, d_p; rw_c = 0.72, rh_c = 1.128, epsilon_p=1e-6)
    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = C_p
    d = d_p
    N_poly, N_faces = dims
    S(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    R(x, y) = [x ; y]
    epsilon = epsilon_p

    @variables(model, begin
        ω[t=0:nh, i=0:(N_poly-1), j=0:(4+N_faces-1)] >= 0, (start = 0)
    end)

    x = model[:x]
    y = model[:y]
    ϕ = model[:ϕ]

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose([A * transpose(S(ϕ[t])); C[N_faces*i+1 : N_faces*(i+1), :]]) * ω[t, i, :] == 0
        con_c2[t=1:nh, i=0:(N_poly-1)], transpose([b .+ A * transpose(S(ϕ[t])) * R(x[t], y[t]); d[N_faces*i+1 : N_faces*(i+1)]]) * ω[t, i, :] <= -epsilon
    end)
end

function Npoly_rect_2017_penetration_collisions!(model, nh, dims, C_p, d_p; kappa=10)
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot

    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = C_p
    d = d_p
    N_poly, N_faces = dims
    R(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    tr(x, y) = [x ; y]

    @variables(model, begin
        s[t=0:nh, i=0:(N_poly-1)]                  >= 0, (start = 0)
        μ[t=0:nh, i=0:(N_poly-1), j=0:3]           >= 0, (start = 0)
        λ[t=0:nh, i=0:(N_poly-1), j=0:(N_faces-1)] >= 0, (start = 0)
    end)

    T = model[:T]
    #=v = model[:v]
    e_t=model[:e_t]
    ϕdot=model[:ϕdot]
    F = model[:F]=#

    @objective(model, Min, T + kappa * sum(s))# + (1e-3)*sum(abs.((transpose.(v) .* e_t) .* ϕdot .- F)))

    x = model[:x]
    y = model[:y]
    ϕ = model[:ϕ]

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose(transpose(C[N_faces*i+1 : N_faces*(i+1), :]) * λ[t, i, :]) * (transpose(C[N_faces*i+1 : N_faces*(i+1), :]) * λ[t, i, :]) == 1
        con_c2[t=1:nh, i=0:(N_poly-1)], transpose(A) * μ[t, i, :] + transpose(C[N_faces*i+1 : N_faces*(i+1), :] * R(ϕ[t])) * λ[t, i, :] == 0
        con_c3[t=1:nh, i=0:(N_poly-1)], -transpose(b) * μ[t, i, :] + transpose(C[N_faces*i+1 : N_faces*(i+1), :] * tr(x[t], y[t]) - d[N_faces*i+1 : N_faces*(i+1)]) * λ[t, i, :] >= -s[t, i]
    end)
end

function Npoly_rect_2017_collisions!(model, nh, dims, C_p, d_p; d_min=1e-4)
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot

    # contraintes de collisions
    A = [1 0; 0 1; -1 0; 0 -1]
    b = [rh_c/2; rw_c/2; rh_c/2; rw_c/2]
    C = C_p
    d = d_p
    N_poly, N_faces = dims
    R(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    tr(x, y) = [x ; y]

    @variables(model, begin
        μ[t=0:nh, i=0:(N_poly-1), j=0:3]           >= 0, (start = 0)
        λ[t=0:nh, i=0:(N_poly-1), j=0:(N_faces-1)] >= 0, (start = 0)
    end)

    x = model[:x]
    y = model[:y]
    ϕ = model[:ϕ]

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose(transpose(C[N_faces*i+1 : N_faces*(i+1), :]) * λ[t, i, :]) * (transpose(C[N_faces*i+1 : N_faces*(i+1), :]) * λ[t, i, :]) <= 1
        con_c2[t=1:nh, i=0:(N_poly-1)], transpose(A) * μ[t, i, :] + transpose(C[N_faces*i+1 : N_faces*(i+1), :] * R(ϕ[t])) * λ[t, i, :] == 0
        con_c3[t=1:nh, i=0:(N_poly-1)], -transpose(b) * μ[t, i, :] + transpose(C[N_faces*i+1 : N_faces*(i+1), :] * tr(x[t], y[t]) - d[N_faces*i+1 : N_faces*(i+1)]) * λ[t, i, :] >= d_min
    end)
end

function Npoly_rect_2023_collisions!(model, nh, poly; d_min=1e-4)
    rw_c = 0.72 # largeur du robot
    rh_c = 1.128 # longueur du robot

    # contraintes de Collisions
    N_poly = length(poly)
    Ve = [(-rh_c/2) (rh_c/2) (rh_c/2) (-rh_c/2); (rw_c/2) (rw_c/2) (-rw_c/2) (-rw_c/2)]
    Vo = [reduce(hcat, vertices) for vertices in poly]

    R(ϕ) = [cos(ϕ) -sin(ϕ); sin(ϕ) cos(ϕ)]
    tr(x, y) = [x ; y]

    @variables(model, begin
        ξ[t=0:nh, i=0:(N_poly-1), j=0:1], (start = 0.0)
        μ_r[t=0:nh, i=0:(N_poly-1)], (start = 0.0)
        μ_o[t=0:nh, i=0:(N_poly-1)], (start = 0.0)
    end)
    
    x = model[:x]
    y = model[:y]
    ϕ = model[:ϕ]

    # Collisions
    @constraints(model, begin
        con_c1[t=1:nh, i=0:(N_poly-1)], transpose(ξ[t, i, :]) * ξ[t, i, :] / 4 + μ_r[t, i] + μ_o[t, i] + d_min^2 <= 0
        con_c2[t=1:nh, i=0:(N_poly-1)], -transpose(R(ϕ[t]) * Ve) * ξ[t, i, :] .- transpose(tr(x[t], y[t]))*ξ[t, i, :] .- μ_r[t, i] <= 0
        con_c3[t=1:nh, i=0:(N_poly-1)], transpose(Vo[i+1]) * ξ[t, i, :] .- μ_o[t, i] <= 0
    end)
end



function solve!(model; max_iter=500, verbose=5)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "max_iter", max_iter)
    JuMP.set_optimizer_attribute(model, "print_level", verbose)
    JuMP.optimize!(model)
end
using JuMP
using Ipopt
using LinearAlgebra

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

function basic_model!(model, nh, cdt_0, cdt_f; u_max_p=2.75)
    u_max = u_max_p # vitesse maximale d'une roue

    # conditions initiales
    x_0, y_0, ϕ_0, u_0, r_0 = cdt_0

    # conditions finales
    x_f, y_f, ϕ_f, u_f, r_f = cdt_f

    # etat initial
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_straight(nh, [x_0, y_0], [x_f, y_f], u_max/2)
    step = 1 / nh

    @variables(model, begin
        x[i=0:nh], (start = x_i[i+1])
        y[i=0:nh], (start = y_i[i+1])
        ϕ[i=0:nh], (start = ϕ_i[i+1])
        u[i=0:nh], (start = u_i[i+1])
        r[i=0:nh], (start = r_i[i+1])
        0.0 <= T,  (start = T_i)
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

    return x_i, y_i, ϕ_i, u_i, r_i, T_i
end

function speed_model!(model, nh, cdt_0, cdt_f; u_max=2.7, rw_c=0.72)
    basic_model!(model, nh, cdt_0, cdt_f; u_max_p=u_max)

    d_c = rw_c*0.4;

    u = model[:u]
    r = model[:r]

    # Dynamics
    @constraints(model, begin
        con_s1[i=1:nh], -u_max <= u[i] + r[i]*d_c <= u_max # cinématique
        con_s2[i=1:nh], -u_max <= u[i] - r[i]*d_c <= u_max # cinématique
    end)
end

function dynamic_model!(model, nh, cdt_0, cdt_f; U_max=48, rw=0.72, rh=1.128, m=105, epsilon = 1e-4)
    Izz = m * (rh^2 + rw^2) / 12
    ρ = 0.10;
    d = rw*0.4; # espacement des roues
    l = 0.1; # hauteur du centre de gravit au-dessus de l'axe des roues
    b = rh*0.4 # distance des roulettes
    f = 0.5 # coefficient de frottements secs

    R = 1.43; # résistance du moteur
    k = 0.107; # constante de couple

    Meq = m # masse équivalente du robot avec ses roues
    Ieq = Izz # intertie équivalente du robot avec ses roues
    g = 9.81

    _, _, _, u0, r0, _ = basic_model!(model, nh, cdt_0, cdt_f; u_max_p=2.7)


    u = model[:u]
    r = model[:r]
    δt = model[:δt]

    @variables(model, begin
        -U_max <= ul[i=0:nh] <= U_max, (start = k * (u0[i+1] - d*r0[i+1]) / ρ)
        -U_max <= ur[i=0:nh] <= U_max, (start = k * (u0[i+1] - d*r0[i+1]) / ρ)
        
        # forces (pour les limites de frottements)
        δN[i=0:nh] >= 0, (start = 0)
        τplus[i=0:nh] >= 0, (start = 0)
        τmoins[i=0:nh] >= 0, (start = 0)
        Tl[i=0:nh], (start = 0)
        Tr[i=0:nh], (start = 0)
    end)

    @expressions(model, begin
        τl[i=0:nh], k * (ul[i] - k/ρ*(u[i] - d*r[i])) / R
        τr[i=0:nh], k * (ur[i] - k/ρ*(u[i] + d*r[i])) / R

        # forces
        τ[i=0:nh], τl[i] + τr[i] 
        nl[i=0:nh], -((l+ρ) * r[i]*u[i] - d*Meq*g) / (2*d)
        nr[i=0:nh],  ((l+ρ) * r[i]*u[i] + d*Meq*g) / (2*d)
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
        con_comp_T[i=0:nh],  Tl[i] + Tr[i] == -r[i]*u[i]
        con_comp_N1[i=0:nh], Tl[i]^2 + (τl[i]/ρ)^2 <= (f * (nl[i] - (2*τplus[i]-τ[i]) / b))^2
        con_comp_N2[i=0:nh], Tl[i]^2 + (τl[i]/ρ)^2 <= (f * (nl[i] - (2*τmoins[i]+τ[i]) / b))^2
        con_comp_N3[i=0:nh], Tr[i]^2 + (τr[i]/ρ)^2 <= (f * (nr[i] - (2*τplus[i]-τ[i]) / b))^2
        con_comp_N4[i=0:nh], Tr[i]^2 + (τr[i]/ρ)^2 <= (f * (nr[i] - (2*τmoins[i]+τ[i]) / b))^2
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

function Npoly_rect_2017_collisions!(model, nh, dims, C_p, d_p; kappa=10)
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
    @objective(model, Min, T + kappa * sum(s))

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
        con_c2[t=1:nh, i=0:(N_poly-1)], -transpose(R(ϕ[t]) * Ve .+ tr(x[t], y[t])) * ξ[t, i, :] .- μ_r[t, i] <= 0
        con_c3[t=1:nh, i=0:(N_poly-1)], transpose(Vo[i+1]) * ξ[t, i, :] .- μ_o[t, i] <= 0
    end)
end


function solve!(model; max_iter=500)
    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "max_iter", max_iter)
    JuMP.optimize!(model)
end
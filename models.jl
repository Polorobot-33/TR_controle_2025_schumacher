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
    x_i, y_i, ϕ_i, u_i, r_i, T_i = init_arc(nh, -3, 3, u_max/2, 1)
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
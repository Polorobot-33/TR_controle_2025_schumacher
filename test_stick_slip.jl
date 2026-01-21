using JuMP
using CairoMakie


# simulation params
Tf = 5
nh = 100
epsilon = 1e-4

v_ref = 1
x_ref = LinRange(0, Tf * v_ref, nh+1)
δt = Tf / nh


# physical params
m = 1
k = 10
μ = 0.5

model = Model()

@variables(model, begin
    x[i=0:nh], (start = x_ref[i+1])
    v[i=0:nh], (start = v_ref)

    xplus[i=0:nh] >= 0, (start = max(v_ref, 0))
    xmoins[i=0:nh] >= 0, (start = max(0, -v_ref))
    -1 <= sgn[i=0:nh] <= 1, (start = sign(v_ref))
end)

@expressions(model, begin
    f[i=0:nh], k*(x_ref[i+1]-x[i]) - μ*9.81*m*sgn[i]
end)

@objective(model, Min, 1)

# Dynamics
@constraints(model, begin
    con_x[i=1:nh], x[i] == x[i-1] + δt * (v[i-1] + v[i]) / 2
    con_v[i=1:nh], v[i] == v[i-1] + δt * (f[i] + f[i-1]) / (2 * m)
    
    # frottements
    con_xsgn[i=0:nh], xplus[i] - xmoins[i] == v[i]
    con_orth[i=0:nh], xplus[i] * xmoins[i] <= epsilon
    con_sgn1[i=0:nh], (1 - sgn[i]) * xplus[i] <= epsilon
    con_sgn2[i=0:nh], (1 + sgn[i]) * xmoins[i] <= epsilon
end)

# Boundary constraints
@constraints(model, begin
    x_ic, x[0] == x_ref[begin]
    v_ic, v[0] == 0
end)

JuMP.set_optimizer(model, Ipopt.Optimizer)
JuMP.set_optimizer_attribute(model, "max_iter", 1000)
JuMP.optimize!(model)

function var(model, val)
    return JuMP.value.(model[val]).data
end


# rendering
f = Figure(size = (560, 920))#(512, 920))
times = LinRange(0, Tf, nh+1)


ax_x = CairoMakie.Axis(f[1, :], xlabel="tps (s)", ylabel="position (m)", alignmode=CairoMakie.Inside())
lines!(ax_x, times, var(model, :x))
lines!(ax_x, times, x_ref)

ax_v = CairoMakie.Axis(f[2, :], xlabel="tps (s)", ylabel="vitesse (m/s)", alignmode=CairoMakie.Inside())
lines!(ax_v, times, var(model, :v))

ax_f = CairoMakie.Axis(f[3, :], xlabel="tps (s)", ylabel="force (N)", alignmode=CairoMakie.Inside())
lines!(ax_f, times, var(model, :f))
lines!(ax_f, times, var(model, :sgn))

f

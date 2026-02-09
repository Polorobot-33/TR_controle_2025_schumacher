using JuMP
using Ipopt
using LinearAlgebra
using CairoMakie

include("../collocation.jl")

# simulation params
Tf = 5
nh = 99
rk = CrankNicolson()
nc = length(rk.c)

epsilon = 1e-3


δt = Tf / nh

# physical params
m = 1.0
g = 9.81
μ = 0.5

Fx(i) = max(0, min(i/nh * 14 - 1, 7))
Fy(i) = max(i/nh * 12 - 9, 0)


function solve_inner_problem(sol)
    qx_0, qy_0, vx_0, vy_0, rx_0, ry_0, s = sol

    model = Model()
    @objective(model, Min, 1)

    @variables(model, begin
        qx[i=0:nh, j=0:nc], (start = qx_0[i+1, j+1])
        qy[i=0:nh, j=0:nc], (start = qy_0[i+1, j+1])

        vx[i=0:nh, j=0:nc], (start = vx_0[i+1, j+1])
        vy[i=0:nh, j=0:nc], (start = vy_0[i+1, j+1])

        rx[i=0:nh], (start = rx_0[i+1])
        ry[i=0:nh], (start = ry_0[i+1])
    end)

    @expressions(model, begin
        dqx[i=0:nh-1, j=0:nc-1], vx[i, j]
        dqy[i=0:nh-1, j=0:nc-1], vy[i, j]

        dvx[i=0:nh-1, j=0:nc-1], (rx[i] + Fx(i)) / m
        dvy[i=0:nh-1, j=0:nc-1], (ry[i] + Fy(i)) / m
    end)

    # Dynamics
    @constraints(model, begin
        con_dqx[i=0:nh-1, j=1:nc], qx[i, j] == qx[i, 0] + δt * sum(rk.a[j, k+1] * dqx[i, k] for k in 0:nc-1)
        con_dqy[i=0:nh-1, j=1:nc], qy[i, j] == qy[i, 0] + δt * sum(rk.a[j, k+1] * dqy[i, k] for k in 0:nc-1)
        con_dvx[i=0:nh-1, j=1:nc], vx[i, j] == vx[i, 0] + δt * sum(rk.a[j, k+1] * dvx[i, k] for k in 0:nc-1)
        con_dvy[i=0:nh-1, j=1:nc], vy[i, j] == vy[i, 0] + δt * sum(rk.a[j, k+1] * dvy[i, k] for k in 0:nc-1)

        [i=0:nh-1], qx[i+1, 0] == qx[i, 0] + δt * sum(rk.b[j+1] * dqx[i, j] for j in 0:nc-1)
        [i=0:nh-1], qy[i+1, 0] == qy[i, 0] + δt * sum(rk.b[j+1] * dqy[i, j] for j in 0:nc-1)
        [i=0:nh-1], vx[i+1, 0] == vx[i, 0] + δt * sum(rk.b[j+1] * dvx[i, j] for j in 0:nc-1)
        [i=0:nh-1], vy[i+1, 0] == vy[i, 0] + δt * sum(rk.b[j+1] * dvy[i, j] for j in 0:nc-1)
    end)

    # Friction
    @expression(model, vTr[i=0:nh], vx[i, 0]*rx[i] + vy[i, 0]*ry[i])
    @constraints(model, begin
        [i=0:nh], -epsilon <= vTr[i] + μ*s[i+1] * g*m <= epsilon
        [i=0:nh], (rx[i]^2 + ry[i]^2) <= (g*m*μ)^2
    end)

    fix(qx[0, 0], 0.0; force=true)
    fix(qy[0, 0], 0.0; force=true)
    fix(vx[0, 0], 0.0; force=true)
    fix(vy[0, 0], 0.0; force=true)

    JuMP.set_optimizer(model, Ipopt.Optimizer)
    JuMP.set_optimizer_attribute(model, "max_iter", 100)
    JuMP.set_optimizer_attribute(model, "print_level", 4)
    JuMP.optimize!(model)

    return  value.(model[:qx]).data, value.(model[:qy]).data,
            value.(model[:vx]).data, value.(model[:vy]).data,
            value.(model[:rx]).data, value.(model[:ry]).data,
            s
end

function solve_outer_problem(sol)
    qx, qy, vx, vy, rx, ry, _ = sol
    return qx, qy, vx, vy, rx, ry, [norm(vx[i+1, 1], vy[i+1, 1]) for i in 0:nh]
end

max_iter = 10
qx_0 = zeros(Float64, nh+1, nc+1)
qy_0 = zeros(Float64, nh+1, nc+1)
vx_0 = zeros(Float64, nh+1, nc+1)
vy_0 = zeros(Float64, nh+1, nc+1)
rx_0 = zeros(Float64, nh+1)
ry_0 = zeros(Float64, nh+1)
s_0  = zeros(Float64, nh+1)
sol = qx_0, qy_0, vx_0, vy_0, rx_0, ry_0, s_0
for i in 1:max_iter
    sol1 = solve_inner_problem(sol)
    global sol = solve_outer_problem(sol1)

    println(i / max_iter * 100)
end

function var(model, val)
    return JuMP.value.(model[val]).data
end
_, _, vx, vy, rx, ry, s = sol


# rendering
f = Figure(size = (560, 920))#(512, 920))
times = LinRange(0, Tf, nh+1)


ax_v = CairoMakie.Axis(f[1, :], xlabel="tps (s)", ylabel="vitesse (m)", alignmode=CairoMakie.Inside())
#lines!(ax_v, times, var(model, :v), label="||v||")
lines!(ax_v, times, vx[:, 1], label="vx")
lines!(ax_v, times, vy[:, 1], label="vy")
axislegend(ax_v)

ax_f = CairoMakie.Axis(f[2, :], xlabel="tps (s)", ylabel="force (N)", alignmode=CairoMakie.Inside())
lines!(ax_f, times, [Fx(i) for i in 0:nh], label="Fx")
lines!(ax_f, times, [Fy(i) for i in 0:nh], label="Fy")
lines!(ax_f, times, rx, label="rx")
lines!(ax_f, times, ry, label="ry")
lines!(ax_f, times, s, label="s")
axislegend(ax_f)

f

using JuMP
using Ipopt
using LinearAlgebra
using CairoMakie

include("../collocation.jl")

# simulation params
Tf = 5
nh = 99
rk = ImplicitEuler()
nc = length(rk.c)

epsilon = 1e-3


δt = Tf / nh

# physical params
m = 1.0
g = 9.81
μ = 0.5

Fx(i) = max(0, min(i/nh * 14 - 1, 7))
Fy(i) = max(i/nh * 12 - 9, 0)

model = Model()

@objective(model, Min, 1)

@variables(model, begin
    qx[i=0:nh, j=0:nc], (start = 0)
    qy[i=0:nh, j=0:nc], (start = 0)

    vx[i=0:nh, j=0:nc], (start = 0)
    vy[i=0:nh, j=0:nc], (start = 0)

    s[i=0:nh],  (start=0)
    rx[i=0:nh], (start=0)
    ry[i=0:nh], (start=0)
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
    vx == 0 && vy == 0 && return
    H[1, 1] = (1 - vx^2 / (vx^2+vy^2)) / (vx^2 + vy^2)^0.5
    H[2, 1] = -vx*vy / (vx^2 + vy^2)^(3/2)
    H[2, 2] = (1 - vy^2 / (vx^2+vy^2)) / (vx^2 + vy^2)^0.5
    return
end
@operator(model, reg_norm, 2, h, ∇h, ∇²h)

# Expressions
@expressions(model, begin
    v2[i=0:nh], vx[i, 0]^2 + vy[i, 0]^2
    v[i=0:nh], v2[i] ^ 0.5
    rt[i=0:nh], (rx[i]^2 + ry[i]^2) ^ 0.5
end)

# Friction
@expression(model, vTr[i=0:nh], vx[i, 0]*rx[i] + vy[i, 0]*ry[i])
@constraints(model, begin
    [i=0:nh], -epsilon <= vTr[i] + μ*s[i] * g*m <= epsilon
    [i=0:nh], (rx[i]^2 + ry[i]^2) <= (g*m*μ)^2
    
    [i=0:nh], s[i] == reg_norm(vx[i, 0], vy[i, 0])
    #=[i=0:nh], s[i]^2 == v2[i]
    [i=0:nh], s[i] >= 0=#
end)


# Boundary constraints
@constraints(model, begin
    qx[0, 0] == 0
    qy[0, 0] == 0
    vx[0, 0] == 0
    vy[0, 0] == 0
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


ax_v = CairoMakie.Axis(f[1, :], xlabel="tps (s)", ylabel="vitesse (m)", alignmode=CairoMakie.Inside())
lines!(ax_v, times, var(model, :v), label="||v||")
lines!(ax_v, times, var(model, :vx)[:, 1], label="vx")
lines!(ax_v, times, var(model, :vy)[:, 1], label="vy")
axislegend(ax_v)

ax_f = CairoMakie.Axis(f[2, :], xlabel="tps (s)", ylabel="force (N)", alignmode=CairoMakie.Inside())
lines!(ax_f, times, [Fx(i) for i in 0:nh], label="Fx")
lines!(ax_f, times, [Fy(i) for i in 0:nh], label="Fy")
lines!(ax_f, times, var(model, :rx), label="rx")
lines!(ax_f, times, var(model, :ry), label="ry")
lines!(ax_f, times, var(model, :rt), label="||rt||")
axislegend(ax_f)

f

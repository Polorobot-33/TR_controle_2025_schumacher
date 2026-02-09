using JuMP
using Ipopt # Clarabel # 
using LinearAlgebra
using CairoMakie

include("../collocation.jl")

# simulation params
Tf = 5
nh = 99

epsilon = 1e-6
kappa = 1e+6

δt = Tf / nh

# physical params
m = 1.0
g = 9.81
μ = 0.5

Fx(i) = max(0, min(i/nh * 14 - 1, 7))
Fy(i) = max(i/nh * 12 - 9, 0)

# conditions initiales
qx_s, qy_s, vx_s, vy_s, rx_s, ry_s = [0.], [0.], [0.], [0.], [0.], [0.]
s_s = [0.]

for i in 1:nh
    println("progress : $(trunc(i / nh * 100, digits=1, base=10)) %");

    s_ref = Ref(0.0);
    model = Model()

    @variables(model, begin
        qx, (start = qx_s[end])
        qy, (start = qy_s[end])

        vx, (start = vx_s[end])
        vy, (start = vy_s[end])
        s,  (start = sqrt(vx_s[end]^2 + vy_s[end]^2))

        rx, (start = rx_s[end])
        ry, (start = ry_s[end])
    end)

    #=f(p...) = (p[1] - (p[2]^2 + p[3]^2)^0.5)^2
    function ∇f(g::AbstractVector{T}, p::T...) where {T}
        g[1] = 2 * (p[1] - (p[2]^2 + p[3]^2)^0.5)
        g[2] = p[2] == 0 ? 0 : 2*p[2] * (1 - p[1] / (p[2]^2 + p[3]^2)^0.5)
        g[3] = p[3] == 0 ? 0 : 2*p[3] * (1 - p[1] / (p[2]^2 + p[3]^2)^0.5)
        return
    end
    #=∇²f(H::AbstractMatrix{T}, x::T...) where {T}
        ...
        return
    end=#
    @operator(model, distObj, 3, f, ∇f)=#
    
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


    @objective(model, Min, 1)#kappa * distObj(s, vx, vy))

    # Dynamics
    @constraints(model, begin
        # Crank Nicolson
        qx == qx_s[end] + δt * (vx + vx_s[end]) / 2
        qy == qy_s[end] + δt * (vy + vy_s[end]) / 2
        vx == vx_s[end] + δt * (rx + rx_s[end] + Fx(i-1) + Fx(i)) / (2*m)
        vy == vy_s[end] + δt * (ry + ry_s[end] + Fy(i-1) + Fy(i)) / (2*m)

        #=# Implicit Euler
        qx == qx_s[end] + δt * vx
        qy == qy_s[end] + δt * vy
        vx == vx_s[end] + δt * (rx + Fx(i)) / m
        vy == vy_s[end] + δt * (ry + Fy(i)) / m=#

        # ~Explicit Euler
        #=qx == qx_s[end] + δt * vx_s[end]
        qy == qy_s[end] + δt * vy_s[end]
        vx == vx_s[end] + δt * (rx + Fx(i)) / m
        vy == vy_s[end] + δt * (ry + Fy(i)) / m=#
    end)

    # Friction
    @expression(model, vTr, vx*rx + vy*ry)

    @constraints(model, begin
        -epsilon <= vTr + μ * s * g*m <= epsilon
        #[g*m*μ, rx, ry] in SecondOrderCone()
        (rx^2 + ry^2) <= (g*m*μ)^2

        #s >= 0
        #s^2 == (vx^2 + vy^2)
        #[s, vx, vy] in SecondOrderCone()
        s == reg_norm(vx, vy)
    end)

    #s_ref[] = sqrt(value(vx)^2 + value(vy)^2)

    JuMP.set_optimizer(model, Ipopt.Optimizer)# Clarabel.Optimizer) #
    JuMP.set_optimizer_attribute(model, "max_iter", 100)
    JuMP.set_optimizer_attribute(model, "print_level", 0)
    JuMP.optimize!(model)

    push!(qx_s, JuMP.value.(model[:qx]))
    push!(qy_s, JuMP.value.(model[:qy]))
    push!(vx_s, JuMP.value.(model[:vx]))
    push!(vy_s, JuMP.value.(model[:vy]))
    push!(rx_s, JuMP.value.(model[:rx]))
    push!(ry_s, JuMP.value.(model[:ry]))

    push!(s_s, JuMP.value.(model[:s]))
end


# rendering
f = Figure(size = (560, 920))#(512, 920))
times = LinRange(0, Tf, nh+1)


ax_v = CairoMakie.Axis(f[1, :], xlabel="tps (s)", ylabel="vitesse (m)", alignmode=CairoMakie.Inside())
lines!(ax_v, times, ((x, y) -> sqrt(x^2 + y^2)).(vx_s, vy_s), label="||v||")
lines!(ax_v, times, s_s, label="s")
lines!(ax_v, times, vx_s, label="vx")
lines!(ax_v, times, vy_s, label="vy")
axislegend(ax_v, position=:lt)

ax_f = CairoMakie.Axis(f[2, :], xlabel="tps (s)", ylabel="force (N)", alignmode=CairoMakie.Inside())
lines!(ax_f, times, [Fx(i) for i in 0:nh], label="Fx")
lines!(ax_f, times, [Fy(i) for i in 0:nh], label="Fy")
lines!(ax_f, times, rx_s, label="rx")
lines!(ax_f, times, ry_s, label="ry")
lines!(ax_f, times, ((x, y) -> sqrt(x^2 + y^2)).(rx_s, ry_s), label="||rt||")
axislegend(ax_f, position=:lt)

f

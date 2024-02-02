using Pkg 
Pkg.add("Plots")

using Pkg 
Pkg.add("DifferentialEquations")

using Plots
using DifferentialEquations

function forbit(U::Vector)
    x, y, vx, vy = U
    xdot = vx
    ydot = vy
    vxdot = -x/(x^2 + y^2)^(3/2)
    vydot = -y/(x^2 + y^2)^(3/2)
    Udot = [xdot, ydot, vxdot, vydot]
    return Udot 
end 

function euler(f, y, h)
    k = f(y)
    y = y + h * k
    return y
end

function rk2(f, y, h)
    k = f(y)
    y = y + h/2 * k
    k = f(y)
    y = y + h * k
    return y 
end 

N = 100

h = 0.01
U0 = [1.0, 0.0, 0, 1.0]
Uall_euler = [U0]
Uall_rk2 = [U0]

for _ in 1:N
    U_enext = euler(forbit, last(Uall_euler), h)
    push!(Uall_euler, U_enext)

    U_rnext = rk2(forbit, last(Uall_rk2), h)
    push!(Uall_rk2, U_rnext)
end


X_euler = []
Y_euler = []

X_rk2 = []
Y_rk2 = []

for i in 1:N
    push!(X_euler, Uall_euler[i][1])
    push!(Y_euler, Uall_euler[i][2])
    push!(X_rk2, Uall_rk2[i][1])
    push!(Y_rk2, Uall_rk2[i][2])
end

prob = ODEProblem((U,p,t) -> forbit(U), U0, (0.0, 10.0))
sol = solve(prob)

sol_ODEP = collect(eachrow(reduce(hcat, sol.u)))

plot()
plot!(X_euler, Y_euler)
plot!(X_rk2, Y_rk2)
plot!(sol_ODEP[1], sol_ODEP[2])

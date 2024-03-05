using Pkg 
Pkg.add("DifferentialEquations")

# Stability : error changes in uncontrolled ways

# Harmonic oscillator: x, v

# Can numerical error model friction (can you make the error simulate friction?)

function fho(U::Vector)
    x, v = U
    xdot = v
    ydot = -x
    Udot = [xdot, ydot]
    return Udot 
end 

function euler(f, y, h)
    k = f(y)
    y = y + h*k
    return y
end


U = [1.0, 0.0]
h = 0.1
Uall = [U]

for n in 1:10
    U = euler(fho, U, h)
    push!(Uall, U)
end

Uall


using DifferentialEquations
U0 = [1.0, 0.0]
prob = ODEProblem((U, p, t) -> fho(U), U0, (0, 3.14))

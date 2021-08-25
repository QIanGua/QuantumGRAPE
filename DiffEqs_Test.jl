using DifferentialEquations
using Plots

# Example 1: lorzen curve
function lorenz!(du,u,p,t)
    du[1] = 10.0*(u[2]-u[1])
    du[2] = u[1]*(28.0-u[3]) - u[2]
    du[3] = u[1]*u[2] - (8/3)*u[3]
end

u0 = [1.0;0.0;0.0]
tspan = (0.0,100.0)
prob = ODEProblem(lorenz!,u0,tspan)
sol = solve(prob)

plot(sol,vars=(1,2,3))


# Example 2: circle
# dx/dt = - y
# dy/dt = x
# x0 = 1.0
# y0 = 0.0
function circle!(du,u,p,t)
    du[1] = -u[2]
    du[2] = u[1]
end

u0 = [1.0;0.0]
tspan = (0.0,2*pi)
prob = ODEProblem(circle!,u0,tspan)
sol = solve(prob,RK4())

plot(sol,vars = (1,2))

# matrix form
A = [0.0 -1.0
     1.0 0.0]
f2(u,p,t) = A*u
u0 = [1.0;0.0]
tspan = (0.0,2*pi)
prob2 = ODEProblem(f2,u0,tspan)
sol2 = solve(prob2,RK4())
plot(sol2,vars = (1,2))

using SparseArrays
f3(u,p,t) = sparse(A)*u
prob3 = ODEProblem(f3,u0,tspan)
sol3 = solve(prob3,RK4())

plot(sol3,vars = (1,2))
# value check
sol.t[1]
sol.u[1][1:2]

sol2.t[1]
sol2.u[1][1:2]
sol.t
sol2.t
using Pkg
Pkg.activate(".")

using DifferentialEquations

function derivatives!(du, u, p, t)
  du[1] = -u[1]
end

u0 = [1.0]
tspan = (0.0, 1e0)
dt=1e-2
prob = ODEProblem(derivatives!, u0, tspan)

observations = Float64[]
function condition(u, t, integrator)
  push!(observations, u[1])
  return false
end
cb = DiscreteCallback(condition, nothing)

sol = solve(prob, Tsit5(), dt=dt, adaptive=false, callback=cb, save_everystep=true)
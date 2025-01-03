using Pkg
Pkg.activate(".")

using DifferentialEquations
using GLMakie

### Thomas system.

function thomas_derivatives!(du, u, p, t)
  x, y, z = u
  b = p[1]
  du[1] = sin(y) - b * x
  du[2] = sin(z) - b * y
  du[3] = sin(x) - b * z
end

u0 = [1.0, 1.0, 0.5]
p = [0.208186]
tspan = (0.0, 1e6)

prob = ODEProblem(thomas_derivatives!, u0, tspan, p)

### Lorenz system (debug).

# function lorenz_derivatives!(du, u, p, t)
#   x, y, z = u
#   σ = p[1]
#   ρ = p[2]
#   β = p[3]
#   du[1] = σ * (y - x)
#   du[2] = x * (ρ - z) - y
#   du[3] = x * y - β * z
# end

# u0 = [10.0, 0.0, 20.0]
# p = [10.0, 28.0, 8/3]
# tspan = (0.0, 2e3)
# transient_time = 1e2

# prob = ODEProblem(lorenz_derivatives!, u0, tspan, p)

### Specify Poincare section as a terminator callback.
# Define equilibrium points.
# eq1 = Point3f(2.575647587674765, 2.575647587674765, 2.575647587674765)
# eq2 = Point3f(2.459147755309479, 0.5374650978896967, 3.0294651370145176)
# condition(u, t, integrator) = p[1] * u[3] - sin(u[1])
# affect!(integrator) = terminate!(integrator)
# cb = ContinuousCallback(condition, affect!, nothing)
###

# Solve orbit and discard transient.
@time begin
  sol_full = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10, saveat=1e-2)#, callback=cb)
  # Get index of first point after transient_time
  start_idx = findfirst(t -> t >= transient_time, sol_full.t)
  # Create new solution with transient discarded
  sol = DiffEqBase.build_solution(
    prob,
    sol_full.alg,
    sol_full.t[start_idx:end],
    sol_full.u[start_idx:end]
  )
end

# Compute the weighted ordinal pattern entropies.
include("../return_maps/opn.jl")
# OPN parameters for original timeseries.
w = 3
m = 4
τ = 2
# OPN parameters for each ordinal symbol's constructed timeseries.
w2 = 1
m2 = 3
τ2 = 1
xs = [u[3] for u in sol.u]
ordinal_symbols = encode_chronological_timeseries(xs, w, m, τ)
weighted_entropies = Float64[]
@time for ordinal_symbol in 1:factorial(m)
  ordinal_timeseries = get_ordinal_timeseries(
    xs,
    ordinal_symbols,
    ordinal_symbol,
    w
  )
  if length(ordinal_timeseries) == 0
    push!(weighted_entropies, 0)
  else
    ordinal_symbols2 = encode_chronological_timeseries(
      ordinal_timeseries,
      w2,
      m2,
      τ2
    )
    push!(weighted_entropies, weighted_entropy(length(xs), m, ordinal_symbols2, m2))
  end
end
println("Weighted entropies (sorted): $(reverse(sort(weighted_entropies)))")
ranked_ordinal_symbol_index = 2
ordinal_symbol = sortperm(weighted_entropies, rev=true)[ranked_ordinal_symbol_index]
# println("Highest entropy ordinal symbol: $highest_entropy_ordinal_symbol")

# Get all first-of-consecutive occurrences
# of the highest entropy ordinal symbol.
section_indices = Int[]
for i in 2:length(ordinal_symbols)
  if ordinal_symbols[i] == ordinal_symbol && ordinal_symbols[i-1] != ordinal_symbol
    push!(section_indices, (i - 1)*w + 1)
  end
end

# Get the points in the section.
section_points = [
  Point3f(sol.u[i][1], sol.u[i][2], sol.u[i][3])
  for i in section_indices
]

# Create 3D plot of the solution.
begin
  fig = Figure()
  ax = Axis3(fig[1, 1], 
      xlabel="x", ylabel="y", zlabel="z",
      title="Thomas Attractor")

  # Extract solution points.
  points = [Point3f(u[1], u[2], u[3]) for u in sol.u]

  # Create line plot.
  lines!(ax, points, color=:blue, linewidth=0.5)

  # Draw a line between two equilibria.
  # lines!(ax, [eq1, eq2], color=:red, linewidth=2)

  # Plot points of section.
  scatter!(ax, section_points, color=:red, markersize=10)

  # Display the figure.
  display(fig)

  # Optionally save the final figure.
  # save("thomas_attractor.png", fig)
end
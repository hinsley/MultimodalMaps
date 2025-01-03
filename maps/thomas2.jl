# Constructs a return map for the Thomas attractor with
# spherical section.
# This method is inferior to the OPN method in thomas.jl
# because it requires manually tuning the radius of the
# section and produces a less interpretable return map.

using Pkg
Pkg.activate(".")

using DifferentialEquations
using LinearAlgebra
using Plots

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
transient_time = 1e2

prob = ODEProblem(thomas_derivatives!, u0, tspan, p)

# Track locations of sphere exits.
radius = 4.47
sphere_exits = []
condition(u, _, integrator) = norm(u) - radius
function affect!(integrator)
  if integrator.t > transient_time
    push!(sphere_exits, collect(integrator.u))
  end
end
cb = ContinuousCallback(condition, affect!, nothing)

# Solve orbit and discard transient.
@time begin
  sol_full = solve(prob, Tsit5(), reltol=1e-10, abstol=1e-10, saveat=1e-2, callback=cb)
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

# The timeseries function used for the coordinate in
# the return map.
function return_map_function(u)
  return u[3]
end

# Construct the return map graph.
return_map_points = return_map_function.(sphere_exits)
section_iterate = 3
section_iterate_offset = 1 # Should be from 1 to section_iterate.
return_map_points = [
  return_map_points[i]
  for i in 1:section_iterate:length(return_map_points)
]
return_map_endpoints = [
  minimum(return_map_points),
  maximum(return_map_points)
]
return_map_xs = return_map_points[1:end-1]
return_map_ys = return_map_points[2:end]
cobweb_xs = [return_map_xs[1]]
cobweb_ys = [return_map_xs[1]]
for i in 1:length(return_map_xs)-1
    push!(cobweb_xs, cobweb_xs[end], return_map_xs[i+1])
    push!(cobweb_ys, return_map_ys[i], return_map_ys[i])
end
push!(cobweb_xs, cobweb_xs[end])
push!(cobweb_ys, return_map_ys[end])

# Plot the return map.
begin
  fig = Figure()
  ax = Axis(fig[1, 1],
    aspect=1,
    title="Return map for Thomas attractor",
    titlesize=30,
    xlabel=L"\hat{x}_i",
    ylabel=L"\hat{x}_{i+1}",
    xlabelsize=30,
    ylabelsize=30
  )
  lines!(ax, cobweb_xs, cobweb_ys, color=:black, linewidth=0.03)
  lines!(ax, return_map_endpoints, return_map_endpoints, color=:gray)
  scatter!(ax, return_map_xs, return_map_ys, color=:blue, markersize=5)
  display(fig)

  # Optionally save the final figure.
  # save("thomas_return_map.png", fig)
end

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

  # Plot points of section.
  scatter!(ax, [
    Point3f(point...) for point in sphere_exits[section_iterate_offset:section_iterate:end]
  ], color=:red, markersize=10)

  # Plot section-pruning sphere.
  u = range(0, 2π, length=100)
  v = range(0, π, length=100)
  
  x = [radius * cos(u) * sin(v) for u in u, v in v]
  y = [radius * sin(u) * sin(v) for u in u, v in v]
  z = [radius * cos(v) for u in u, v in v]

  surface!(ax, x, y, z, color=:blue, alpha=0.5, transparency=true)

  # Set equal aspect ratio
  ax.aspect = :data

  # Display the figure.
  display(fig)

  # Optionally save the final figure.
  # save("thomas_attractor.png", fig)
end
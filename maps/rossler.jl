using Pkg
Pkg.activate(".")

using Colors
using DifferentialEquations
using GLMakie

function rossler_derivatives!(du, u, p, t)
    x, y, z = u
    a, b, c = p
    du[1] = -y - z
    du[2] = x + a * y
    du[3] = b + z * (x - c)
end

u0 = [0.1, 0.0, 0.0]
p = [0.38, 0.3, 5.5]

prob = ODEProblem(rossler_derivatives!, u0, (0.0, 1e5), p)
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

# Capture all x minima in the solution where z is small.
# Find indices where x has a local minimum and z is small.
x_minima = Int[]
for i in 2:(length(sol.u)-1)
    x_prev = sol.u[i-1][1]
    x_curr = sol.u[i][1] 
    x_next = sol.u[i+1][1]
    z_curr = sol.u[i][3]
    
    # Check if this is a local minimum and z is small enough.
    if x_curr < x_prev && x_curr < x_next && z_curr < 0.1
        push!(x_minima, i)
    end
end

begin # Plot 3D state space with return section points.
  fig = Figure()
  ax = Axis3(fig[1, 1], xlabel="x", ylabel="y", zlabel="z")

  xs = [u[1] for u in sol.u]
  ys = [u[2] for u in sol.u]
  zs = [u[3] for u in sol.u]

  lines!(ax, xs, ys, zs, color=:blue)
  scatter!(ax, xs[x_minima], ys[x_minima], zs[x_minima], color=:red)

  display(fig)
end

begin # Plot return map of x values for x_minima.
  fig = Figure()
  ax = Axis(fig[1, 1], xlabel="x", ylabel="x", aspect=1)
  
  # Get min and max x values to set limits for identity line
  xmin = minimum(xs[x_minima])
  xmax = maximum(xs[x_minima])
  
  # Plot identity line
  lines!(ax, [xmin, xmax], [xmin, xmax], color=:black, linestyle=:dash)
  
  # Plot return map points
  scatter!(ax, xs[x_minima][1:end-1], xs[x_minima][2:end], color=:red)
  display(fig)
end

# Calculate return map from ordinal partition network.
include("../return_maps/return_map.jl")
return_map_itinerary = calculate_return_itinerary(
  rossler_derivatives!,
  p,
  u0,
  200.0,
  0.01,
  (0.0, 1e3),
  (0.0, 1e5),
  u -> u[1],
  u -> u[1],
  u -> u[3] < 0.1,
  3,
  4,
  2,
  1,
  3,
  1
)

iterate = 1
dbscan_ε = 3e-2
y_scaling = 3e-2
filter_ε = 1e-1
critical_points, boundary_indices, orientations = kneading_increment_idxs(
  return_map_itinerary[1:iterate:end-iterate],
  dbscan_ε,
  y_scaling,
  filter_ε
)

include("../kneading/sequences.jl")
K = 6 # Number of iterates to use in the truncated kneading matrix calculation.
kneading_sequences = [
  Int8.(seq) for seq in
  calculate_kneading_sequences(
    boundary_indices,
    critical_points,
    return_map_itinerary,
    K
  )
]
include("../kneading/matrix.jl")
kneading_matrix = allocate_kneading_matrix(critical_points, K)
kneading_matrix_from_kneading_sequences!(kneading_matrix, kneading_sequences, orientations, K)

include("../kneading/encodings.jl")
encoded_kneading_matrix = matrix_encoding(kneading_matrix)
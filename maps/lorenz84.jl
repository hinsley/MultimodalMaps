using Pkg
Pkg.activate(".")

using Colors

include("../return_maps/return_map.jl")

function lorenz84_derivatives!(du, u, p, t)
  x, y, z = u
  a, b, F, G = p
  du[1] = -y^2-z^2+a*(F-x)
  du[2] = x*y-b*x*z-y+G
  du[3] = b*x*y+z*(x-1.0)
end

u0 = [0.5, 0.0, 0.0]
p = [0.25, 4.0, 8.0, 1.0]
transient_time = 3e0
dt = 1e-3
tspan1 = (0.0, 1e3)
tspan2 = (1e3, 3e4)
ordinal_function = u -> u[3]
return_map_coordinate = ordinal_function
section_constraint = u -> ordinal_function(u) < 1
iterate = 1 # Compositional power of the return map.
w1 = 1
m1 = 4
τ1 = 5
w2 = 1
m2 = 3
τ2 = 1
ranked_ordinal_symbol_index = 1

return_map_itinerary = calculate_return_itinerary(
  lorenz84_derivatives!,
  p,
  u0,
  transient_time,
  dt,
  tspan1,
  tspan2,
  ordinal_function,
  return_map_coordinate,
  section_constraint,
  w1,
  w2,
  m1,
  m2,
  τ1,
  τ2,
  ranked_ordinal_symbol_index
)

# include("../return_maps/return_map.jl")
# dbscan_ε = 3e-5
# y_scaling = 3e-2
# filter_ε = 2.5e-5
# clusters, simplified_clusters = kneading_increment_idxs(
#   return_map_itinerary[1:iterate:end-iterate],
#   dbscan_ε,
#   y_scaling,
#   filter_ε,
# )

iterate = 1 # Compositional power of the return map.

# Construct the return map graph.
return_map_xs = return_map_itinerary[1:iterate:end-iterate]
return_map_ys = return_map_itinerary[1+iterate:iterate:end]
return_map_endpoints = [
  minimum(return_map_xs),
  maximum(return_map_xs)
]
cobweb_xs = [return_map_xs[1]]
cobweb_ys = [return_map_ys[1]]
for i in 1:length(return_map_xs)-1
    push!(cobweb_xs, cobweb_xs[end], return_map_xs[i+1])
    push!(cobweb_ys, return_map_ys[i], return_map_ys[i])
end
push!(cobweb_xs, cobweb_xs[end])
push!(cobweb_ys, return_map_ys[end])

# Plot the return map.
begin
  fig = Figure()
  ax = Axis(
    fig[1, 1],
    aspect=1,
    title="Return map",
    titlesize=30,
    xlabel=L"x_i",
    ylabel=L"x_{i+1}",
    xlabelsize=30,
    ylabelsize=30
  )
  lines!( # Plot cobweb.
    ax,
    cobweb_xs,
    cobweb_ys,
    color=:black,
    linewidth=0.01
  )
  lines!( # Plot identity line.
    ax,
    return_map_endpoints,
    return_map_endpoints,
    color=:gray
  )
  scatter!( # Plot return map graph.
    ax,
    return_map_xs,
    return_map_ys,
    color=:blue,
    markersize=2
  )
  # cluster_colors = []
  # for cluster in clusters
  #   push!(cluster_colors, rand(RGB))
  #   scatter!(
  #     ax,
  #     cluster,
  #     markersize=5,
  #     color=(cluster_colors[end], 0.2)
  #   )
  # end
  # # Plot decimated return map pieces.
  # for (simplified_cluster, cluster_color) in zip(simplified_clusters, cluster_colors)
  #   lines!(
  #     ax,
  #     simplified_cluster,
  #     linewidth=3,
  #     color=cluster_color
  #   )
  # end
  # Plot circles of radius dbscan_ε around each point
  # for (x, y) in zip(return_map_xs, return_map_ys)
  #   points = [Point2f(x + dbscan_ε * cos(t), y + dbscan_ε * sin(t) / y_scaling) for t in range(0, 2π, length=100)]
  #   lines!(ax, points, color=(:black, 0.2))
  # end
  
  # Plot extrema.
  # for cluster_extrema in extrema
  #   for (x, y) in cluster_extrema
  #     scatter!(
  #       ax,
  #       [x],
  #       [y],
  #       color=:black,
  #       markersize=10
  #     )
  #   end
  # end

  # Plot sampled point curves.
  # for segment in points
  #   lines!(
  #     ax,
  #     segment[1],
  #     segment[2],
  #     color=:red,
  #     linewidth=2
  #   )
  #   println("Plotting segment from x=$(segment[1][1]) to x=$(segment[1][end])")
  # end

  display(fig)
end

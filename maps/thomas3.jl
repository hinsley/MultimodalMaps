using Pkg
Pkg.activate(".")

using Colors

include("../return_maps/return_map.jl")

function thomas_derivatives!(du, u, p, t)
  x, y, z = u
  b = p[1]
  du[1] = sin(y) - b * x
  du[2] = sin(z) - b * y
  du[3] = sin(x) - b * z
end

u0 = [1.0, 1.0, 0.5]
p = [0.208186]
transient_time = 2e2
dt = 3e-3
tspan1 = (0.0, 1e3)
tspan2 = (1e3, 7e5)
ordinal_function = u -> sqrt(u[1]^2+u[2]^2+u[3]^2)
return_map_coordinate = ordinal_function
section_constraint = u -> ordinal_function(u) > 4.41
iterate = 4 # Compositional power of the return map.
w1 = 1
m1 = 4
τ1 = 5
w2 = 1
m2 = 3
τ2 = 1
ranked_ordinal_symbol_index = 2

return_map_itinerary = calculate_return_itinerary(
  thomas_derivatives!,
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

include("../return_maps/return_map.jl")
dbscan_ε = 3e-5
y_scaling = 3e-2
filter_ε = 2.5e-5
clusters, simplified_clusters = laps(
  return_map_itinerary[1:iterate:end-iterate],
  dbscan_ε,
  y_scaling,
  filter_ε,
)

# return_map, critical_points = interpolate_return_map(
# include("../return_maps/return_map.jl")
# return_map, points = interpolate_return_map(
#   return_map_itinerary[1:iterate:end-iterate],
#   4,
#   200,
#   2,
#   # true
# )

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
    linewidth=0.03
  )
  lines!( # Plot identity line.
    ax,
    return_map_endpoints,
    return_map_endpoints,
    color=:gray
  )
  # scatter!( # Plot return map graph.
  #   ax,
  #   return_map_xs,
  #   return_map_ys,
  #   color=:blue,
  #   markersize=5
  # )
  cluster_colors = []
  for cluster in clusters
    push!(cluster_colors, rand(RGB))
    scatter!(
      ax,
      cluster,
      markersize=5,
      color=(cluster_colors[end], 0.2)
    )
  end
  # Plot decimated return map pieces.
  for (simplified_cluster, cluster_color) in zip(simplified_clusters, cluster_colors)
    lines!(
      ax,
      simplified_cluster,
      linewidth=3,
      color=cluster_color
    )
  end
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
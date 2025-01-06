using Pkg
Pkg.activate(".")

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
dt = 1e-2
tspan1 = (0.0, 1e3)
tspan2 = (1e3, 1e5)
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
  scatter!( # Plot return map graph.
    ax,
    return_map_xs,
    return_map_ys,
    color=:blue,
    markersize=5
  )

  display(fig)
end
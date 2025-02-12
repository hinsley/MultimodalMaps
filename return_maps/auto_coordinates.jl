using LinearAlgebra

function auto_coordinates(
  points::Vector{Vector{Float64}}, # Point cloud of connected return map section.
  decimation_threshold::Float64, # Douglas-Peucker ε.
  distance_function, # Function which takes a vector of length twice the dimension of the space, and returns a scalar.
  truncation_quantity::Int64 = 0; # Number of longest segments to discard.
  return_curve::Bool = false # Whether to return the curve itself instead of a function.
)::Union{Function, Vector{Vector{Float64}}}
  # Create a doubly linked list representing a piecewise linear curve passing once through each point.
  curve_points = Vector{Float64}[points[1]]
  # Keep track of points not yet added to the curve.
  remaining_points = Vector{Float64}[points[2:end]...]

  # Find the nearest among the remaining points to the specified point.
  function find_nearest_point(
    point::Vector{Float64},
    remaining_points::Vector{Vector{Float64}}
  )::Tuple{Int, Float64}
    min_distance = Inf
    nearest_point_idx = nothing
    for (index, candidate) in enumerate(remaining_points)
      distance = distance_function(vcat(point, candidate))
      if distance < min_distance
        min_distance = distance
        nearest_point_idx = index
      end
    end
    return nearest_point_idx, min_distance
  end

  # Keep track of index of first-added point as the curve is reconstructed.
  first_point_idx = 1

  # Iterate over remaining points until all are added to the curve.
  while length(remaining_points) > 0
    # Find the nearest point to the start and end of curve_points.
    nearest_idx_start, nearest_dist_start = find_nearest_point(curve_points[1], remaining_points)
    nearest_idx_end, nearest_dist_end = find_nearest_point(curve_points[end], remaining_points)
    # Check if the nearest point to the start is closer than the nearest point to the end;
    # extend the curve at the appropriate endpoint by the nearer of the two.
    if nearest_dist_start < nearest_dist_end
      # Update the index of the first-added point.
      first_point_idx += 1
      # Add the nearest point to the start of the curve.
      pushfirst!(curve_points, remaining_points[nearest_idx_start])
      # Remove the nearest point from the remaining points.
      splice!(remaining_points, nearest_idx_start)
    else
      # Add the nearest point to the end of the curve.
      push!(curve_points, remaining_points[nearest_idx_end])
      # Remove the nearest point from the remaining points.
      splice!(remaining_points, nearest_idx_end)
    end
  end

  for _ in 1:truncation_quantity
    # Figure out which line segment is the longest.
    max_distance = 0.0
    chop_idx = 0
    for i in 1:length(curve_points)-1
      distance = distance_function(vcat(curve_points[i], curve_points[i+1]))
      if distance > max_distance
        max_distance = distance
        chop_idx = i
      end
    end
    # Chop the curve at the largest segment.
    if chop_idx < first_point_idx
      curve_points = curve_points[chop_idx+1:end]
      first_point_idx -= chop_idx
    else
      curve_points = curve_points[1:chop_idx]
    end
  end

  # Decimate the constructed curve.
  decimated_curve = decimate(curve_points, decimation_threshold, distance_function)

  # If we don't care about returning a coordinate function, we're done.
  if return_curve
    return decimated_curve
  end

  # Build projection function using decimated curve.
  function project_to_curve(q::Vector{Float64})::Float64
    min_dist = Inf
    nearest_t = 0.0
    total_length = 0.0
    
    # Precompute segment lengths and total path length.
    segment_lengths = [
      distance_function(vcat(decimated_curve[i+1], decimated_curve[i]))
      # 1.0 # Use same length for each segment.
      for i in 1:length(decimated_curve)-1
    ]
    total_length = sum(segment_lengths)
    
    # Find nearest segment.
    for i in 1:length(decimated_curve)-1
      a, b = decimated_curve[i], decimated_curve[i+1]
      v = b - a
      t = clamp(dot(q - a, v) / dot(v, v), 0.0, 1.0)
      proj = a + t .* v
      dist = distance_function(vcat(q, proj))
      
      if dist < min_dist
        min_dist = dist
        # Calculate accumulated length up to this segment.
        accumulated_length = sum(segment_lengths[1:i-1])
        current_segment_progress = t * segment_lengths[i]
        nearest_t = (accumulated_length + current_segment_progress) / total_length
      end
    end
    return nearest_t
  end

  return project_to_curve
end
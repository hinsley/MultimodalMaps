function decimate_graph(
  xs::Vector{Float64},
  ys::Vector{Float64},
  ε::Float64,
  exact_extrema::Bool = false
  )::Tuple{Vector{Float64}, Vector{Float64}}
  @assert length(xs) == length(ys) "xs and ys must have the same length"
  
  # A helper function to compute the perpendicular distance of the point (x, y).
  # to the line defined by (x1, y1) and (x2, y2).
  distance_to_segment(x, y, x1, y1, x2, y2) = begin
    # If the two end points are the same, return the distance to that point.
    dx = x2 - x1
    dy = y2 - y1
    if dx == 0.0 && dy == 0.0
      return sqrt((x - x1)^2 + (y - y1)^2)
    end
    
    # Compute the projection of the point onto the line parameter t.
    t = ((x - x1)*dx + (y - y1)*dy) / (dx*dx + dy*dy)
    # Constrain t to [0, 1] so we stay within the segment.
    t_clamped = max(min(t, 1.0), 0.0)
    
    # Find projected point on the line.
    proj_x = x1 + t_clamped*dx
    proj_y = y1 + t_clamped*dy
    # Distance from the point to the projected point.
    return sqrt((x - proj_x)^2 + (y - proj_y)^2)
  end
  
  # Recursive function that performs the Douglas-Peucker logic.
  function douglas_peucker(xcoords, ycoords, start::Int, finish::Int, eps)
    # If there are only two points, we cannot simplify further
    if finish - start < 2
      return [start, finish]
    end
    
    # 1. Find the point of maximum distance from the segment [start, finish].
    x1, y1 = xcoords[start], ycoords[start]
    x2, y2 = xcoords[finish], ycoords[finish]
    
    max_dist = -1.0
    max_idx = start
    for i in start+1:finish-1
      dist = distance_to_segment(xcoords[i], ycoords[i], x1, y1, x2, y2)
      if dist > max_dist
        max_dist = dist
        max_idx = i
      end
    end
    
    # 2. If the max distance is above our tolerance, we split and recurse.
    if max_dist > eps
      left_indices  = douglas_peucker(xcoords, ycoords, start, max_idx, eps)
      right_indices = douglas_peucker(xcoords, ycoords, max_idx, finish, eps)
      # Combine them, removing one duplicate to avoid repeating the middle.
      return vcat(left_indices[1:end-1], right_indices)
    else
      # Otherwise keep only the endpoints of this segment.
      return [start, finish]
    end
  end
  
  # Run the main recursion on the entire array.
  idxs = douglas_peucker(xs, ys, 1, length(xs), ε)
  
  # Build the new arrays from the indices.
  new_xs = xs[idxs]
  new_ys = ys[idxs]

  if exact_extrema
    refined_new_xs = new_xs
    refined_new_ys = new_ys
    # Find the extrema of the decimated array.
    # Note: This can introduce artifacts if two extrema
    # are adjacent in the decimated array.
    for i in length(new_xs)-1:-1:2
      # If point is a local extremum, refine.
      if (
        new_ys[i] > new_ys[i-1]
        && new_ys[i] > new_ys[i+1]
      )
        # Local maximum.
        candidate_idxs = idxs[i-1]+1:idxs[i+1]-1
        max_idx = argmax(ys[candidate_idxs]) + candidate_idxs[1] - 1
        # Insert refined maximum at correct position.
        if xs[max_idx] < new_xs[i]
          insert!(refined_new_xs, i, xs[max_idx])
          insert!(refined_new_ys, i, ys[max_idx])
        elseif xs[max_idx] > new_xs[i]
          insert!(refined_new_xs, i+1, xs[max_idx]) 
          insert!(refined_new_ys, i+1, ys[max_idx])
        end
      elseif (
        new_ys[i] < new_ys[i-1]
        && new_ys[i] < new_ys[i+1]
      )
        # Local minimum.
        candidate_idxs = idxs[i-1]+1:idxs[i+1]-1
        min_idx = argmin(ys[candidate_idxs]) + candidate_idxs[1] - 1
        # Insert refined minimum at correct position.
        if xs[min_idx] < new_xs[i]
          insert!(refined_new_xs, i, xs[min_idx])
          insert!(refined_new_ys, i, ys[min_idx])
        elseif xs[min_idx] > new_xs[i]
          insert!(refined_new_xs, i+1, xs[min_idx])
          insert!(refined_new_ys, i+1, ys[min_idx])
        end
      end
    end
    new_xs = refined_new_xs
    new_ys = refined_new_ys
  end
  
  return (new_xs, new_ys)
end

function decimate(
  points::Vector{Vector{Float64}},
  ε::Float64,
  distance_function
)::Vector{Vector{Float64}}
  # Helper to compute distance from point to line segment in R^n
  function distance_to_segment(p::Vector{Float64}, a::Vector{Float64}, b::Vector{Float64})
    v = b .- a
    if all(v .== 0)
      return distance_function(vcat(p, a))
    end
    
    t = clamp(dot(p .- a, v) / dot(v, v), 0.0, 1.0)
    projection = a .+ t .* v
    return distance_function(vcat(p, projection))
  end

  # Recursive DP implementation for R^n
  function douglas_peucker_rn(pts::Vector{Vector{Float64}}, first::Int, last::Int, eps::Float64)
    if last - first < 2
      return [first, last]
    end

    a, b = pts[first], pts[last]
    max_dist = -1.0
    max_idx = first

    for i in first+1:last-1
      dist = distance_to_segment(pts[i], a, b)
      if dist > max_dist
        max_dist = dist
        max_idx = i
      end
    end

    if max_dist > eps
      left = douglas_peucker_rn(pts, first, max_idx, eps)
      right = douglas_peucker_rn(pts, max_idx, last, eps)
      return vcat(left[1:end-1], right)
    else
      return [first, last]
    end
  end

  idxs = douglas_peucker_rn(points, 1, length(points), ε)
  return points[idxs]
end
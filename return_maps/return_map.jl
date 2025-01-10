using Pkg
Pkg.activate(".")

using Clustering
using DataStructures
using DifferentialEquations
using GLMakie
using Peaks
using Statistics

include("../return_maps/opn.jl")
include("../return_maps/douglas_peucker.jl")

# Constructs a sequence of function values 
function calculate_return_itinerary(
  derivatives!::Function,
  p::Vector{Float64},
  u0::Vector{Float64},
  transient_time::Float64,
  dt::Float64,
  tspan1::Tuple{Float64, Float64},
  tspan2::Tuple{Float64, Float64}, # Should start where tspan1 ends.
  ordinal_function::Function, # Should map state space to reals.
  return_map_coordinate::Function, # Should map state space to reals.
  section_constraint::Function, # Takes state vector as argument, returns Bool. If false, discards points from Poincare section.
  w1::Int, # OPN parameters for original timeseries of each trajectory.
  w2::Int,
  m1::Int,
  m2::Int, # OPN parameters for each ordinal symbol's constructed timeseries.
  τ1::Int,
  τ2::Int,
  ranked_ordinal_symbol_index::Int = 1 # If 1, use highest weighted entropy ordinal pattern. If higher, use lower entropy. Maximum of factorial(m1).
)::Vector{Float64}
  return_map_itinerary = Float64[]

  prob1 = ODEProblem(derivatives!, u0, tspan1, p)

  println("Integrating trajectory 1.")
  @time begin
    sol1_full = solve(
      prob1,
      Tsit5(),
      reltol=1e-8,
      abstol=1e-8,
      dt=dt,
      adaptive=false,
      maxiters=1e10,
      save_everystep=true
    )
    # Get index of first point after transient_time.
    start_idx = findfirst(t -> t >= transient_time, sol1_full.t)
    # Create new solution with transient discarded.
    sol1 = DiffEqBase.build_solution(
      prob1,
      sol1_full.alg,
      sol1_full.t[start_idx:end],
      sol1_full.u[start_idx:end]
    )
  end

  # Get weighted entropy for each ordinal symbol.
  xs = [ordinal_function(u) for u in sol1.u]
  ordinal_symbols1 = encode_chronological_timeseries(xs, w1, m1, τ1)
  weighted_entropies = Float64[]
  println("Computing weighted entropies for trajectory 1.")
  @time for ordinal_symbol in 1:factorial(m1)
    # Construct the ordinal pattern's timeseries.
    ordinal_timeseries1 = get_ordinal_timeseries(
      xs,
      ordinal_symbols1,
      ordinal_symbol,
      w1
    )
    # If the ordinal pattern never occurs, the weighted entropy is zero.
    if length(ordinal_timeseries1) == 0
      push!(weighted_entropies, 0)
    else
      # Construct the OPN for the ordinal pattern's timeseries.
      ordinal_symbols2 = encode_chronological_timeseries(
        ordinal_timeseries1,
        w2,
        m2,
        τ2
      )
      # Calculate and store this ordinal pattern's weighted entropy.
      push!(
        weighted_entropies,
        weighted_entropy(
          length(xs),
          w1,
          m1,
          τ1,
          ordinal_symbols2,
          m2
        )
      )
    end
  end
  println("Weighted entropies (sorted): $(reverse(sort(weighted_entropies)))")

  # Select an ordinal pattern to use for constructing the return map.
  # Filter out NaN values and sort remaining indices
  valid_indices = findall(!isnan, weighted_entropies)
  sorted_valid_indices = sort(valid_indices, by=i->weighted_entropies[i], rev=true)
  section_ordinal_symbol = sorted_valid_indices[ranked_ordinal_symbol_index]

  ### Don't be wasteful! Use this first trajectory to construct some of the return map.
  
  # Get indices of points in the Poincare section.
  section_idxs = Int[]
  for i in 2:length(ordinal_symbols1)
    sol1_idx = (i-1)*w1+1
    if (
      ordinal_symbols1[i] == section_ordinal_symbol
      && ordinal_symbols1[i-1] != section_ordinal_symbol
      && section_constraint(sol1.u[sol1_idx])
    )
      push!(section_idxs, sol1_idx)
    end
  end

  # Add the points to the return map.
  append!(
    return_map_itinerary,
    [return_map_coordinate(sol1.u[i]) for i in section_idxs]
  )

  ### Construct the return map during integration.
  prob2 = ODEProblem(derivatives!, sol1.u[end], tspan2, p)

  # Set up integrator callback for capturing Poincare section hits.
  ordinal_indices = get_ordinal_indices(m1) # Precompute a lookup table for ordinal symbols.
  prev_ordinal_symbol = ordinal_symbols1[end] # Previous ordinal symbol for detecting changes in ordinal pattern.
  l1 = (m1-1)*τ1+1 # Capture window size.
  capture_window_delay_counter = 0
  capture_window_buffer = CircularBuffer{Float64}(l1)
  u_buffer = CircularBuffer{Vector{Float64}}(l1) # State vector buffer.
  # Populate the buffers with values from the first trajectory.
  append!(capture_window_buffer, xs[end-l1+1:end])
  append!(u_buffer, sol1.u[end-l1+1:end])

  function condition(u, t, integrator)
    # Increment counters.
    capture_window_delay_counter += 1
    # Capture ordinal value for current state.
    push!(capture_window_buffer, ordinal_function(u))
    # Capture current state.
    push!(u_buffer, collect(u))
    # Check if we are ending a capture window.
    if capture_window_delay_counter == w1
      # Reset the delay counter.
      capture_window_delay_counter = 0
      # Determine the current ordinal symbol.
      ordinal_symbol = encode_chronological(
        collect(capture_window_buffer)[1:τ1:end],
        ordinal_indices
      )
      # Check if we are in a valid region of state space and
      # the ordinal symbol has changed to the ordinal symbol of interest.
      if (
        ordinal_symbol == section_ordinal_symbol
        && prev_ordinal_symbol != section_ordinal_symbol
        && section_constraint(u)
      )
        # Add the first point in the window to the return map.
        push!(return_map_itinerary, return_map_coordinate(u_buffer[1]))
      end
      # Save the current ordinal symbol for next step.
      prev_ordinal_symbol = ordinal_symbol
    end

    return false
  end
  cb = DiscreteCallback(condition, nothing)

  println("Integrating trajectory 2.")
  @time sol2 = solve(
    prob2,
    Tsit5(),
    reltol=1e-8,
    abstol=1e-8,
    dt=dt,
    adaptive=false,
    callback=cb,
    save_everystep=false,
    maxiters=1e10
  )

  return return_map_itinerary
end

# Finds the lap endpoints of a return map by clustering the points
# in the graph reconstructed from the itinerary.
function laps(
  itinerary::Vector{Float64},
  dbscan_ε::Float64,
  y_scaling::Float64,
  filter_ε::Float64
)
# )::Vector{Vector{Tuple{Float64, Float64}}}
  xs = itinerary[1:end-1]
  ys = itinerary[2:end]
  order = sortperm(xs)
  point_matrix = hcat(xs[order], ys[order] * y_scaling)'
  @time dbscan_result = dbscan(point_matrix, dbscan_ε)
  cluster_endpoints = Tuple{Float64, Float64}[]
  clusters = Vector{Tuple{Float64, Float64}}[]
  simplified_clusters = Vector{Tuple{Float64, Float64}}[]
  extrema = Vector{Tuple{Float64, Float64}}[]

  for cluster in dbscan_result.clusters
    cluster_xs = xs[order][cluster.core_indices]
    cluster_ys = ys[order][cluster.core_indices]
    cluster_points = collect(zip(cluster_xs, cluster_ys))
    push!(clusters, cluster_points)
    push!(cluster_endpoints, (cluster_xs[1], cluster_xs[end]))
    println("Decimating cluster.")
    @time simplified_cluster_xs, simplified_cluster_ys = decimate_function(
      cluster_xs,
      cluster_ys,
      filter_ε,
      true
    )
    simplified_cluster_points = collect(zip(simplified_cluster_xs, simplified_cluster_ys))
    push!(simplified_clusters, simplified_cluster_points)
  end


  return clusters, simplified_clusters
end

# Interpolate a return map function.
function interpolate_return_map(
  itinerary::Vector{Float64},
  pieces::Int,
  piece_resolution::Int,
  sparsity::Int=1,
  return_critical_points::Bool=false
)::Union{Tuple{Function, Vector{Tuple{Vector{Float64}, Vector{Float64}}}}, Tuple{Function, Vector{Float64}}}
  xs = itinerary[1:end-1]
  ys = itinerary[2:end]
  order = sortperm(xs)[1:sparsity:end]
  # Calculate differences between consecutive y values after sorting by x.
  y_diffs = diff(ys[order])
  # Discontinuities are where the graph suddenly changes drastically
  # in y.
  discontinuities = sort(sortperm(abs.(y_diffs), rev=true)[1:pieces-1])
  # Segment the domain into points lying on each continuous piece.
  segments = Vector{Tuple{Int, Int}}()
  start_idx = 1
  for d in discontinuities
    push!(segments, (start_idx, d))
    start_idx = d + 1
  end
  push!(segments, (start_idx, length(xs))) # Last segment.
  println([xs[order][seg_end] for (seg_start, seg_end) in segments[1:end-1]])
  println(segments[1][2])
  # For each continuous domain segment, sample `piece_resolution` points.
  piece_data = Vector{Tuple{Vector{Float64}, Vector{Float64}}}(undef, length(segments))
  for (seg_i, (seg_start, seg_end)) in enumerate(segments)
    point_idxs = [round(Int, x) for x in range(seg_start, seg_end, length=piece_resolution)]
    point_idxs = [argmin(abs.(x .- (1:sparsity:length(xs)))) for x in point_idxs]
    piece_data[seg_i] = (xs[order][point_idxs], ys[order][point_idxs])
  end
  # Construct the piecewise continuous function which is linearly
  # interpolated between sampled points.
  function map(x::Float64)::Float64
    # Determine to which segment x belongs.
    for (seg_i, (seg_start, seg_end)) in enumerate(segments)
      if x >= xs[order][seg_start] && x <= xs[order][seg_end]
        # Find the two sampled points to interpolate between.
        for i in 1:piece_resolution-1
          left_x = piece_data[seg_i][1][i]
          right_x = piece_data[seg_i][1][i+1]
          if left_x <= x <= right_x
            # Linearly interpolate between the sampled points.
            left_y = piece_data[seg_i][2][i]
            right_y = piece_data[seg_i][2][i+1]
            return left_y + (right_y - left_y) * (x - left_x) / (right_x - left_x)
          end
        end
      end
    end
    # If x is outside the domain, throw a DomainError.
    throw(DomainError(x, "input outside of domain"))
  end

  if return_critical_points
    # Leftmost and rightmost points of domain are critical points.
    critical_points = [xs[order][1], xs[order][end]]
    # Push the discontinuities between segments.
    for seg_i in 1:length(segments)-1
      push!(
        critical_points,
        [
          # Right side of left segment.
          xs[order][segments[seg_i][2]],
          # Left side of right segment.
          xs[order][segments[seg_i+1][1]]
        ]
      )
    end
    # Push the minima and maxima in the interior of each segment.
    for segment_points in piece_data
      # Find local extrema by checking points against neighbors.
      for i in 2:length(segment_points[2])-1
          y_prev = segment_points[2][i-1]
          y_curr = segment_points[2][i]
          y_next = segment_points[2][i+1]
          
          # Local maximum
          if y_curr > y_prev && y_curr > y_next
              push!(critical_points, [segment_points[1][i]])
          end
          # Local minimum 
          if y_curr < y_prev && y_curr < y_next
              push!(critical_points, [segment_points[1][i]])
          end
      end
    end

    return map, sort(critical_points)
  else
    return map, piece_data
  end
end
using Pkg
Pkg.activate(".")

using DataStructures
using DifferentialEquations
using GLMakie

include("../return_maps/opn.jl")

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
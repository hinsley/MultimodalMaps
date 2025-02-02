function calculate_kneading_sequences(
  boundary_indices::Vector{Tuple{Int, Int}},
  critical_points::Vector{Float64},
  return_map_itinerary::Vector{Float64},
  K::Int
)::Vector{Vector{Int}}
  kneading_sequences = Vector{Int}[]
  for (i, (left_idx, right_idx)) in enumerate(boundary_indices)
    left_kneading_sequence = [i]
    right_kneading_sequence = [i+1]
    for k in 1:K # Iterate over iterates of the map.
      left_iterate_lap = findlast(
        critical_point -> critical_point < return_map_itinerary[left_idx + k],
        critical_points
      )
      if isnothing(left_iterate_lap) # Iterate lies within first lap.
        left_iterate_lap = 1
      else
        left_iterate_lap += 1
      end
      right_iterate_lap = findlast(
        critical_point -> critical_point < return_map_itinerary[right_idx + k],
        critical_points
      )
      if isnothing(right_iterate_lap) # Iterate lies within first lap.
        right_iterate_lap = 1
      else
        right_iterate_lap += 1
      end
      push!(left_kneading_sequence, left_iterate_lap)
      push!(right_kneading_sequence, right_iterate_lap)
    end
    append!(kneading_sequences, [left_kneading_sequence, right_kneading_sequence])
  end
  return kneading_sequences
end
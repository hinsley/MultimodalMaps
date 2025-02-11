# Calculate the truncated Milnor-Thurston kneading matrix.
function kneading_matrix!(matrix::Array{Int8}, map, critical_points, p)
  # TODO: Allow for discontinuous maps.
  orientation = map(p, critical_points[1]) > map(p, critical_points[2]) ? 1 : -1
  m = length(critical_points) # Number of critical points.
  # n = size(matrix, 2)         # Number of laps.
  K = size(matrix, 3) - 1     # Number of forward iterates.
  for i in 1:m
    matrix[i, i,   1] =  1 * orientation * (isodd(i) ? 1 : -1)
    matrix[i, i+1, 1] = -1 * orientation * (isodd(i) ? 1 : -1)
    point = critical_points[i]
    cumulative_orientation = orientation
    for k in 2:K+1
      # Iterate point forward.
      point = map(p, point)
      # Find the lap index of the iterate.
      lap = 1
      for j in 1:m
        if point >= critical_points[j]
          lap = j + 1
        else
          break
        end
      end
      # Add coefficient to kneading matrix.
      cumulative_orientation *= orientation * (isodd(lap) ? 1 : -1)
      matrix[i, lap, k] += 2 * cumulative_orientation
    end
  end
end

function kneading_matrix_from_kneading_sequences!(
  matrix::Array{Int8},
  kneading_sequences::Vector{Vector{Int8}},
  orientations::Vector{Bool},
  K::Int
)
  m = length(orientations) - 1 # Number of critical points.
  @assert length(kneading_sequences) == 2 * m
  
  for i in 1:m
    # Process left kneading sequence (odd index).
    left_seq = kneading_sequences[2*i-1]
    cumulative_orientation = orientations[i] ? 1 : -1
    lap = i # Start in current lap.
    
    # Initialize zeroth iterate (k=0).
    matrix[i, lap, 1] += cumulative_orientation
    
    for k in 2:K
      lap = abs(left_seq[k])
      cumulative_orientation *= orientations[lap] ? 1 : -1
      matrix[i, lap, k] -= cumulative_orientation
    end
    
    # Process right kneading sequence (even index).
    right_seq = kneading_sequences[2*i]
    cumulative_orientation = orientations[i+1] ? 1 : -1
    lap = i + 1 # Start in next lap.
    
    # Initialize zeroth iterate (k=0).
    matrix[i, lap, 1] += cumulative_orientation
    
    for k in 2:K+1
      lap = abs(right_seq[k])
      cumulative_orientation *= orientations[lap] ? 1 : -1
      matrix[i, lap, k] += cumulative_orientation
    end
  end
end

function allocate_kneading_matrix(critical_points, K)
  m = length(critical_points) # Number of critical points.
  n = m+1 # Number of laps.
  # K is the degree of the truncated power series entries
  # of the matrix.
  return zeros(Int8, m, n, K+1)
end
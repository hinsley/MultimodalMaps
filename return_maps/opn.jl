# Ordinal partition networks
# doi: 10.1063/5.0141438

# The parameters for the ordinal partition network are:
# - w: The number of points in a window before beginning the
#   successive window ("lag").
# - m: The number of points to sample per window ("embedding dimension").
# - τ: The sampling delay between points within a window. 1 is no delay.

using Pkg
Pkg.activate(".")

using Combinatorics

# Get a dictionary mapping each ordinal partition to its index.
function get_ordinal_indices(m::Int)::Dict{Vector{Int}, Int}
    # Generate all possible permutations of indices 1:m.
    perms = collect(permutations(1:m))
    
    # Create a dictionary mapping each permutation to its index.
    lookup = Dict{Vector{Int}, Int}()
    for (i, p) in enumerate(perms)
        lookup[p] = i
    end
    
    return lookup
end

# Encode the samples from a time series window into an ordinal
# symbol using the chronological index ranking method.
function encode_chronological(
  samples::Vector{Float64},
  ordinal_indices::Dict{Vector{Int}, Int}
)::Int
  sorted_indices = collect(sortperm(samples))
  return ordinal_indices[sorted_indices]
end

# Encode a time series into a sequence of ordinal symbols using
# the chronological index ranking method.
function encode_chronological_timeseries(
  xs::Vector{Float64},
  w::Int,
  m::Int,
  τ::Int
)::Vector{Int}
  ordinal_indices = get_ordinal_indices(m)
  l = (m-1)*τ
  return [
    encode_chronological(
      xs[starting_index:τ:starting_index+l],
      ordinal_indices
    )
    for starting_index in 1:w:length(xs)-l
  ]
end

# Construct a time series for a particular ordinal symbol.
function get_ordinal_timeseries(
  xs::Vector{Float64},
  ordinal_symbols::Vector{Int},
  symbol::Int,
  w::Int, # Window size for original timeseries.
)::Vector{Float64}
  return [
    xs[(i-1)*w + 1]
    for i in 1:length(ordinal_symbols)
    if ordinal_symbols[i] == symbol
  ]
end

# Calculate the Markov transition matrix from a sequence of ordinal
# symbols. The entry in the ith row and jth column is a_{i, j}.
function get_markov_matrix(
  ordinal_symbols::Vector{Int},
  m::Int
)::Matrix{Int}
  n = factorial(m)
  matrix = zeros(Int, n, n)
  for i in 1:length(ordinal_symbols)-1
    matrix[ordinal_symbols[i], ordinal_symbols[i+1]] += 1
  end
  return matrix
end

# Calculate p_i, the stationary probability of the ith ordinal symbol.
function ordinal_probability(markov_matrix::Matrix{Int}, i::Int)::Float64
  return sum(markov_matrix[i, :]) / sum(markov_matrix)
end

# Calculate weighted entropy for a particular ordinal symbol timeseries.
function weighted_entropy(
  N::Int, # Length of original time series.
  w::Int, # Window lag for original time series.
  m::Int, # Samples per window for original time series.
  τ::Int, # Sampling delay for original time series.
  ordinal_symbols2::Vector{Int},
  m2::Int # Samples per window for this ordinal symbol's time series.
)::Float64
  # O is the number of occurrences of the ordinal symbol.
  O = length(ordinal_symbols2)
  # T is the total number of windows into which the original
  # time series is divided.
  T = floor(float(N - (m - 1)*τ - 1) / float(w)) + 1
  # K is the probability weighting for each ordinal symbol.
  K = O / T

  markov_matrix = get_markov_matrix(ordinal_symbols2, m2)

  ordinal_probabilities = [
    K * ordinal_probability(markov_matrix, i)
    for i in 1:size(markov_matrix, 1)
  ]
  return -sum(
    ordinal_probabilities[i] == 0.0
      ? 0.0
      : ordinal_probabilities[i] * log2(ordinal_probabilities[i])
    for i in 1:size(markov_matrix, 1)
  )
end

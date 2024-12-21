include("./power_series.jl")

function determinant(A::Array{Integer, 3})::Vector{Integer}
  # Calculate determinant of A using Laplace expansion.
  m, n, K = size(A)
  @assert m == n
  if m == 1
    return A[1, 1, :]
  end
  det = zeros(Integer, K)
  sgn = 1
  for i in 1:m
    add!(det, scale(sgn, determinant(minor(A, i, 1))))
    i *= -1
  end
  return det
end

function minor(A::Array{Integer, 3}, i, j)::Array{Integer, 3}
  # View first minor of matrix without i row and j column.
  rows, cols, _ = size(A)
  
  row_idxs = vcat(1:i-1, i+1:rows)
  col_idxs = vcat(1:j-1, j+1:cols)
  return @view A[row_idxs, col_idxs, :]
end
include("./power_series.jl")

function determinant(
  A::Array{Integer, 3},
  increasing_first_lap=true
)::Vector{Integer}
  # Calculate determinant of A using Laplace expansion.
  m, n, K = size(A)
  @assert m == n
  if m == 1
    return A[1, 1, :]
  end
  det = Integer[]
  sgn = 1
  for i in 1:m
    det = add(
      det,
      scale(
        sgn,
        multiply(
          A[1, i, :],
          determinant(minor(A, 1, i))
        )
      )
    )
    sgn *= -1
  end

  rational_factor = ones(Integer, K)
  if !increasing_first_lap
    for k in 2:2:K
      rational_factor[k] = -1
    end
  end
  det = multiply(det, rational_factor)[1:K]

  return det
end

function minor(A::Array{Integer, 3}, i, j)::Array{Integer, 3}
  # View first minor of matrix without i row and j column.
  rows, cols, _ = size(A)
  
  row_idxs = vcat(1:i-1, i+1:rows)
  col_idxs = vcat(1:j-1, j+1:cols)
  return @view A[row_idxs, col_idxs, :]
end
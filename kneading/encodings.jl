function coefficient_encoding(coeff)
  return abs(coeff) * 2 - (coeff < 0 ? 1 : 0)
end

function determinant_encoding(det)
  K = length(det) - 1
  coeffs_encoded = Vector{Integer}(undef, K + 1)
  for k in 1:K+1
    coeffs_encoded[k] = coefficient_encoding(det[k])
  end
  alphabet_size = maximum(coeffs_encoded) + 1
  
  encoding = 0
  power = 1
  for k in 1:K+1
    encoding += coeffs_encoded[k] / alphabet_size^power
    power += 1
  end

  return encoding
end

function matrix_encoding(matrix)
  m = size(matrix, 1)
  n = size(matrix, 2)
  K = size(matrix, 3) - 1

  encoding = 0
  power = 1
  for k in 1:K+1
    for i in 1:m
      for j in 1:n
        encoding += coefficient_encoding(matrix[i, j, k]) * 5^power
        power += 1
      end
      power += 1
    end
    power += 1
  end

  return encoding
end
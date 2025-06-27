# A double tent map on [0, 1].
# The map is piecewise linear on [0, 0.25], [0.25, 0.75], and [0.75, 1].
# f(0) = 0, f(1) = 1.
# The parameters p = [a, b] are the heights at x = 0.25 and x = 0.75, respectively.
# So, f(0.25) = a and f(0.75) = b.

# Define the map.
function double_tent_map(p, x, n_iter=1)
    a, b = p
    for _ in 1:n_iter
        if 0 <= x < 0.25
            x = 4 * a * x
        elseif 0.25 <= x < 0.75
            slope = (b - a) / 0.5
            x = slope * (x - 0.25) + a
        elseif 0.75 <= x <= 1
            x = (4 - 4 * b) * (x - 0.75) + b
        end
    end
    return x
end

# Define the critical points of the map.
function critical_points(p)
  return [0.25, 0.75]
end

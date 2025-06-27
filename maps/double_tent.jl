# A double tent map on [0, 1].
# The map is piecewise linear on [0, 1/3], [1/3, 2/3], and [2/3, 1].
# f(0) = 0, f(1) = 1.
# The parameters p = [a, b] are the heights at x = 1/3 and x = 2/3, respectively.
# So, f(1/3) = a and f(2/3) = b.

# Define the map.
function double_tent_map(p, x, n_iter=1)
    a, b = p
    for _ in 1:n_iter
        if 0 <= x < 1/3
            x = 3 * a * x
        elseif 1/3 <= x < 2/3
            slope = 3 * (b - a)
            x = slope * (x - 1/3) + a
        elseif 2/3 <= x <= 1
            slope = 3 * (1 - b)
            x = slope * (x - 2/3) + b
        end
    end
    return x
end

# Define the critical points of the map.
function critical_points(p)
  return [1/3, 2/3]
end

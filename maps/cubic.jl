# Cubic with a root at x = 0.
# Two parameters:
# p[1] = negative root.
# p[2] = positive root.
function map(p, x)
  return -x * (x - p[1]) * (x - p[2])
end

function derivative(p, x)
  return -3 * x^2 - 2 * (p[1] + p[2]) * x + p[1] * p[2]
end

function critical_points(p)
  a = p[1]
  b = p[2]
  average = (a+b)/3
  difference = sqrt((a+b)^2-3*a*b)/3
  return [
    average - difference,
    average + difference
  ]
end

is_continuous = true
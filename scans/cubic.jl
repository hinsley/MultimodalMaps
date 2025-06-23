using Pkg
Pkg.activate(".")

using Plots
using Random

include("../kneading/encodings.jl")
include("../kneading/matrix.jl")
include("../maps/cubic.jl")

# Define parameter ranges.
ε = 1.4
a_vals = range(-1.8, stop=-ε, length=1000)
b_vals = range(ε, stop=1.8, length=1000)

# a_vals = range(-1.83, stop=-1.78, length=10)
# b_vals = range(1.32, 1.43, length=10)

# a_vals = range(-1.43, stop=-1.32, length=5)
# b_vals = range(1.78, stop=1.83, length=5)

# Allocate a matrix to store encoding values for each (a, b) pair.
Z = zeros(length(a_vals), length(b_vals))

# Sweep over parameter space.
# @gif for iterates in 1:20
  iterates = 20
  @time for i in 1:length(a_vals)
    for j in 1:length(b_vals)
      a = a_vals[i]
      b = b_vals[j]
      p = [a, b]

      # Compute kneading matrix for given parameters.
      crit_points = critical_points(p)
      matrix = allocate_kneading_matrix(crit_points, iterates)
      kneading_matrix!(matrix, map, crit_points, p)
      # println("Matrix for parameters a=$a, b=$b:")
      # for k in 1:size(matrix,3)
      #     println("Layer $k:")
      #     display(matrix[:,:,k])
      #     println()
      # end
      # println("----------------------------------------")

      # Get encoding.
      encoding = matrix_encoding(matrix)
      # println(encoding)
      # Random.seed!(encoding)
      # Z[j, i] = rand()
      Z[j, i] = encoding
    end
  end

  # Create a heatmap with exact region boundaries.
  heatmap(
    a_vals, b_vals, Z,
    xlabel="a", ylabel="b", 
    title="Matrix encoding with region boundaries: f(x)=-x(x-a)(x-b)",
    color=:thermal,#osloS,
    size=(1200, 1200), # Increase plot resolution.
    dpi=300, # Increase DPI for better quality.
    aspect_ratio=:equal # Make axes symmetric.
  )

  # Find all unique values in Z to draw exact boundaries.
  # We only take every second level to avoid visual clutter.
  unique_levels = sort(unique(vec(Z)))
  contour!(
    a_vals, b_vals, Z,
    levels=unique_levels,
    color=:white,
    linewidth=1
  )
# end fps=4

# Optionally save the figure.
savefig("matrix_encoding_heatmap_with_boundaries.png")

include("../kneading/determinant.jl")
# @gif for iterates in 1:20
  iterates = 20
  @time for i in 1:length(a_vals)
    for j in 1:length(b_vals)
      a = a_vals[i]
      b = b_vals[j]
      p = [a, b]

      # Compute kneading matrix for given parameters.
      crit_points = critical_points(p)
      matrix = allocate_kneading_matrix(crit_points, iterates)
      kneading_matrix!(matrix, map, crit_points, p)

      # Convert the matrix entries to Integers so they can
      # be used in the determinant calculation.
      matrix = convert(Array{Integer, 3}, matrix)
      # println("Matrix for parameters a=$a, b=$b:")
      # for k in 1:size(matrix,3)
      #     println("Layer $k:")
      #     display(matrix[:,:,k])
      #     println()
      # end
      # println("----------------------------------------")

      # Calculate determinant.
      det = determinant(matrix[:, 2:end, :], false)

      # Get encoding.
      encoding = determinant_encoding(det)
      # println(det)
      # println(encoding)
      Random.seed!(Int(floor(encoding * typemax(Int))))
      Z[j, i] = rand()
      # Z[j, i] = encoding
    end
  end

  # Create heatmap.
  heatmap(
    a_vals, b_vals, Z,
    xlabel="a", ylabel="b",
    title="Determinant encoding: f(x)=-x(x-a)(x-b)",
    color=:thermal,#osloS,
    size=(1200, 1200), # Increase plot resolution.
    dpi=300, # Increase DPI for better quality.
    aspect_ratio=:equal # Make axes symmetric.
  )
# end fps=3

# Optionally save the figure.
savefig("determinant_encoding_heatmap.png")

include("../kneading/smallest_root.jl")
  iterates = 30
  @time for i in 1:length(a_vals)
    for j in 1:length(b_vals)
      a = a_vals[i]
      b = b_vals[j]
      p = [a, b]

      # Compute kneading matrix for given parameters.
      crit_points = critical_points(p)
      matrix = allocate_kneading_matrix(crit_points, iterates)
      kneading_matrix!(matrix, map, crit_points, p)

      # Convert the matrix entries to Integers so they can
      # be used in the determinant calculation.
      matrix = convert(Array{Integer, 3}, matrix)
      # println("Matrix for parameters a=$a, b=$b:")
      # for k in 1:size(matrix,3)
      #     println("Layer $k:")
      #     display(matrix[:,:,k])
      #     println()
      # end
      # println("----------------------------------------")

      # Calculate determinant.
      det = determinant(matrix[:, 2:end, :], true)

      # Get smallest root.
      root = smallest_root(convert(Vector{Int}, det))
      Z[j, i] = root < 1 ? log(1/root) : 0
    end
  end

  # Create heatmap.
  heatmap(
    a_vals, b_vals, Z,
    xlabel="a", ylabel="b",
    title="Topological entropy: f(x)=-x(x-a)(x-b)",
    color=:thermal,
    size=(1200, 1200), # Increase plot resolution.
    dpi=300, # Increase DPI for better quality.
    aspect_ratio=:equal # Make axes symmetric.
  )

savefig("topological_entropy_heatmap.png")

# Calculate maximum Lyapunov exponent.
iterates = 50
burn_in = 10
@time for i in 1:length(a_vals)
  for j in 1:length(b_vals)
    a = a_vals[i]
    b = b_vals[j]
    p = [a, b]
    
    crit_points = critical_points(p)
    lyap_exponents = []

    for crit_point in crit_points
      x = crit_point
      # Burn in.
      for k in 1:burn_in
        x = map(p, x)
      end

      # Calculate exponent.
      sum_val = 0
      for k in 1:iterates
        sum_val += log(abs(derivative(p, x)))
        x = map(p, x)
      end
      push!(lyap_exponents, sum_val / iterates)
    end
    Z[j, i] = max(lyap_exponents...)
  end
end

# Create heatmap.
heatmap(
  a_vals, b_vals, Z,
  xlabel="a", ylabel="b",
  title="Maximum Lyapunov Exponent: f(x)=-x(x-a)(x-b)",
  color=:thermal,
  size=(1200, 1200), # Increase plot resolution.
  dpi=300, # Increase DPI for better quality.
  aspect_ratio=:equal # Make axes symmetric.
)

# Optionally save the figure.
savefig("lyapunov_exponent_heatmap.png")
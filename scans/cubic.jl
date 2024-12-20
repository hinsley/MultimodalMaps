using Pkg
Pkg.activate(".")

using Plots
using Random

include("../kneading/encodings.jl")
include("../kneading/matrix.jl")
include("../maps/cubic.jl")

# Define parameter ranges.
ε = 0.8
a_vals = range(-2, stop=-ε, length=2000)
b_vals = range(ε, stop=2, length=2000)

# Allocate a matrix to store encoding values for each (a, b) pair.
Z = zeros(length(a_vals), length(b_vals))

# Sweep over parameter space.
@gif for iterates in 1:20
  # iterates = 2
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
          Random.seed!(encoding)
          Z[j, i] = rand()
      end
  end

  # Create a heatmap.
    heatmap(
      a_vals, b_vals, Z,
      xlabel="a", ylabel="b", 
      title="f(x)=-x(x-a)(x-b)",
      color=:osloS,
      size=(1200, 1200), # Increase plot resolution
      dpi=300, # Increase DPI for better quality
      aspect_ratio=:equal # Make axes symmetric
  )
end fps=4

# Optionally save the figure
savefig("matrix_encoding_heatmap.png")
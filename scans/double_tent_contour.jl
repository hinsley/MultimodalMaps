using Pkg
Pkg.activate(".")

using Plots
using Random

include("../kneading/encodings.jl")
include("../kneading/matrix.jl")
include("../kneading/determinant.jl")
include("../scans/contours.jl")
include("../maps/double_tent.jl")

# Define parameter ranges.
a_vals = range(0.5, stop=1, length=2000)
b_vals = range(0, stop=0.5, length=2000)

# Allocate a matrix to store encoding values for each (a, b) pair.
Z = zeros(length(a_vals), length(b_vals))

# Choose the scan type: :matrix or :determinant.
scan_type = :matrix

# Choose the number of iterates.
iterates = 50

# Choose the color exponent for separation of iterates.
color_exp = 2

# Calculate kneading diagram.
fig = plot(
  aspect_ratio=:equal,
  colorbar=false,
  xlims=(minimum(a_vals), maximum(a_vals)),
  ylims=(minimum(b_vals), maximum(b_vals)),
  xlabel="a",
  ylabel="b",
  legend=false,
  size=(1000, 1000),
  xguidefontsize=14,
  yguidefontsize=14,
  xtickfontsize=12,
  ytickfontsize=12,
  left_margin=3Plots.mm,
  bottom_margin=3Plots.mm,
  right_margin=3Plots.mm
)
for iterate in iterates:-1:2
  @time for i in 1:length(a_vals)
    for j in 1:length(b_vals)
      a = a_vals[i]
      b = b_vals[j]
      p = [a, b]
      
      # Compute kneading matrix.
      crit_points = critical_points(p)
      matrix = allocate_kneading_matrix(crit_points, iterate)
      kneading_matrix!(matrix, double_tent_map, crit_points, p)

      if scan_type == :matrix
        encoding = matrix_encoding(matrix)
      elseif scan_type == :determinant
        matrix = convert(Array{Integer, 3}, matrix)
        det = determinant(matrix[:, 2:end, :], false)
        encoding = determinant_encoding(det)
      end
      Z[j, i] = encoding
    end
  end

  # Create a contour plot.
  contour_xs, contour_ys = march_squares_simple(
    Z,
    a_vals,
    b_vals
  )
  plot!(
    fig,
    contour_xs,
    contour_ys,
    color=RGB(
      ((iterate-2)/iterates)^(1/color_exp),
      ((iterate-2)/iterates)^(1/color_exp),
      ((iterate-2)/iterates)^(1/color_exp)
    ),
    linewidth=1,
    label="Iterate $iterate"
  )
end

# Add a title to the plot.
title!(
  fig,
  "Kneading diagram: $(uppercasefirst(string(scan_type))) encoding"
)

# Display the plot.
display(fig)

# Save the plot to a file.
savefig(fig, "kneading_diagram_$(uppercasefirst(string(scan_type)))_$(iterates).png")
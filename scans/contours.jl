function march_squares_simple(
  Z::Matrix{Float64},
  x_vals::AbstractVector{Float64},
  y_vals::AbstractVector{Float64}
)::Tuple{Vector{Float64}, Vector{Float64}}
  """
  March over the matrix Z, and return the x and y coordinates of the
  contours.
  Splits line segments via NaN values instead of returning
  multiple line segments to make plotting more efficient.
  """

  # Initialize the list of contours.
  contour_xs = []
  contour_ys = []

  # Iterate over the squares defined by the matrix Z.
  for i in 1:size(Z, 1)-1
    for j in 1:size(Z, 2)-1
      # Get the four corners of the square.
      x_tl, y_tl = x_vals[i], y_vals[j]
      x_tr, y_tr = x_vals[i+1], y_vals[j]
      x_br, y_br = x_vals[i+1], y_vals[j+1]
      x_bl, y_bl = x_vals[i], y_vals[j+1]

      # Get the locations of the side midpoints.
      x_tm, y_tm = (x_tl + x_tr) / 2, (y_tl + y_tr) / 2
      x_rm, y_rm = (x_tr + x_br) / 2, (y_tr + y_br) / 2
      x_bm, y_bm = (x_br + x_bl) / 2, (y_br + y_bl) / 2
      x_lm, y_lm = (x_bl + x_tl) / 2, (y_bl + y_tl) / 2

      # Get the values of the corners.
      z_tl = Z[i, j]
      z_tr = Z[i+1, j]
      z_br = Z[i+1, j+1]
      z_bl = Z[i, j+1]

      # Determine contour type.
      # Both constant squares and double-diagonal contour squares are left
      # undrawn, the latter based on the assumption that the resolution of the
      # plot is sufficiently fine to avoid double-diagonals.
      if z_tl != z_tr == z_br == z_bl # Diagonal contour at top left.
        append!(contour_xs, [
          x_tm,
          x_lm,
          NaN
        ])
        append!(contour_ys, [
          y_tm,
          y_lm,
          NaN
        ])
      elseif z_tr != z_tl == z_br == z_bl # Diagonal contour at top right.
        append!(contour_xs, [
          x_rm,
          x_tm,
          NaN
        ])
        append!(contour_ys, [
          y_rm,
          y_tm,
          NaN
        ])
      elseif z_br != z_tr == z_tl == z_bl # Diagonal contour at bottom right.
        append!(contour_xs, [
          x_bm,
          x_rm,
          NaN
        ])
        append!(contour_ys, [
          y_bm,
          y_rm,
          NaN
        ])
      elseif z_bl != z_br == z_tr == z_tl # Diagonal contour at bottom left.
        append!(contour_xs, [
          x_lm,
          x_bm,
          NaN
        ])
        append!(contour_ys, [
          y_lm,
          y_bm,
          NaN
        ])
      elseif z_tl == z_tr != z_bl == z_br # Horizontal contour.
        append!(contour_xs, [
          x_lm,
          x_rm,
          NaN
        ])
        append!(contour_ys, [
          y_lm,
          y_rm,
          NaN
        ])
      elseif z_tl == z_bl != z_tr == z_br # Vertical contour.
        append!(contour_xs, [
          x_tm,
          x_bm,
          NaN
        ])
        append!(contour_ys, [
          y_tm,
          y_bm,
          NaN
        ])
      end
    end
  end

  # Discard the last NaN values and return the contour polyisolines.
  return contour_xs[1:end-1], contour_ys[1:end-1]
end

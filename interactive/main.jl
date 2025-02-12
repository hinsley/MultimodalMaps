using Pkg
Pkg.activate(".")

using GLMakie, DifferentialEquations

include("../systems/lorenz.jl") # Import system definition.

# Initialize system parameters from module.
sys = LorenzSystem
u0 = sys.initial_condition
p0 = sys.default_parameters
opn_params = sys.default_opn_params
param_ranges = sys.parameter_ranges

# Create main figure.
fig = Figure(resolution=(1600, 1200))

# Create 3D axis for trajectory plot.
ax3d = Axis3(fig[1, 1],
  title="System Trajectory",
  xlabel="x", ylabel="y", zlabel="z"
)

# Create 2D axis for return map plot.
ax2d = Axis(fig[1, 2],
  title="Return Map",
  xlabel="xₖ", ylabel="xₖ₊₁",
  aspect=DataAspect()
)

# Create control panel.
settings = GridLayout(fig[2, 1:2])
buttongrid = GridLayout(settings[1, 1])

# Parameter controls.
param_controls = [
  (; label="Trajectory 1 Time",    range=100:100:10000,    start=1000,          type=:float),
  (; label="Trajectory 2 Time",    range=0:100:100000,     start=10000,         type=:float),
  (; label="Transient Time",       range=100:100:1000,     start=200,           type=:float),
  (; label="Time Step",            range=1e-4:1e-4:1,      start=1e-2,          type=:float),
  (; label="m1 (Ordinal Pattern)", range=1:50,             start=opn_params[1], type=:natural),
  (; label="n1",                   range=1:50,             start=opn_params[3], type=:natural),
  (; label="l1",                   range=1:50,             start=opn_params[5], type=:natural),
  (; label="m2",                   range=1:50,             start=opn_params[2], type=:natural),
  (; label="n2",                   range=1:50,             start=opn_params[4], type=:natural),
  (; label="l2",                   range=1:50,             start=opn_params[6], type=:natural),
  (; label="Ordinal Function:",    range=nothing,          start="u[3]",        type=:expr),
  (; label="Ordinal Symbol Index", range=1:99999,          start=1,             type=:natural),
  (; label="Return Coordinate:",   range=nothing,          start="u[3]",        type=:expr),
  (; label="Section Condition:",   range=nothing,          start="true",        type=:expr)
]

controls = []
# Split controls into two columns.
for (i, pc) in enumerate(param_controls)
  if i <= 10  # Original parameters in first column
    row = i
    lbl_col = 1
    ctrl_col = 2
  else  # New parameters in second column
    row = i - 10
    lbl_col = 3
    ctrl_col = 4
  end
  
  lbl = Label(buttongrid[row, lbl_col], pc.label)
  
  if pc.type == :slider
    control = Slider(buttongrid[row, ctrl_col], range=pc.range, startvalue=pc.start)
  elseif pc.type == :float
    control = Textbox(buttongrid[row, ctrl_col], placeholder=string(pc.start), width=240)
    on(control.stored_string) do s
      val = tryparse(Float64, s)
      if !isnothing(val) && val in pc.range
        # Valid input
      else
        @async begin
          sleep(0.1)
          control.stored_string[] = string(pc.start)
        end
      end
    end
  elseif pc.type == :natural
    control = Textbox(buttongrid[row, ctrl_col], placeholder=string(pc.start), width=240)
    on(control.stored_string) do s
      val = tryparse(Int, s)
      if !isnothing(val) && val in pc.range
        # Valid input.
      else
        @async begin
          sleep(0.1)
          control.stored_string[] = string(pc.start)
        end
      end
    end
  elseif pc.type == :expr
    control = Textbox(buttongrid[row, ctrl_col], placeholder=string(pc.start), width=240)
  end
  push!(controls, control)
end

# Prepopulate controls with default values.
for (i, control) in enumerate(controls)
  if param_controls[i].type != :slider
    control.stored_string[] = string(param_controls[i].start)
  end
end

# Set fixed widths for labels and controls.
buttongrid.colsizes[1] = Fixed(150)  # First column labels.
buttongrid.colsizes[2] = Fixed(240)  # First column controls.
buttongrid.colsizes[3] = Fixed(150)  # Second column labels.
buttongrid.colsizes[4] = Fixed(240)  # Second column controls.

# Solve button.
# Place solve button at bottom spanning both columns.
solve_btn = Button(buttongrid[11, 1:4], label="Solve & Update")

# Storage for plot elements.
traj_lines = []
scatter_plots = []
return_scatters = []

include("../return_maps/return_map.jl")
include("../return_maps/auto_coordinates.jl")

# Create function generators for the dynamic expressions.
function create_expr_function(expr_str::String)::Function
  parsed = Meta.parse(expr_str)
  f = @eval (u) -> $parsed
  # Force compilation by calling the function once with a dummy input.
  try
    Base.invokelatest(f, zeros(3))
  catch e
    Base.invokelatest(f, zeros(6)) # In case of a distance function.
  end
  return (u) -> Base.invokelatest(f, u)
end

# Update plots function.
function update_plots()
  # Clear previous plots.
  empty!(ax3d)
  empty!(ax2d)
  
  # Get current parameter values.
  trajectory_time1 = tryparse(Float64, controls[1].stored_string[]) |> v -> isnothing(v) ? param_controls[1].start : v
  trajectory_time2 = tryparse(Float64, controls[2].stored_string[]) |> v -> isnothing(v) ? param_controls[2].start : v
  transient_time = tryparse(Float64, controls[3].stored_string[]) |> v -> isnothing(v) ? param_controls[3].start : v
  timestep = tryparse(Float64, controls[4].stored_string[]) |> v -> isnothing(v) ? param_controls[4].start : v
  opn_params = Int[tryparse(Int, s.stored_string[]) |> v -> isnothing(v) ? param_controls[i+4].start : v for (i,s) in enumerate(controls[5:10])]
  
  # Create compiled functions from expressions.
  ordinal_fn = create_expr_function(controls[11].stored_string[])
  ranked_ordinal_symbol_index = tryparse(Int, controls[12].stored_string[]) |> v -> isnothing(v) ? 1 : v
  return_coord_string = controls[13].stored_string[]
  # Check if return coordinate function string starts with "auto".
  if startswith(return_coord_string, "auto")
    # Split string into words and parse parameters.
    parts = split(return_coord_string)
    if length(parts) >= 4
      threshold = parse(Float64, parts[2])
      truncation = parse(Int64, parts[3])
      distance_fn = create_expr_function(String(parts[4]))
      # Store both parameters as a tuple.
      return_coord = (threshold, truncation, distance_fn)
    else
      # Throw an error.
      error("Invalid auto return coordinate string. Should be Douglas-Peucker threshold, truncation quantity, distance function.")
    end
  else
    # Create compiled function from expression.
    return_coord = create_expr_function(controls[13].stored_string[])
  end
  section_condition_fn = create_expr_function(controls[14].stored_string[])
  
  # Solve ODE.
  prob = ODEProblem(sys.derivatives!, u0, (0.0, trajectory_time1), p0)
  sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)
  
  # Plot 3D trajectory (after transient time).
  if trajectory_time2 > 0.0
    t_idx = findfirst(t -> t >= transient_time, sol.t)
    lines!(ax3d, sol[1,t_idx:end], sol[2,t_idx:end], sol[3,t_idx:end], color=:blue)
  end
  
  # Calculate return map using OPN parameters.
  return_vals, return_points = calculate_return_itinerary(
    sys.derivatives!,
    p0,
    u0,
    transient_time,
    timestep,
    (0.0, trajectory_time1),
    (0.0, trajectory_time2),
    ordinal_fn,
    return_coord,
    section_condition_fn,
    opn_params...;
    ranked_ordinal_symbol_index=ranked_ordinal_symbol_index,
    return_points=true
  )
  
  # Plot return map points in state space.
  scatter!(
    ax3d,
    [point[1] for point in return_points],
    [point[2] for point in return_points],
    [point[3] for point in return_points],
    color=:red,
    markersize=3
  )

  # If automatic coordinates used for return map, plot coordinate curve.
  if trajectory_time2 == 0.0 && length(return_coord) > 1
    return_curve = auto_coordinates(
      return_points,
      return_coord[1],
      return_coord[3],
      return_coord[2],
      return_curve=true
    )
    lines!(
      ax3d,
      [point[1] for point in return_curve],
      [point[2] for point in return_curve],
      [point[3] for point in return_curve],
      color=:green,
      linewidth=3
    )
  end
  
  # Plot return map.
  scatter!(ax2d, return_vals[1:end-1], return_vals[2:end], color=:red, markersize=5)
  
  # Set axis limits.
  xmin, xmax = extrema(return_vals)
  lines!(ax2d, [xmin, xmax], [xmin, xmax], color=:black, linestyle=:dash)
  
  autolimits!(ax3d)
  autolimits!(ax2d)
end

# Connect solve button - use lift instead of on to prevent multiple handlers
lift(solve_btn.clicks) do _
    update_plots()
end

# Initial plot.
update_plots()

# Display figure.
display(fig)
using Pkg
Pkg.activate(".")

using DifferentialEquations
using GLMakie

# 3D Morris-Lecar model from https://www.izhikevich.org/publications/nesb.pdf
# Equations 25, 26 and parameters from Figure 59.
function ml3_derivatives!(du, u, p, t)
  V, w, z = u
  gl, gk, gCa, El, Ek, ECa, μ, V1, V2, V3, V4 = p
  m∞(V) = 0.5 * (1 + tanh((V - V1) / V2))
  w∞(V) = 0.5 * (1 + tanh((V - V3) / V4))
  λ(V) = 1.0 / 3.0 * cosh((V - V3) / 2.0 / V4)
  du[1] = -z - gl * (V - El) - gk * w * (V - Ek) - gCa * m∞(V) * (V - ECa)
  du[2] = λ(V) * (w∞(V) - w)
  du[3] = μ * (0.2 + V)
end

u0 = [-0.3, 0.1, 0.1]
p = [
  0.5,   # gl
  2.0,   # gk
  1.2,   # gCa
  -0.5,  # El
  -0.7,  # Ek
  1.0,   # ECa
  0.005, # μ
  -0.01, # V1
  0.15,  # V2
  0.1,   # V3
  0.05   # V4
]

prob = ODEProblem(ml3_derivatives!, u0, (0.0, 1e5), p)
sol = solve(prob, Tsit5(), abstol=1e-8, reltol=1e-8)

fig = Figure()
ax = Axis3(fig[1,1], xlabel="V", ylabel="w", zlabel="z")

xs = [u[1] for u in sol.u]
ys = [u[2] for u in sol.u]
zs = [u[3] for u in sol.u]

lines!(ax, xs, ys, zs, color=:blue)

display(fig)

# Plot voltage timeseries
fig2 = Figure()
ax2 = Axis(fig2[1,1], xlabel="t", ylabel="V")
lines!(ax2, sol.t, xs, color=:blue)
display(fig2)

# Find the equilibrium of the system.
gl, gk, gCa, El, Ek, ECa, μ, V1, V2, V3, V4 = p
m∞(V) = 0.5 * (1 + tanh((V - V1) / V2))
w∞(V) = 0.5 * (1 + tanh((V - V3) / V4))
λ(V) = 1.0 / 3.0 * cosh((V - V3) / 2.0 / V4)
V0 = -0.2
w0 = (1.0 + tanh(-6.0)) / 2.0
z0 = -gl * (V0 - El) - gk * w0 * (V0 - Ek) - gCa * m∞(V0) * (V0 - ECa)

# Compute the eigenvalues of the Jacobian of the system at the equilibrium.
using ForwardDiff

# Define the system as a function for autodiff
function ml3_system(u)
    V, w, z = u
    [
        -z - gl * (V - El) - gk * w * (V - Ek) - gCa * m∞(V) * (V - ECa),
        λ(V) * (w∞(V) - w),
        μ * (0.2 + V)
    ]
end

# Compute Jacobian at equilibrium point
J = ForwardDiff.jacobian(ml3_system, [V0, w0, z0])

# Calculate eigenvalues
using LinearAlgebra
eigenvals = eigvals(J)
println("Eigenvalues at equilibrium: ", eigenvals)
# Calculate and print eigenvectors
eigenvecs = eigvecs(J)
println("\nEigenvectors at equilibrium:")
for i in 1:3
    println("\nEigenvector $i (corresponding to λ = $(eigenvals[i])):")
    println(eigenvecs[:,i])
end


# Plot eigenvectors from equilibrium point in 3D
scale_factor = 1e-1  # Adjust this to make eigenvectors more/less visible
for i in 1:3
    # Get real parts of eigenvalue and eigenvector
    lambda = real(eigenvals[i])
    v = real(eigenvecs[:,i])
    
    # Normalize eigenvector and scale by eigenvalue magnitude
    v = v / norm(v) * scale_factor
    
    # Choose color based on eigenvalue sign
    arrow_color = lambda > 0 ? :red : :blue
    
    # Plot eigenvector as arrow from equilibrium point
    arrows!(ax, 
        [V0], [w0], [z0],  # Start point (equilibrium)
        [v[1]], [v[2]], [v[3]],  # Direction and magnitude
        color=arrow_color,
        linewidth=1e-3,
        arrowsize=3e-3
    )
end

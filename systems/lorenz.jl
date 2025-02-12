module LorenzSystem
export derivatives!, default_parameters, parameter_ranges, default_opn_params, initial_condition

function derivatives!(du, u, p, t)
    x, y, z = u
    σ, ρ, β = p
    du[1] = σ*(y - x)
    du[2] = x*(ρ - z) - y
    du[3] = x*y - β*z
end

const default_parameters = [10.0, 75.9, 8/3]
const parameter_ranges = [(0.1, 50.0), (0.1, 50.0), (0.1, 10.0)]
const default_opn_params = (3, 4, 2, 1, 3, 1)
const initial_condition = [1.0, 1.0, 1.0]

end
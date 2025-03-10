using JuMP
using Ipopt

# model = Model(Ipopt.Optimizer)
model = Model()
set_optimizer(model, Ipopt.Optimizer)

@variable(model, x >= 0)
@variable(model, 0 <= y <= 30)

@objective(model, Min, 12x + 20y)

@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)

# for i in 1:3
#     @constraint(model, 6x + 4y >= 5i)
# end

optimize!(model)

if !is_solved_and_feasible(model)
    error("Solver did not find an optimal solution")
end

println("optimal value : ", objective_value(model))
println("optimal x : ", value(x), ", optimal y : ", value(y))

# https://jump.dev/JuMP.jl/stable/tutorials/getting_started/getting_started_with_JuMP/
# varialbe basics - sparseaxisarrays
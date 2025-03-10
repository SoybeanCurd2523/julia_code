# all variables must be integer, and they must typically have finite bounds

using JuMP
using HiGHS

model = Model(HiGHS.Optimizer)
set_silent(model)
@variable(model, 1 <= x[1:4] <= 4, Int)

# The MOI.AllDifferent set ensures that every element in a list takes a different integer value.
@constraint(model, x in MOI.AllDifferent(4))
optimize!(model)
is_solved_and_feasible(model)
print(value.(x))

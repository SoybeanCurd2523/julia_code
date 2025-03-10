using JuMP
import HiGHS

# standard form example
# https://en.wikipedia.org/wiki/Basic_feasible_solution

# A = [1 5 3 4 6; 0 1 3 5 6]
# b = [14, 7]

# n = size(A, 2)
# model = Model(HiGHS.Optimizer)
# set_silent(model)
# @variable(model, x[1:n] >= 0)
# @constraint(model, A * x == b)
# optimize!(model)
# is_solved_and_feasible(model)

# println("optimal value of x : ", value.(x))

# # basic 변수는 최적 해를 구성하는 데 사용되고, nonbasic 변수는 상한 또는 하한에 고정
# get_attribute.(x, MOI.VariableBasisStatus())

# # indices : index 의 복수형
# indices = get_attribute.(x, MOI.VariableBasisStatus()) .== MOI.BASIC

# B = A[:, indices]
# println(B\b)
# println(value.(x[indices]))


# A more complicated example
model = Model(HiGHS.Optimizer)
set_silent(model)

@variable(model, x >= 0)
@variable(model, 0 <= y <= 3)
@variable(model, z <= 1)

@objective(model, Min, 12x + 20y - z)

@constraint(model, c1, 6x + 8y >= 100)
@constraint(model, c2, 7x + 12y >= 120)
@constraint(model, c3, x + y <= 20)

optimize!(model)
is_solved_and_feasible(model)

get_attribute.(x, MOI.VariableBasisStatus())
get_attribute.(y, MOI.VariableBasisStatus())
get_attribute.(z, MOI.VariableBasisStatus())

# 모든 변수 확인
# all_variables(model)

########### Computing the duals of a mixed-integer program
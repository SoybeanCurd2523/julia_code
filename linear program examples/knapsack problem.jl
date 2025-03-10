#=
using JuMP
import HiGHS

n = 5
capacity = 10.0

profit = [5.0, 3.0, 2.0, 7.0, 4.0]
weight = [2.0, 8.0, 4.0, 2.0, 5.0]

model = Model(HiGHS.Optimizer)

@variable(model, x[1:n], Bin) # binary 변수. 0 또는 1

weight_sum = 0
for i in 1:n
    weight_sum += weight[i] * x[i]
end
@constraint(model, weight_sum <= capacity)
# @constraint(model, sum(weight[i] * x[i] for i in 1:n) <= capacity)

profit_sum = 0
for i in 1:n
    profit_sum += profit[i] * x[i]
end

@objective(model, Max, profit_sum)

print(model)

optimize!(model)

solution_summary(model)


item_chosen = [] # 빈 리스트 설정

for i in 1:n
    if value(x[i]) > 0.5
        push!(item_chosen, i) # push!는 기존 배열에 요소 추가(수정함)
    end
end

# items_chosen = [i for i in 1:n if value(x[i]) > 0.5]
=#

######################### 

# 함수로 만들기
using JuMP
using HiGHS

function solve_knapsack_problem(profit, weight, capacity)

    n = length(weight)
    @assert length(profit) == n # @assert는 디버깅과 오류 확인을 위해 사용되는 Julia의 매크로입니다. 조건이 false일 경우, 에러를 발생시킵니다.

    model = Model(HiGHS.Optimizer)
    set_silent(model) # 모델이 실행 중에 출력 메시지를 숨깁니다
    @variable(model, x[1:n], Bin)
    @objective(model, Max, profit' * x)
    @constraint(model, weight' * x <= capacity)

    optimize!(model)
    if !is_solved_and_feasible(model)
        error("Solver did not find an optimal solution")
    end
    println("objective is : ", objective_value(model))
    println("solution is " )
    for i in 1:n
        println("x[$i] : ", value(x[i]))
    end
    chosen_items = [i for i in 1:n if value(x[i]) > 0.5]
    return chosen_items
end

solve_knapsack_problem([5.0, 3.0, 2.0, 7.0, 4.0], [2.0, 8.0, 4.0, 2.0, 5.0], 10)


# indicator constraint
# 이진변수 z가 1이면 제약조건 활성화, 0이면 비활성화
# @constraint(model, indicator_constraint, z => {x <= 10})

# big-M trick : indicator constraint 대체
# M은 매우 큰 상수
# @constraint(model, g(x) <= 1 + M * (1 - z))
# z가 1이면 제약 조건 적용
# z가 0이면 g(x) <= 1+M 으로 사실상 g(x)는 비제약 상태

# @constraint(model, [i = 1:n], sum(x[i,:]) == 1)
# 각 행의 모든 원소의 합이 1이다
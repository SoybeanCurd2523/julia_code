using JuMP
import HiGHS
import MultiObjectiveAlgorithms as MOA
import Plots

profit = [77, 94, 71, 63, 96, 82, 85, 75, 72, 91, 99, 63, 84, 87, 79, 94, 90]
desire = [65, 90, 90, 77, 95, 84, 70, 94, 66, 92, 74, 97, 60, 60, 65, 97, 93]
weight = [80, 87, 68, 72, 66, 77, 99, 85, 70, 93, 98, 72, 100, 89, 67, 86, 91]
capacity = 900
N = length(profit)

# scatter 동사, 흩뿌리다
# ; 뒤에 오는 인자는 키워드 인자로, 플롯의 스타일을 결정하는데 사용

Plots.scatter(profit, desire ; xlabel = "profit", ylabel = "desire", legend = false)

model = Model()
@variable(model, x[1:N], Bin)
# @constraint(model, sum(weight[i]*x[i] for i in 1:N) <= capacity )
@constraint(model, weight'*x <= capacity)

# JuMP에서 변수, 상수, 수학적 계산식을 간단히 정의하기 위해 사용
@expression(model, profit_expr, sum(profit[i]*x[i] for i in 1:N))
@expression(model, desire_expr, sum(desire[i]*x[i] for i in 1:N))
@objective(model, Max, [profit_expr, desire_expr])

set_optimizer(model, () -> MOA.Optimizer(HiGHS.Optimizer))
set_silent(model)
set_attribute(model, MOA.Algorithm(), MOA.EpsilonConstraint())

optimize!(model)
if !is_solved_and_feasible(model)
    error("error")
end

solution_summary(model)

# multi-objective optimization problem은 여러 개의 최적 해를 가진다.
# result_count(model) 이 9가 출력된다.

plot = Plots.scatter(
    [value(profit_expr ; result = i) for i in 1:result_count(model)],
    [value(desire_expr ; result = i) for i in 1:result_count(model)],
    xlabel = "profit",
    ylabel = "desire",
    title = "objective space",
    label = false,
    xlim = (915, 960),
)

for i in 1:result_count(model)
    y = objective_value(model ; result = i)
    # plot에 주석을 추가
    # y[1]-1 : x축 위치. profit 값보다 약간 왼쪽 에 표시
    # y[2] : y축 위치. desire 값에 표시
    # (i,10) : 주석 내용. 결과 i의 인덱스 표시
    Plots.annotate!(y[1]-1, y[2], (i, 10))
end

# 모델의 이상적 해(각 목적 함수에서 가능한 최대/최소 값)
# (955, 983)
ideal_point = objective_bound(model)
Plots.scatter!([ideal_point[1]], [ideal_point[2]] ; label = "ideal point")

items_chosen = [i for i in 1:N if value(x[i]; result = 7) > 0.9]
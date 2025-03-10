# nested : 중첩된

using JuMP
import Ipopt

# x...  은 가변 인자로, 함수가 임의의 개수의 인자를 받을수 있는 것을 나타냄
function solve_lower_level(x...)
    # println("solve_lower_level function")
    # [println("x[$i] : $i") for i in 1:length(x)]
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    @variable(model, y[1:2])
    @objective(model, Max, x[1]^2*y[1] + x[2]^2*y[2] - x[1]*y[1]^4 - 2*x[2]*y[2]^4)
    @constraint(model, (y[1] - 10)^2 + (y[2] - 10)^2 <= 25)
    optimize!(model)
    solution_summary(model)
    # println("result_count : ", result_count(model))
    # println("objective_value : ", objective_value(model))
    # println("optimal y : ", value.(y))
    return objective_value(model), value.(y)
end

function V(x...)
    # println("V function")
    # [println("x[$i] : $i") for i in 1:length(x)]

    # f,_ 는 solve_lower_level 함수가 반환하는 튜플의 두 번째 값을 무시하고 첫첫 번째 값만 사용
    f,_ = solve_lower_level(x...)
    return f
end

# ∇ \nabla : gradient 기호
# g::AbstractVector 는 함수의 첫 번째 인수가 추상벡터 타입. g의 입력 타입을 제한함
function ∇V(g::AbstractVector, x...)
    # println("∇V function")
    _,y = solve_lower_level(x...)

    # V(x1,x2)의 x1에 대한 편미분
    g[1] = 2*x[1] *y[1] - y[1]^4

    # V(x1,x2)의 x2에 대한 편미분
    g[2] = 2*x[2]*y[2] - 2*y[2]^4

    return
end
    
# ² \^2 제곱 기호  
function ∇²V(H::AbstractMatrix, x...)
    # println("∇²V function")
    _,y = solve_lower_level(x...)

    # Hessian 헤시안 행렬
    H[1,1] = 2*y[1]
    H[2,2] = 2*y[2]
    # H[1,2] = H[2,1] = 0
    return
end

model = Model(Ipopt.Optimizer)
set_silent(model)
@variable(model, x[1:2] >= 0)
@operator(model, op_V, 2, V, ∇V, ∇²V)
# @operator(model, op_V, 2, V, ∇V)

# IPOPT is designed to exploit 1st derivative (gradient) and 2nd derivative (Hessian) information if provided 
# (usually via automatic differentiation routines in modeling environments such as AMPL). 
# If no Hessians are provided, IPOPT will approximate them using a quasi-Newton methods, specifically a BFGS update.

# V가 최적화 문제를 내부적으로 풀기 때문에, 1차,2차 도함수를 계산하기 위한 자동 미분을 사용할 수 없다.
# 그래서 V의 gradient(와 hessian 함수)를 만들어 줘야 한다.
@objective(model, Min, x[1]^2 + x[2]^2 + op_V(x[1], x[2]))
optimize!(model)
solution_summary(model)
# println("result_count : ", result_count(model))
# println("objective_value : ", objective_value(model))
# println("optimal x : ", value.(x))
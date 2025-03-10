using JuMP
import Ipopt
import Test

# ctrl + L : repl 터미널 텍스트 지우기
# alt +j and alt +r : workspace 초기화
function_calls = 0
function foo(x, y)
    global function_calls += 1 # 글로벌 변수
    common_term = x^2 + y^2
    term_1 = sqrt(1 + common_term)
    term_2 = common_term
    return term_1, term_2
end

# # 함수 foo 가 두 개의 return 값을 갔고 있기에 두 개의 함수로 나눔
foo_1(a, b) = foo(a, b)[1]
foo_2(a, b) = foo(a, b)[2]

# function foo_1(x, y)
#     a, _ = foo(x, y)
#     return a
# end

# function foo_2(x, y)
#     _, b = foo(x, y)
#     return b
# end



model = Model(Ipopt.Optimizer)
@variable(model, x[1:2] >= 0, start = 0.1)

# @operator로 사용자 정의 함수를 모델에 추가. 
# op_foo_1 : JuMP에서 정의된 연산자의 이름, 2 : 연산자가 받는 인수의 개수, foo_1 : 사용자 정의 함수
@operator(model, op_foo_1, 2, foo_1)
@operator(model, op_foo_2, 2, foo_2)
@objective(model, Max, op_foo_1(x[1], x[2]))
@constraint(model, op_foo_2(x[1], x[2]) <= 2)
# function_calls = 0
optimize!(model)

Test.@test objective_value(model) ≈ √3 atol = 1e-4
Test.@test value.(x) ≈ [1.0, 1.0] atol = 1e-4

println("result_count : ", result_count(model))
println("objective_value : ", objective_value(model))
println("optimal value : ", value.(x))
println("Naive approach : function_calls = $(function_calls)")
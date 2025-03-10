using JuMP, Ipopt, LinearAlgebra, Plots
model = Model(Ipopt.Optimizer)
@variable(model, x[1:10])
# @variable(model, y)
# @variable(model, z)

# println("x = ", x) # x = VariableRef[x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10]]
# prime_number = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
sin_func = [sin(pi/10*i) for i in 1:10] # 대괄호는 결과를 벡터로 만들어 준다.
# println(sin_func)
# println("===============")
# println(x) # x[1] ~ x[10] 가능
# println(typeof(x)) # Vector{Float64}

# 다변수 함수의 최적화
# 목적함수에서 벡터연산은 직접 사용할 수 없고, 각 요소에 대한 연산을 명시적으로 작성해야 한다.
@objective(model, Min, sum( (x[i] - sin_func[i])^2 for i in 1:10))
# @objective(model, Min, x - y + 2 * z + 1)
# @objective(model, Min, sum(x_vector[i]-i for i in 1:10 ) )# Σi=1 to 10, (x-i)^2  

# # 구속조건 설정
# @constraint(model, x >= 0)
# @constraint(model, y >= 1)
# @constraint(model, z >= 2)

optimize!(model)

println("Optimal x: ",value.(x)) # value. 을 사용하여 벡터의 각 원소 출력

plot(1:10, value.(x))
# println("Optimal y: ", value(y))    
# println("Optimal z: ", value(z))
println("optimal value of cost function: ", objective_value(model))

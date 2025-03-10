using JuMP
using Ipopt
using Random
using Statistics
using Test

function exmaple_qcp()
    model = Model(Ipopt.Optimizer)
    @variable(model, x)
    @variable(model, y >= 0)
    @variable(model, z >= 0)
    @objective(model, Max, x)

    @constraint(model, x+y+z == 1)

    # 구속조건이 이차항으로 표현됨
    @constraint(model, x*x + y*y - z*z <= 0)
    @constraint(model, x*x - y*z <= 0)
    optimize!(model)
    
    println("result count : ", result_count(model))
    println("objective value : ", objective_value(model))
    println("x = ", value(x))
    println("y = ", value(y))
    println("z = ", value(z))

    @test objective_value(model) ≈ 0.32699 atol = 1e-5
    @test value(x) ≈ 0.32699 atol = 1e-5
    @test value(y) ≈ 0.25707 atol = 1e-5
    return
end

exmaple_qcp()
using JuMP
using Ipopt
using Random
using Statistics
using Test

function example_rosenbrock()
    model = Model(Ipopt.Optimizer)
    # set_silent(model)
    @variable(model, x)
    @variable(model, y)
    @objective(model, Min, (1-x)^2 + 100(y-x^2)^2)
    optimize!(model)
    is_solved_and_feasible(model)
    println("result_count : ", result_count(model))
    println("optimal objective value : ", objective_value(model))

    # Test.@test는 테스트 목적을 위해 설계된 매크로로, 실패한 조건에 대한 상세한 메시지를 제공하며 프로그램이 종료되지는 않음.
    @test objective_value(model) ≈ 0.0 atol = 1e-10
    @test value(x) ≈ 1.0
    @test value(y) ≈ 1.0
    for i in all_variables(model)
        println("optimal value of $i : ", value(i))
    end
    return
end

example_rosenbrock()



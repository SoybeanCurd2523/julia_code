using DelimitedFiles, Plots, Statistics, JuMP, Ipopt, Random

sim_data_file_path = "C:\\Users\\super\\OneDrive\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\new_th5_values.txt" # 최적화 전
human_data_file_path = "C:\\Users\\super\\OneDrive\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\subjmean.txt" # 논문의 값
optimal_sim_data_file_path = "C:\\Users\\super\\OneDrive\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\optimal_sim_data.txt" # 최적화 이후

sim_data = readdlm(sim_data_file_path) # 101×1 Matrix{Float64}, (101, 1)
human_data = readdlm(human_data_file_path)
optimal_sim_data = readdlm(optimal_sim_data_file_path)

function RMSE(data1, data2)
     return sqrt(mean((data1 .- data2).^2))
end

model = Model(Ipopt.Optimizer)

@variable(model, r1, start=14) # 초기값을 설정하지 않으면 0으로 시작
@variable(model, r2, start=3)
@variable(model, r5, start=12)
@variable(model, r6, start=21)
@variable(model, th1, start=deg2rad(13)) # radian
# @variable(model, w2, start=rand(Float64))
@variable(model, offset, start = deg2rad(10))

th2_init = deg2rad(248)

function calc_theta_sim(r1, r2, r5, r6, th1, th2)
     rs = sqrt(max(r1^2 + r2^2 + 2*r1*r2*cos(th1-th2), 0)) # sqrt 정의역 조건
     if rs ≈ 0 
          print("error! rs is 0 \n")
          println("r1 : ", r1, ", r2 : ", r2, ", r5 : ", r5, ", r6 : ", r6, ", th1(rad) : ", th1, ", w2(rad/s) : ", w2)
          return 0.0
     end

     ths = atan(r1*sin(th1) + r2*sin(th2), r1*cos(th1) + r2*cos(th2))
     input_to_acos = (-r6^2+r5^2+rs^2)/ (2*r5*rs) # acos 정의역 : -1과 1 사이여야 함
     if input_to_acos < -1
          print("error! input_to_acos < -1 \n")
          println("r1 : ", r1, ", r2 : ", r2, ", r5 : ", r5, ", r6 : ", r6, ", th1(rad) : ", th1, ", w2(rad/s) : ", w2)
          input_to_acos = -1
     elseif input_to_acos > 1
          print("error! input_to_acos > 1 \n")
          println("r1 : ", r1, ", r2 : ", r2, ", r5 : ", r5, ", r6 : ", r6, ", th1(rad) : ", th1, ", w2(rad/s) : ", w2)
          input_to_acos = 1
     end
     # valid_to_acos = clamp(input_to_acos, -1, 1)
     theta_sim = π/2 - (ths - acos(input_to_acos) + π)
     return theta_sim
end

# @NLconstraint(model, (-r6^2+r5^2+sqrt(r1^2 + r2^2 + 2*r1*r2*cos(th1-th2))^2) / (2*r5*sqrt(r1^2 + r2^2 + 2*r1*r2*cos(th1-th2))) >= -1)
# @NLconstraint(model, (-r6^2+r5^2+sqrt(r1^2 + r2^2 + 2*r1*r2*cos(th1-th2))^2) / (2*r5*sqrt(r1^2 + r2^2 + 2*r1*r2*cos(th1-th2))) <= 1)

#=
 grashof's law : min(r1,r2,r5,r6)+max(r1,r2,r5,r6) <= sum([r1,r2,r5,r6])-min(r1,r2,r5,r6)-max(r1,r2,r5,r6)
 @constraint(model, sort([value(r1),value(r2),value(r5),value(r6)])[1] + sort([value(r1),value(r2),value(r5),value(r6)])[4] 
                     <= sort([value(r1),value(r2),value(r5),value(r6)])[2] + sort([value(r1),value(r2),value(r5),value(r6)])[3])
=#

@constraint(model, r2 + r6 <= r1 + r5)

@constraint(model, r2 <= r6)
@constraint(model, r2 <= r1)
@constraint(model, r2 <= r5)

# r2 <= r6은 위에 이미 있음
@constraint(model, r1 <= r6)
@constraint(model, r5 <= r6)


@constraint(model, r1 <= 14*1.2)
@constraint(model, r1 >= 14*0.8)
@constraint(model, r2 <= 6)
@constraint(model, r2 >= 3*0.8)
@constraint(model, r5 <= 50)
@constraint(model, r5 >= 12*0.8)
@constraint(model, r6 <= 50*sqrt(2)) # 70.71
@constraint(model, r6 >= 21*0.8)
@constraint(model, th1 <= deg2rad(90))
@constraint(model, th1 >= 0)
# @constraint(model, w2 <= deg2rad(180))
# @constraint(model, w2 >= deg2rad(36))
# @constraint(model, offset >= deg2rad(10))
# @constraint(model, offset <= deg2rad(20))

@NLobjective(model, Min, sum( (deg2rad.(human_data[i]) - offset - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101))
optimize!(model)

println("Optimal r1: ", value(r1)) # value(x)는 JuMP에서 최적화 후 변수를 평가하는 표준 방법
println("Optimal r2: ", value(r2))
println("Optimal r5: ", value(r5))
println("Optimal r6: ", value(r6))
print("Optimal th1(radian): ", value(th1))
println(", th1(degree): ", rad2deg(value(th1)))
# print("Optimal w2(rad/s): ", value(w2))
# println(", w2(deg/s): ", rad2deg(value(w2)))
println("offset(degree) : ", rad2deg(value(offset)))
println("==============================")
println("Optimal value of cost function: ", objective_value(model))
println("correlation coefficient : ", cor(human_data, optimal_sim_data));
println("RMSE : ", RMSE(human_data, optimal_sim_data));

######################################

plot(human_data, label="human_data", xlabel="% gait cycle", ylabel="hip angle(degree)", linewidth=2, title="optimize!")
plot!(sim_data, label="sim_data", color="red", linewidth=2)
# plot!(optimal_sim_data, label="optimal_sim_data", color="green", linewidth=2, line=:dash)
# plot!()
output = [value(offset) + calc_theta_sim(value(r1), value(r2), value(r5), value(r6), value(th1), value(th2_init) + (i-1)*2*pi/100) for i in 1:101]
plot!(rad2deg.(output), label="output", color="green", linewidth=2, line=:dash)

# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999897919621
# Optimal r2: 2.921607944324002
# Optimal r5: 50.000000455715806
# Optimal r6: 55.64809060474908
# Optimal th1(radian): -9.95459981922385e-9, th1(degree): -5.703565563832188e-7
# offset(degree) : -2.5709380969743596
# ==============================
# Optimal value of cost function: 0.24526076013157022
# correlation coefficient : [-0.02265958547583218;;]
# RMSE : 15.758127268834084
using DelimitedFiles, Plots, Statistics, JuMP, Ipopt, Random

sim_data_file_path = "C:\\Users\\Jehyeon\\OneDrive - GIST\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\optimal_sim_data2.txt" # 최적화 전. matlab에서
print(isfile(sim_data_file_path))
human_data_file_path = "C:\\Users\\Jehyeon\\OneDrive - GIST\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\subjmean.txt" # 논문의 값

sim_data = readdlm(sim_data_file_path) # 101×1 Matrix{Float64}, (101, 1)
human_data = readdlm(human_data_file_path)

function RMSE(data1, data2)
     return sqrt(mean((data1 .- data2).^2))
end

model = Model(Ipopt.Optimizer)

@variable(model, r1, start=14) # 초기값을 설정하지 않으면 0으로 시작
@variable(model, r2, start=3)
@variable(model, r5, start=12) 
@variable(model, r6, start=21)
@variable(model, th1, start=deg2rad(13)) # radian

th2_init = deg2rad(248) # 248도??

function calc_theta_sim(r1, r2, r5, r6, th1, th2)

     #rs 계산 시 음수 내부 제곱근 방지
     rs = sqrt(max(r1^2 + r2^2 + 2*r1*r2*cos(th1-th2), 0.0)) # sqrt 정의역 조건2

     # MATLAB과 동일하게 두 좌표의 사잇각은 atan2를 사용하여 계산
     ths = atan(r1*sin(th1) + r2*sin(th2), r1*cos(th1) + r2*cos(th2))

     # acos 입력값 계산 및 클램핑 (정의역 [-1, 1])
     val = (-r6^2+r5^2+rs^2)/ (2*r5*rs)
     val = clamp(val, -1.0, 1.0)

     # 두 후보값 계산: MATLAB과 동일하게 min()을 사용
     th5_1 = ths + acos(val)
     th5_2 = ths - acos(val)

     # MATLAB에서는 두 값 중 min()을 취한 후 pi를 더함
     th5 = min(th5_1, th5_2) + π

     # 수직(90°) 기준 변환: MATLAB은 new_th5 = 90 - rad2deg(th5)
     theta_sim = deg2rad(90) - th5
     return theta_sim
end

# grashof's condition
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
@constraint(model, r5 <= 43) #상한이 성인남성의 평균 thigh 길이인 43cm
@constraint(model, r5 >= 12*0.8)
@constraint(model, r6 <= 43*sqrt(2)) 
@constraint(model, r6 >= 21*0.8)
@constraint(model, th1 <= deg2rad(90))
@constraint(model, th1 >= 0)


@NLobjective(model, Min, sum( (deg2rad.(human_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101))
optimize!(model)


# 최적화 결과
# rad임
output = [calc_theta_sim(value(r1), value(r2), value(r5), value(r6), value(th1), value(th2_init) + (i-1)*2*pi/100) for i in 1:101]


println("Optimal r1: ", value(r1)) # value(x)는 JuMP에서 최적화 후 변수를 평가하는 표준 방법
println("Optimal r2: ", value(r2))
println("Optimal r5: ", value(r5))
println("Optimal r6: ", value(r6))
print("Optimal th1(radian): ", value(th1))
println(", th1(degree): ", rad2deg(value(th1)))
println("==============================")
println("Optimal value of cost function: ", objective_value(model))
println("correlation coefficient : ", cor(human_data, rad2deg.(output)));
println("RMSE : ", RMSE(human_data, rad2deg.(output)));

######################################

plot(human_data, label="human_data", color="red", xlabel="% gait cycle", ylabel="hip angle(degree)", linewidth=2, title="optimize!!")
plot!(sim_data, label="sim_data", color="blue", linewidth=2)
plot!(rad2deg.(output), label="output", color="green", linewidth=2, line=:dash)

# plot(rad2deg.(output), label="output", color="green", linewidth=2, line=:dash, xlabel="% gait cycle", ylabel="hip angle(degree)")

writedlm("C:\\Users\\Jehyeon\\OneDrive - GIST\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\optimal_theta5_output.txt", output)


# new_julia_angle =  rad2deg.( [calc_theta_sim(14.0, 3.0, 12.0, 21.0, deg2rad(13), deg2rad(248) + (i-1)*2*pi/100) for i in 1:101])
# sim_data .- rad2deg.( [calc_theta_sim(14.0, 3.0, 12.0, 21.0, deg2rad(13), deg2rad(248) + (i-1)*2*pi/100) for i in 1:101])
# 차이가 0.28% 이긴 한데...... 좀 그렇긴 하다.. offset 빼도 차이가 남

# EXIT: Optimal Solution Found.
# Optimal r1: 11.422747791577113
# Optimal r2: 2.9934522256178693
# Optimal r5: 42.27568890339048
# Optimal r6: 47.81152056826159
# Optimal th1(radian): -9.953306729289433e-9, th1(degree): -5.70282467787446e-7
# ==============================
# Optimal value of cost function: 0.24788896401983262
# correlation coefficient : [0.9788393679174849;;]
# RMSE : 2.8385107397410234
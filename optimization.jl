using DelimitedFiles, Plots, Statistics, JuMP, Ipopt, Random

subject_number = 11

mean_cycle_file_path ="C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\result\\subject$(subject_number)\\mean_cycle.txt"
sim_data_file_path = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\optimal_sim_data2.txt" # 최적화 전. matlab에서
human_data_file_path = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\subjmean.txt" # 논문의 값

mean_data = readdlm(mean_cycle_file_path)
sim_data = readdlm(sim_data_file_path) # 101×1 Matrix{Float64}, (101, 1)
human_data = readdlm(human_data_file_path)

function RMSE(data1, data2)
     return sqrt(mean((data1 .- data2).^2))
end

function MSE(data1, data2)
     return mean((data1 .- data2).^2)
end

function MAE(data1, data2)
     return mean(abs.(data1 .- data2))
end

model = Model(Ipopt.Optimizer)

@variable(model, r1, start=14) # 초기값을 설정하지 않으면 0으로 시작
@variable(model, r2, start=3)
@variable(model, r5, start=12) 
@variable(model, r6, start=21)
# @variable(model, th1, start=deg2rad(13)) # radian

th2_init = deg2rad(248) # 248도??
th1_init = deg2rad(13)

function calc_theta_sim(r1, r2, r5, r6, th1, th2) # θ₂ 에서 θ₅ 를 구하는 함수

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
# @constraint(model, th1 <= deg2rad(90))
# @constraint(model, th1 >= 0)

# SSE
@NLobjective(model, Min, sum((deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1_init, th2_init + (i-1)*2*pi/100))^2 for i in 1:101))

# MSE
# @NLobjective(model, Min, (1/101)*sum((deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101))

# RMSE
# @NLobjective(model, Min, sqrt( (1/101)*sum((deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101) ))

# MAE
# @NLobjective(model, Min, (1/101)*sum(abs(deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100)) for i in 1:101))


optimize!(model)


# 최적화 결과
# rad임
output = [calc_theta_sim(value(r1), value(r2), value(r5), value(r6), th1_init, th2_init + (i-1)*2*pi/100) for i in 1:101]

println("Given th1(deg) : ", rad2deg(th1_init))
println("Given th2(deg) : ", rad2deg(th2_init))
println("==============================")
println("Optimal r1: ", value(r1)) # value(x)는 JuMP에서 최적화 후 변수를 평가하는 표준 방법
println("Optimal r2: ", value(r2))
println("Optimal r5: ", value(r5))
println("Optimal r6: ", value(r6))
println("==============================")
println("Optimal value of cost function: ", objective_value(model))
println("Pearson correlation coefficient : ", cor(mean_data, rad2deg.(output)));
println("RMSE : ", RMSE(mean_data, rad2deg.(output)));

######################################

p2 = plot(mean_data, label="mean_data", color="red", xlabel="% gait cycle", ylabel="thigh angle(degree)", linewidth=2, title="optimization")
plot!(sim_data, label="sim_data", color="blue", linewidth=2)
plot!(rad2deg.(output), label="output", color="green", linewidth=2, line=:dash)

display(p2)
# plot(rad2deg.(output), label="output", color="green", linewidth=2, line=:dash, xlabel="% gait cycle", ylabel="hip angle(degree)")

# writedlm("C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\data\\optimal_theta5_output.txt", output)

# 1) 최적화된 링크 길이 & 각도 꺼내서 벡터로
optimal_r = [
    value(r1),
    value(r2),
    value(r5),
    value(r6),
]

# 2) 파일에 저장 (한 줄로 콤마 구분)
writedlm("result/subject$(subject_number)/optimal_r_values.txt", optimal_r, ',')
println("Saved optimal r to optimal_r_values.txt")

# hugadb 데이터셋 left thigh 정상 : 
# subject 1, 3, 5, 6, 7, 11
# 01_00, 03_00, 05_00, ** 06_00, 07_00, 11_00


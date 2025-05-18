using DelimitedFiles, Plots, Statistics, JuMP, Ipopt, Random

subject_number = 18
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
@variable(model, th1, start=deg2rad(13)) # radian

th2_init = deg2rad(248) # 248도??

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
@constraint(model, th1 <= deg2rad(90))
@constraint(model, th1 >= 0)

# SSE
@NLobjective(model, Min, sum((deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101))

# MSE
# @NLobjective(model, Min, (1/101)*sum((deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101))

# RMSE
# @NLobjective(model, Min, sqrt( (1/101)*sum((deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100))^2 for i in 1:101) ))

# MAE
# @NLobjective(model, Min, (1/101)*sum(abs(deg2rad.(mean_data[i]) - calc_theta_sim(r1, r2, r5, r6, th1, th2_init + (i-1)*2*pi/100)) for i in 1:101))


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
println("Pearson correlation coefficient : ", cor(mean_data, rad2deg.(output)));
println("RMSE : ", RMSE(mean_data, rad2deg.(output)));
println("MSE : ", MSE(mean_data, rad2deg.(output)));
println("MAE : ", MAE(mean_data, rad2deg.(output)));

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
    value(th1), # radian
]

# 2) 파일에 저장 (한 줄로 콤마 구분)
writedlm("result/subject$(subject_number)/optimal_r_values.txt", optimal_r, ',')
println("Saved optimal r&th1 to optimal_r_values.txt")

# 1, 3, 5, 6, 7, 8, 11
# subject 1번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888249094
# Optimal r2: 2.3999999760812027
# Optimal r5: 43.00000042942297
# Optimal r6: 51.2090098900802
# Optimal th1(radian): 0.43118300471720133, th1(degree): 24.704966368065104
# ==============================
# Optimal value of cost function: 1.4451617911743495
# Pearson correlation coefficient : [0.9175245471764104;;]
# RMSE : 6.853622518952445
# MSE : 46.972141632292065
# MAE : 5.9850450718482495

# subject 3번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888031206
# Optimal r2: 6.000000059971762
# Optimal r5: 43.00000042970637
# Optimal r6: 43.000000420177535
# Optimal th1(radian): 0.011849103029943454, th1(degree): 0.6789035946314358
# ==============================
# Optimal value of cost function: 2.2183561624331993
# Pearson correlation coefficient : [0.953706092883805;;]
# RMSE : 8.491366434618415
# MSE : 72.10330392696426
# MAE : 7.106196564389459

# subject 4번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888004385
# Optimal r2: 6.000000059995855
# Optimal r5: 43.000000429961155
# Optimal r6: 48.19948757791987
# Optimal th1(radian): -9.997668644492766e-9, th1(degree): -5.728242182997142e-7
# ==============================
# Optimal value of cost function: 140.79662391573504
# Pearson correlation coefficient : [0.7949466571782945;;]
# RMSE : 67.64848362776907
# MSE : 4576.317337136539   <<< 값이 매우 이상함.
# MAE : 56.624122672885655

# subject 5번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888038406
# Optimal r2: 6.00000005996432
# Optimal r5: 43.00000042965175
# Optimal r6: 46.75200214152945
# Optimal th1(radian): 0.4102976073179206, th1(degree): 23.50832124363281
# ==============================
# Optimal value of cost function: 3.369037630294835
# Pearson correlation coefficient : [0.9464003003344161;;]
# RMSE : 10.464413115085444
# MSE : 109.50394184317223
# MAE : 8.799902800668386

# subject 6번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888080551
# Optimal r2: 2.7476503174405607
# Optimal r5: 43.00000042969074
# Optimal r6: 51.32223152477398
# Optimal th1(radian): 0.6617309678267027, th1(degree): 37.91439162957733
# ==============================
# Optimal value of cost function: 4.37746634795146
# Pearson correlation coefficient : [0.8683298335688551;;]
# RMSE : 11.928157187980299
# MSE : 142.28093390116607
# MAE : 10.480136221456252

# subject 7번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888167754
# Optimal r2: 6.000000059910131
# Optimal r5: 28.32669489204983
# Optimal r6: 29.9501163536105
# Optimal th1(radian): 0.18442207748084505, th1(degree): 10.566606688687083
# ==============================
# Optimal value of cost function: 1.3633734592590832
# Pearson correlation coefficient : [0.9710086881056553;;]
# RMSE : 6.656859076716242
# MSE : 44.31377276725942
# MAE : 5.711737438763244

# subject 8번 2번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888135592
# Optimal r2: 2.502042642260585
# Optimal r5: 43.00000042947943
# Optimal r6: 51.32930169419344
# Optimal th1(radian): 0.4084035525809562, th1(degree): 23.399799901037987
# ==============================
# Optimal value of cost function: 2.2456034596943235
# Pearson correlation coefficient : [0.9036210209439969;;]
# RMSE : 8.543355531264792
# MSE : 72.98892373359271
# MAE : 7.3462703230599775

# subject 9번 1번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.19999988800253
# Optimal r2: 6.000000059998645
# Optimal r5: 19.076908248624285
# Optimal r6: 24.276908086627852
# Optimal th1(radian): 0.40833481155574375, th1(degree): 23.395861330413915
# ==============================
# Optimal value of cost function: 156.34472426443628
# Pearson correlation coefficient : [0.9269787267857365;;]
# RMSE : 71.28588751272233
# MSE : 5081.677758476501 << 값이 매우 이상함
# MAE : 62.14780392064945

# subject 10번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 16.800000167975433
# Optimal r2: 2.3999999760024955
# Optimal r5: 43.000000429845315
# Optimal r6: 57.331370939135084
# Optimal th1(radian): -9.999456273659294e-9, th1(degree): -5.729266419062907e-7
# ==============================
# Optimal value of cost function: 14.630771595141374
# Pearson correlation coefficient : [-0.30602255469684037;;]  << 값이 이상함
# RMSE : 21.80698178191216
# MSE : 475.5444544366489
# MAE : 19.12463986457209

# subject 11번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888048184
# Optimal r2: 4.042114898387318
# Optimal r5: 43.00000042981501
# Optimal r6: 50.04485324756911
# Optimal th1(radian): 0.6346285521925285, th1(degree): 36.36153759912977
# ==============================
# Optimal value of cost function: 5.451964801111082
# Pearson correlation coefficient : [0.8998353138570311;;]
# RMSE : 13.31185159035524
# MSE : 177.2053927636433
# MAE : 11.89451084332695

# subject 12번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888002386
# Optimal r2: 6.000000059998637
# Optimal r5: 43.00000042985239
# Optimal r6: 48.200000267855685
# Optimal th1(radian): 0.2360438994784012, th1(degree): 13.524319219922644
# ==============================
# Optimal value of cost function: 138.68833368174617
# Pearson correlation coefficient : [0.8984683011870114;;]
# RMSE : 67.14008869321304
# MSE : 4507.791509732515 << 값이 매우 이상함
# MAE : 58.358919360054685

# subject 13번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888006076
# Optimal r2: 6.000000059996864
# Optimal r5: 11.599999994229869
# Optimal r6: 16.79999983223855
# Optimal th1(radian): 0.42148623660569134, th1(degree): 24.14938248035854
# ==============================
# Optimal value of cost function: 73.6204048924139
# Pearson correlation coefficient : [0.9205605064451304;;]
# RMSE : 48.917138735739684
# MSE : 2392.886462091604 << 값이 매우 이상함
# MAE : 42.977417667506664

# subject 14번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Invalid number in NLP function or derivative detected. << 최적화 안됨
# Optimal r1: 11.199999888258601
# Optimal r2: 6.000000059855594
# Optimal r5: 43.000000408264256
# Optimal r6: 48.20000024666638
# Optimal th1(radian): 0.36997740917118604, th1(degree): 21.198144060693718
# ==============================
# Optimal value of cost function: 0.0
# Pearson correlation coefficient : [0.9211142668201543;;]
# RMSE : 57.37739923113067
# MSE : 3292.1659425285543 << 값이 매우 이상함
# MAE : 49.41587564389811

# subject 15번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Invalid number in NLP function or derivative detected.  << 최적화 안됨
# Optimal r1: 11.200000015614659
# Optimal r2: 5.999999991119727
# Optimal r5: 42.999988596671116
# Optimal r6: 48.19998863098277
# Optimal th1(radian): 0.3071389845772189, th1(degree): 17.597767540208327
# ==============================
# Optimal value of cost function: 0.0
# Pearson correlation coefficient : [0.9084510149480618;;]
# RMSE : 64.1575996668838
# MSE : 4116.197595016129  << 값이 매우 이상함
# MAE : 55.29105782263635

# subject 16번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888003232
# Optimal r2: 6.000000059998269
# Optimal r5: 12.228071272832086
# Optimal r6: 17.428071110836978
# Optimal th1(radian): 0.4341307586464453, th1(degree): 24.873860227253886
# ==============================
# Optimal value of cost function: 152.4884834007173
# Pearson correlation coefficient : [0.9173353993559658;;]
# RMSE : 70.40126585351697
# MSE : 4956.338233777576 << 값이 매우 이상함
# MAE : 61.43575455563953

# subject 17번 0번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888003193
# Optimal r2: 6.000000059997642
# Optimal r5: 43.000000429955406
# Optimal r6: 48.20000026795844
# Optimal th1(radian): 0.0690182422545107, th1(degree): 3.954453990594947
# ==============================
# Optimal value of cost function: 126.78969218604736
# Pearson correlation coefficient : [0.8342108181410748;;]
# RMSE : 64.19540147501628
# MSE : 4121.049570538523 << 값이 매우 이상함
# MAE : 54.10884635790832

# subject 18번 3번째 파일
# imu데이터 LPF cutoff 3Hz 일 때
# EXIT: Optimal Solution Found.
# Optimal r1: 11.199999888002674
# Optimal r2: 6.0000000599982295
# Optimal r5: 43.00000042994623
# Optimal r6: 48.200000267947146
# Optimal th1(radian): 0.11682011725170242, th1(degree): 6.693299680745967
# ==============================
# Optimal value of cost function: 139.28737965373654
# Pearson correlation coefficient : [0.8525621213250256;;]
# RMSE : 67.28493383239322
# MSE : 4527.2623208295345 << 값이 매우 이상함
# MAE : 57.60334452384531
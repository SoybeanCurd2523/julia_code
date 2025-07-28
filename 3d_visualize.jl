# 1. 패키지 불러오기 및 시드 설정
using Plots
using Random
using Statistics # std() 함수 사용을 위해 추가
using Printf   # 결과 출력을 위해 추가

plotly() # 동적 그래프를 위한 백엔드 설정
Random.seed!(42) # 재현 가능성을 위해 랜덤 시드 고정

# 2. 원본 데이터 정의
subjects_orig = ["S1", "S3", "S5", "S6", "S7", "S11"]
heights_orig = [177, 183, 183, 168, 161, 203]
weights_orig = [75, 65, 63, 54, 52, 85]
ages_orig = [24, 23, 35, 25, 22, 24]

# 원본 r값 데이터 추가
r_values_orig = [
    14.02 4.50 16.54 20.21; # S1
    12.35 5.50 36.37 41.81; # S3
    11.20 5.52 43.00 47.57; # S5
    11.20 4.65 22.15 25.14; # S6
    11.20 4.65 43.00 48.53; # S7
    11.20 5.30 43.00 45.65; # S11
]

# 3. 노이즈 데이터 생성을 위한 설정 및 배열 초기화
num_noise_points = 100 # 각 subject당 생성할 노이즈 데이터 개수
heights_noise = Float64[] # 노이즈 키 데이터
weights_noise = Float64[] # 노이즈 몸무게 데이터
ages_noise = Int[]      # 노이즈 나이 데이터

r_values_noise = zeros(num_noise_points * length(subjects_orig), 4)


# 노이즈 수준 (표준편차)
height_std = 5.0
weight_std = 3.0
age_std = 2.0

r_noise_factor = 0.1 # r값 표준편차의 10% 수준으로 노이즈 설정
r_stds = [std(r_values_orig[:, i]) for i in 1:4] # 각 r 열의(r1끼리, r2끼리 ...) 표준편차 계산


# 4. 데이터 증강 루프
noise_data_index = 1 # 노이즈 데이터 배열의 인덱스
for i in 1:length(subjects_orig) # 6명의 원본 데이터에 대해 반복
    for _ in 1:num_noise_points # 각 1명당 100번씩 반복
        # 원본 데이터에 정규분포 노이즈 추가
        new_h = heights_orig[i] + randn() * height_std
        new_w = weights_orig[i] + randn() * weight_std
        new_a = round(Int, ages_orig[i] + randn() * age_std)
        
        # 생성된 데이터를 배열에 추가
        push!(heights_noise, new_h)
        push!(weights_noise, new_w)
        push!(ages_noise, new_a)

        # r값에 대한 적응형 노이즈 추가 (Y값)
        for j in 1:4
            noise = randn() * r_stds[j] * r_noise_factor
            r_values_noise[noise_data_index, j] = r_values_orig[i, j] + noise
        end
        noise_data_index += 1 # 행 인덱스 인 듯
    end
end



# 3. 예측 실행 파트 (새로운 부분)
# 역거리 가중법(Inverse Distance Weighting, IDW) 알고리즘

# 3-1. 전체 데이터셋 통합
X_full = vcat(hcat(heights_orig, weights_orig, ages_orig), hcat(heights_noise, weights_noise, ages_noise))
Y_full = vcat(r_values_orig, r_values_noise)

# 3-2. 예측할 새로운 사람 정보 정의
new_person = (height=180, weight=80, age=23)


# 3-3. 정규화 준비
all_heights = vcat(X_full[:, 1], new_person.height)
all_weights = vcat(X_full[:, 2], new_person.weight)
all_ages = vcat(X_full[:, 3], new_person.age)
min_h, max_h = minimum(all_heights), maximum(all_heights)
min_w, max_w = minimum(all_weights), maximum(all_weights)
min_a, max_a = minimum(all_ages), maximum(all_ages)
normalize(val, min_val, max_val) = (val - min_val) / (max_val - min_val)

# 3-4. 새로운 사람 정보 정규화

# 정규화 : 공정한 비교를 위한 '단위 통일' 과정입니다. 
# 예를 들어, 두 사람의 키 차이가 10cm이고 나이 차이가 10살이라고 할 때, 어떤 차이가 더 "큰" 차이일까요? 
# 단위와 범위가 달라 직접 비교하기 어렵습니다.

# 정규화는 모든 값을 "전체 범위에서 몇 % 지점에 위치하는가?"로 바꿔줍니다. 
# 예를 들어, 키 180cm는 "가장 작은 키(0)와 가장 큰 키(1) 사이에서 0.7 지점"으로, 
# 나이 25세는 "가장 어린 나이(0)와 가장 많은 나이(1) 사이에서 0.3 지점" 등으로 변환됩니다. 
# 이렇게 단위를 통일해야 키, 몸무게, 나이라는 세 가지 특성을 모두 공평하게 반영하여 거리를 계산할 수 있습니다.



# 예측하려는 새로운 사람의 키, 몸무게, 나이 값을 0과 1 사이의 값으로 변환
norm_h_new = normalize(new_person.height, min_h, max_h)
norm_w_new = normalize(new_person.weight, min_w, max_w)
norm_a_new = normalize(new_person.age, min_a, max_a)

# 3-5. 거리 및 가중치 계산
weights = Float64[]
for i in 1:size(X_full, 1)
    norm_h_orig = normalize(X_full[i, 1], min_h, max_h)

#     julia> normalize(177.0, 150.8908617649663, 215.9334957485121)
#            0.40141575818775543


    norm_w_orig = normalize(X_full[i, 2], min_w, max_w)
    norm_a_orig = normalize(X_full[i, 3], min_a, max_a)
    dist = sqrt((norm_h_new - norm_h_orig)^2 + (norm_w_new - norm_w_orig)^2 + (norm_a_new - norm_a_orig)^2)
    push!(weights, 1 / (dist + 1e-9)) # "가까울수록 더 중요하다". 가중치는 거리의 역수
end

# 3-6. 가중 평균으로 최종 r값 예측
weighted_r_sums = [0.0, 0.0, 0.0, 0.0]
total_weight = sum(weights)
for i in 1:size(Y_full, 1)
    for j in 1:4
        weighted_r_sums[j] += Y_full[i, j] * weights[i]
    end
end
predicted_r_values = weighted_r_sums ./ total_weight

# 4. 최종 결과 출력
println("\n--- 새로운 사람 정보 ---")
println("키: $(new_person.height), 몸무게: $(new_person.weight), 나이: $(new_person.age)")
println("\n--- 최종 예측된 r 값 ---")
@printf "r1: %.2f\n" predicted_r_values[1]
@printf "r2: %.2f\n" predicted_r_values[2]
@printf "r5: %.2f\n" predicted_r_values[3]
@printf "r6: %.2f\n" predicted_r_values[4]

# 5. 시각화
# 먼저 원본 데이터를 크고 뚜렷하게 플로팅
p = scatter3d(
    heights_orig, 
    weights_orig, 
    ages_orig,
    label = "원본 데이터",
    marker = (:circle, 8, 0.9, :red), # 모양: 원, 크기: 8, 투명도: 0.9, 색: 빨강
    series_annotations = text.(subjects_orig, :bottom, 8) # 원본 데이터에만 라벨 추가
)

# 생성된 노이즈 데이터를 작고 반투명하게 동일한 플롯에 추가
scatter3d!(
    p, # p 플롯에 이어서 그림
    heights_noise, 
    weights_noise, 
    ages_noise,
    title = "원본, 증강된 노이즈 데이터, 새로운 사람",
    xlabel = "키 (cm)",
    ylabel = "몸무게 (kg)",
    zlabel = "나이(정수)",
    label = "노이즈 데이터",
    marker = (:circle, 3, 0.4, :blue) # 모양: 원, 크기: 3, 투명도: 0.4, 색: 파랑
)

# 4-2. 새로운 사람 데이터 플롯 추가
scatter3d!(
    p,
    [new_person.height], # 배열 형태로 전달
    [new_person.weight],
    [new_person.age],
    label = "새로운 사람",
    marker = (:octagon, 8, 1.0, :green) # 별 모양, 크기, 투명도, 색상
)
# 최종 플롯을 표시합니다.
display(p)

# ##########################################################
# # ##########################################################
# # cos 함수로 calc_theta_sim을 근사화

# max_val, max_idx = findmax(data)
# min_val, min_idx = findmin(data)

# const A = (max_val - min_val) / 2
# const C = (max_val + min_val) / 2
# const ω = π / 50.0
# const t_peak = Float64(max_idx)

# # 2. 코사인 함수 데이터 생성 (수식 수정)
# function cos_data_func(t)
#     # 수직 이동 C를 밖에서 더해주는 방식으로 수정
#     y = A * cos(ω*(t - t_peak)) + C
#     return y
# end


# # julia> mean(data .- [cos_data_func(t) for t in 1:101] )
# # -0.04445997171613357

# # julia> rad2deg(-0.044)
# # -2.521014298575622


# # 후보 t를 찾는 내부 헬퍼 함수
# function inverse_cos_data_func(y::Float64)
#     if !(C - A <= y <= C + A)
#         error("y값($y)이 범위 [$(C-A), $(C+A)]를 벗어났습니다.") # 입력값 범위 초과 에러
#     end
#     theta = acos((y - C) / A)
#     candidates = Float64[]
#     for n in -2:2
#         t1 = t_peak + (1/ω) * (theta + 2n*π)
#         t2 = t_peak + (1/ω) * (-theta + 2n*π) 
# # 가능한 t가 세 개 일 때는 어떻하냐? << 그래프 y값이 시작과 끝이 같으므로, 해가 3개가 되는 부분은 그 때 뿐이다.
# # 즉, t가 101일 때나 1일 때나 같은 상황이다
#         if 1.0 <= t1 <= 101.0; push!(candidates, t1); end
#         # println("aaa  t1 : ", t1)
#         if 1.0 <= t2 <= 101.0; push!(candidates, t2); end
#         # println("bbb  t2 : ", t2)
#     end
    
#     return sort(unique(candidates))
# end

# # 최종 t를 결정하는 메인 역함수
# function _decide_t(current_y::Float64, previous_y::Float64)
#     candidates = inverse_cos_data_func(current_y)
#     println("t candidates : ", candidates)
#     if length(candidates) < 2
#         return isempty(candidates) ? nothing : first(candidates)
#     end
#     t_small, t_large = candidates[1], candidates[2]

#     # 이전 값과 비교하여 최종 t 결정
#     if current_y < previous_y
#         return t_small # 감소 중
#     else
#         return t_large # 증가 중 (또는 변화 없음)
#     end
# end

# ### 예제 실행 ###
# println("--- 모델 파라미터 ---")
# @printf "진폭(A): %.4f\n" A
# @printf "수직이동(C): %.4f\n" C
# @printf "최대값위치(t_peak): %.1f\n" t_peak
# println("--------------------")


# current_y = data[101]
# previous_y = data[100]

# println("\n--- 역함수 계산 예제 ---")
# @printf "현재 y: %.5f (rad), y: %.5f (deg), index : %.5f \n" current_y rad2deg(current_y) 101
# @printf "이전 y: %.5f (rad), y: %.5f (deg), index : %.5f \n" previous_y rad2deg(previous_y) 100

# final_t = _decide_t(current_y, previous_y)
# @printf "계산된 최종 t값: %.5f\n" final_t
# println("----------------------")

# plot!([rad2deg.(cos_data_func(t) for t in 1:101)], label="cos_data", color="green", linewidth=2)

# function myfunc()
#     i=1
#     while(i<101)
#         if i==1
#             i = 2
#         end
#         current_y = data[i]
#         previous_y = data[i-1]

#         @printf "현재 y: %.5f (rad), y: %.5f (deg), index : %.5f \n" current_y rad2deg(current_y) i
#         @printf "이전 y: %.5f (rad), y: %.5f (deg), index : %.5f \n" previous_y rad2deg(previous_y) i-1

#         final_t = _decide_t(current_y, previous_y)
#         @printf "계산된 최종 t값: %.5f\n" final_t
#         println("----------------------")

#         if i==2
#             i+=9
#         else
#             i+=10
#         end
#     end
# end

# # myfunc()

######################################################
# 이분법

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

data = [ calc_theta_sim(predicted_r_values[1], predicted_r_values[2], predicted_r_values[3], predicted_r_values[4], deg2rad(13), deg2rad(248)+(i-1)*2*pi/100) for i in 1:101 ]
# plot( rad2deg.( data ), 
#     color="red", xlabel="% gait cycle", ylabel="thigh angle(degree)", linewidth=2, label="data" )


# # pseudo code
# function bisection_method(data_array{Float64}, y_level::Float64, y_level_prev::Float64)
#     # 만약 y_level이 data_array 중 정확히 존재하는 값이라면
#         # 해당하는 인덱스(t값)을 반환

#     # 그렇지 않다면
#         # y_level이 포함되어 있을 수 있는 구간 candidates를 구하기
#             # 구간이 두 개 라면
#                 # y_level_prev를 고려해, 구간 하나로 특정
#             # bisection 방법을 통해 구간 안에서 정확한 t값을 계산하기(epsilon)
#                 # t값을 반환
# end

function bisection_method(target_function::Function, data::Vector{Float64}, y::Float64, y_prev::Float64)
    max_idx = findmax(data)[2]
        println("max_idx : ", max_idx)
        println("y>y_prev? ", y > y_prev)
        # println("candidates: ", candidates)

    tol = 1e-5      # 허용 오차
    max_iter = 100  # 최대 반복

    # 1) 샘플 데이터에서 거의 같은 값이 있으면 그 인덱스 반환
    for i in 1:length(data)
        if abs(data[i] - y) < tol
            println("find exact index")
            return float(i)
        end
    end

    # 2) y 레벨을 끼는 구간(부호 변화) 찾기
    candidates = Int[]
    for i in 1:length(data)-1
        if (data[i] - y) * (data[i+1] - y) <= 0
            push!(candidates, i)  # 해는 [i, i+1] 사이
        end
    end
    if length(candidates) == 0
        error("y 레벨을 끼는 구간을 찾지 못했습니다.")
    end

    # 3) 여러 구간이면 간단 규칙으로 선택
    # sort(candidates)

    if(length(candidates) == 1)
        idx = candidates[1]
    else #if(length(candidates) >= 2)
        for i in 1:length(candidates)
            println("candidates[$(i)] : $(candidates[i])")
        end
        if candidates[2] < max_idx
            if y > y_prev # 증가하는 중
                idx = candidates[2]
            else
                idx = candidates[1]
            end
        else # candidates[2] >= max_idx
            if y > y_prev # 증가하는 중
                idx = candidates[1]
            else
                idx = candidates[2]
            end
        end
    end

    # 4) 연속 t 구간에서 이분법 수행 (t는 실수로 계산)
    a = float(idx)
    b = float(idx + 1)

    g(t) = target_function(t) - y

    g_a = g(a)
    g_b = g(b)
    if g_a * g_b > 0
        error("선택한 구간에 근이 없습니다. (부호 동일)")
    end

    iter = 0
    while (b - a) > tol && iter < max_iter
        m = (a + b) / 2
        g_m = g(m)

        if abs(g_m - 0) < tol
            println("iter : $(iter)")
            return m
        end

        if g_a * g_m < 0
            b = m
            g_b = g_m
        else
            a = m
            g_a = g_m
        end

        iter += 1
    end
    println("iter : $(iter)")
    return (a + b) / 2
end

target_function(t) = calc_theta_sim(
    predicted_r_values[1], 
    predicted_r_values[2], 
    predicted_r_values[3], 
    predicted_r_values[4], 
    deg2rad(13), 
    deg2rad(248) + (t - 1) * 2 * pi / 100
)

# what is next step?
for y in rad2deg(findmin(data)[1]) : rad2deg(findmax(data)[1])
    y_prev = y+0.5 # degree # 감소 중

    t_est = bisection_method(target_function, data, deg2rad(y), deg2rad(y_prev))
    println("t ≈ ", t_est)
    println("f(t) ≈ ", target_function(t_est), " rad  (", rad2deg(target_function(t_est)), " deg)")
    println("---------------")
end



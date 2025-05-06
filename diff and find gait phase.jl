using DelimitedFiles, Plots, Random

optimal_theta5_output_file_path = "C:\\Users\\Jehyeon\\OneDrive - GIST\\바탕화면\\GIST\\4-bar linkage\\julia_code\\data\\optimal_theta5_output.txt" # 최적화 전. matlab에서

output = readdlm(optimal_theta5_output_file_path)
output = vec(output)
output = rad2deg.(output)

gait_cycle = 0:1:100

function calc_gradient(x::AbstractVector, y::AbstractVector)
    n = length(y)
    grad = zeros(n)
    # 첫 번째 포인트: 순방향 차분법
    grad[1] = (y[2] - y[1]) / (x[2] - x[1])
    # 중앙 차분법: 2번째부터 n-1번째 포인트까지
    for i in 2:n-1
        grad[i] = (y[i+1] - y[i-1]) / (x[i+1] - x[i-1])
    end
    # 마지막 포인트: 후방 차분법
    grad[n] = (y[n] - y[n-1]) / (x[n] - x[n-1])
    return grad
end

# 기울기 계산
grad_output = calc_gradient(collect(gait_cycle), output)

plot(grad_output, label="grad_output", color="red", linewidth=2, xlabel="% gait cycle", ylabel="hip angular velocity(degree/index)")


# function find_gait_phase(target_angle::Float64, gait_cycle::AbstractVector, output::AbstractVector, grad_output::AbstractVector; tol=1e-3, prev_angle::Union{Float64,Nothing}=nothing)
#     # target_angle와 차이가 tol 이하인 인덱스를 후보로 추출
#     candidate_indices = findall(i -> abs(output[i] - target_angle) < tol, 1:length(output))
    
#     if isempty(candidate_indices)
#         println("Target angle not found within tolerance.")
#         return nothing
#     elseif length(candidate_indices) == 1
#         # 후보가 한 개이면 해당 gait phase 반환
#         return gait_cycle[candidate_indices[1]]
#     else
#         # 후보가 여러 개이면 이전 hip angle을 활용하여 적합한 분기를 선택
#         if prev_angle !== nothing
#             # 상승 중이면 target_angle > prev_angle → 기울기가 양이어야 함
#             # 하강 중이면 target_angle < prev_angle → 기울기가 음이어야 함
#             expected_sign = (target_angle - prev_angle) > 0 ? 1 : -1
#             valid_candidates = [i for i in candidate_indices if sign(grad_output[i]) == expected_sign]
#             if !isempty(valid_candidates)
#                 # 여러 후보 중, target과 가장 가까운 값을 선택
#                 best = valid_candidates[argmin(abs.(output[valid_candidates] .- target_angle))]
#                 return gait_cycle[best]
#             else
#                 # 기대하는 기울기 부호 후보가 없으면, 전체 후보 중 최소 오차 선택
#                 best = candidate_indices[argmin(abs.(output[candidate_indices] .- target_angle))]
#                 return gait_cycle[best]
#             end
#         else
#             # 이전 hip angle 정보가 없다면, 오차가 가장 작은 후보 선택
#             best = candidate_indices[argmin(abs.(output[candidate_indices] .- target_angle))]
#             return gait_cycle[best]
#         end
#     end
# end


# #degree
# target_angle = 10.0
# prev_angle = 28.0
# println("target angle : $(target_angle) degree")

# index = find_gait_phase(target_angle, gait_cycle, output, grad_output; tol=1, prev_angle=prev_angle)
# println("찾은 gait phase index: ", index)


# Random.seed!(123)  # 시드를 고정하면 매 실행마다 같은 난수 시퀀스를 생성합니다.
# # --------------------------
# # find_gait_phase 함수
# function find_gait_phase(target_angle::Float64, gait_cycle::AbstractVector, 
#                          output::AbstractVector, grad_output::AbstractVector;
#                          tol=1e-3, prev_angle::Union{Float64,Nothing}=nothing)
#     candidate_indices = findall(i -> abs(output[i] - target_angle) < tol, 1:length(output))
    
#     if isempty(candidate_indices)
#         return nothing
#     elseif length(candidate_indices) == 1
#         return gait_cycle[candidate_indices[1]]
#     else
#         if prev_angle !== nothing
#             expected_sign = (target_angle - prev_angle) > 0 ? 1 : -1
#             valid_candidates = [i for i in candidate_indices if sign(grad_output[i]) == expected_sign]
#             if !isempty(valid_candidates)
#                 best = valid_candidates[argmin(abs.(output[valid_candidates] .- target_angle))]
#                 return gait_cycle[best]
#             else
#                 best = candidate_indices[argmin(abs.(output[candidate_indices] .- target_angle))]
#                 return gait_cycle[best]
#             end
#         else
#             best = candidate_indices[argmin(abs.(output[candidate_indices] .- target_angle))]
#             return gait_cycle[best]
#         end
#     end
# end

# # --------------------------
# # 랜덤 target angle 50개 생성 (output의 min~max 사이, 단위: degree)
# min_val = minimum(output)
# max_val = maximum(output)
# num_targets = 1000
# target_angles = [min_val + (max_val - min_val) * rand() for _ in 1:num_targets]

# # 각 target에 대해, prev_angle을 50% 확률로 제공
# # 제공되는 경우, target의 ±10% 범위 내에서 랜덤하게 설정
# prev_angles = [ rand() < 0.5 ? nothing : ta * (1 + (rand()*0.2 - 0.1)) for ta in target_angles ]

# # 허용 오차 tol (degree)
# tol = 0.01

# # --------------------------
# # 각 target에 대해 find_gait_phase 호출 후 성공/실패 판별 및 결과 저장
# # 주의: gait_cycle은 0부터 시작하므로, 결과로 나온 gait phase 값은 실제 % 값입니다.
# success_targets = Float64[]
# success_phases = Float64[]
# failure_targets = Float64[]
# failure_phases = Float64[]

# num_success = 0
# num_failure = 0

# for (ta, pa) in zip(target_angles, prev_angles)
#     phase = find_gait_phase(ta, gait_cycle, output, grad_output; tol=tol, prev_angle=pa)
#     if phase === nothing
#         num_failure += 1
#         push!(failure_targets, ta)
#         push!(failure_phases, NaN)
#         println("Target angle $(round(ta,digits=2))° : No gait phase found.")
#     else
#         # gait_cycle의 값는 0~100, 이 값을 그대로 사용 (예, 37이면 37% gait cycle)
#         idx = phase + 1  # Julia 배열은 1-indexed (gait_cycle[1]==0)
#         found_val = output[idx]
#         err = abs(found_val - ta)
#         if err <= tol
#             num_success += 1
#             push!(success_targets, ta)
#             push!(success_phases, phase)
#         else
#             num_failure += 1
#             push!(failure_targets, ta)
#             push!(failure_phases, phase)
#         end
#         println("Target angle: $(round(ta,digits=2))° | Prev_angle: $(pa === nothing ? "none" : string(round(pa,digits=2)) ) | Found phase: $(phase)% | Output: $(round(found_val,digits=2))° | error: $(round(err,digits=2))°")
#     end
# end

# println("\n총 타겟 개수: $(num_targets)")
# println("성공: $(num_success), 실패: $(num_failure)")

# # NaN이 아닌 실패 데이터만 필터링합니다.
# filtered_failure_phases = [p for p in failure_phases if !isnan(p)]
# filtered_failure_targets = [failure_targets[i] for i in 1:length(failure_phases) if !isnan(failure_phases[i])]

# # 기존 플롯(가로: gait phase, 세로: target hip angle)
# plt = scatter(success_phases, success_targets, 
#         marker=:circle, markersize=2, markercolor=:green,
#         label="success (error ≤ $(tol)°)", legend=:bottomright,
#         grid=true, xlabel="% Gait Phase", ylabel="Target Hip Angle (°)",
#         title="Gait Phase vs. Target Hip Angle")

# # 실패 데이터가 존재하면 빨간 다이아몬드로 추가
# if !isempty(filtered_failure_phases)
#     scatter!(plt, filtered_failure_phases, filtered_failure_targets, 
#             marker=:diamond, markersize=10, markercolor=:red,
#             label="fail")
# end

# display(plt)

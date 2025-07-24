using ForwardDiff, Plots, Statistics, DelimitedFiles
plotlyjs()


# 분리된 gait cycles 불러오기
# 1) 파일에서 모든 줄을 문자열로 읽어 들인다

subject_number = 1
gait_cycles_length = 128
lines = readlines("result/subject$(subject_number)/gait_cycles_$(gait_cycles_length)ea.txt")

# 2) 각 줄을 빈칸(또는 탭)으로 split 하고, Float64 로 parse
gait_cycles = [ parse.(Float64, split(line)) for line in lines ]

println("Loaded $(length(gait_cycles)) gait cycles.")
# 확인: 각 cycle 의 길이가 다양함을 볼 수 있음
# for (i,c) in enumerate(gait_cycles)
#     @info "cycle $i length = $(length(c)) samples"
# end


# 최적화된 4-바 링크 파라미터들 불러오기
# 1) 파일에서 읽어와서 1×5 Matrix 로 들어오므로 벡터로 변환
# optimal_r = vec(readdlm("result/optimal_r_values.txt", ','))

# # 2) 변수에 재할당
# r1, r2, r5, r6, θ1 = optimal_r

# println("Loaded: ", "r1=", r1, ", r2=", r2, ", r5=", r5, 
#         ", r6=", r6, ", th1(rad)=", θ1)
    
# 로봇 링크 파라미터
# 이미 최적화 한 값을 사용해야 하나???
# 이미 최적화 된 값을 사용해야, 인간과 유사하게 걷는 4-바 링크의 수식을 이용해서
# thigh angle을 통해 gait phase를 구했다고 말 할 수 있지 않나?



# θ2 → θ5 함수 (라디안 단위)
# 각도 : (+x축) 으로부터 CCW방향으로 증가
function f(θ2)
    θs = atan(r1*sin(θ1)+r2*sin(θ2),
               r1*cos(θ1)+r2*cos(θ2))
    rs = sqrt(r1^2 + r2^2 + 2r1*r2*cos(θ1-θ2))
    arg = clamp((-r6^2 + r5^2 + rs^2) / (2r5*rs), -1.0, 1.0)
    θ5 = θs - acos(arg) + π
    return θ5
end

# 완전 수학적으로 계산
# aaa와 bbb가 계속 무한으로 출력되는 문제 상황
function invert_newton_raphson(θ5_target_rad, θ2₀_rad; tol=1e-8, max_iter=1000)
    θ2 = θ2₀_rad # 초기값 설정
    for k in 1:max_iter 
        Fv = f(θ2) - θ5_target_rad 
        Jv = ForwardDiff.derivative(f, θ2)

        # --- 도함수가 0에 가까운 경우 처리 ---
        # 도함수가 0이면 뉴턴 스텝(Fv / Jv) 계산이 불가능합니다.
        # 이런 경우 일반적으로 수렴 실패로 처리하고 중단하는 것이 논리적입니다.
        # 부동 소수점 계산에서는 0과 정확히 일치하는지 대신 아주 작은 값보다 작은지로 판단합니다.
        if abs(Jv) < 1e-10 # 예: 1e-10보다 작으면 0에 가깝다고 판단
            error("Derivative is too close to zero ($(Jv)) at iteration $(k), θ2 = $(rad2deg(θ2))°. Cannot proceed with Newton-Raphson.")
        end
        println("aaa")
        Δ  = Fv / Jv
        θ2 -= Δ
        println("bbb")
        if abs(Δ) < tol
            return θ2
        end
    end
    error("Newton-Raphson did not converge within $(max_iter) iterations. Final step size: $(abs(Δ))")
end

function calc_forward_data() # θ2 로 θ5를 계산하는 순방향 함수
    θ2_deg = collect(1:360) # array 만들어주는 함수
    θ2_rad = deg2rad.(θ2_deg)
    θ5_rad = f.(θ2_rad) # 함수 f를 브로드캐스트 호출하여 θ5_rad 배열 생성
    θ5_deg = rad2deg.(θ5_rad)

    return θ2_deg, θ5_deg
end

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

function normalization()

    cycle = gait_cycles[1]
    println("size of cycle : $(size(cycle)[1])")

    phase = collect(LinRange(0, 100, length(cycle))) # 0~100을 cycle 길이만큼 균등 분할
    println(phase)
    
    p = plot(
      phase, cycle,
      xlabel = "Gait Cycle (%)",
      ylabel = "Thigh Angle (degree)",
      title  = "Normalized 1st Gait Cycle"
    )

    display(p)
    return p
end



# # 인자
# - `gait_cycles` : Vector{Vector{Float64}} (각 element 가 한 사이클)
# - `N_phase`     : Int  (재샘플링할 위상 분할 개수, 기본 101)

# # 반환
# - `phase`      :: Vector{Float64} — 0…100 까지 균일 분할된 위상
# - `mean_cycle` :: Vector{Float64} — 위상별 평균 thigh-angle
function mean_gait_cycle(gait_cycles::Vector{Vector{Float64}}; N_phase::Int64=101)
    # 1) 공통 위상 축 (0…100% N_phase 점)
    phase = collect(LinRange(0, 100, N_phase))

    # 2) 개별 사이클을 위상 축에 선형 보간 재샘플링
    function resample_cycle(cycle)
        orig_ph = collect(LinRange(0, 100, length(cycle)))
        out = similar(phase)
        for (i, p) in enumerate(phase)
            j = searchsortedfirst(orig_ph, p)
            if j ≤ 1
                out[i] = cycle[1]
            elseif j > length(cycle)
                out[i] = cycle[end]
            else
                p1, p2 = orig_ph[j-1], orig_ph[j]
                y1, y2 = cycle[j-1], cycle[j]
                out[i] = y1 + (p - p1)/(p2 - p1)*(y2 - y1)
            end
        end
        return out
    end

    # 3) 모든 사이클 재샘플링 → 행렬로
    R = hcat(resample_cycle.(gait_cycles)...)

    # 4) 위상별 평균 (각 행의 평균)
    mean_cycle = vec(mean(R; dims=2))
    writedlm("result/subject$(subject_number)/mean_cycle.txt", mean_cycle, ',')
    println("Saved mean_cycle.txt")

    # 5) 결과 플롯
    p = plot(phase, mean_cycle,
         xlabel="Gait Cycle (%)",
         ylabel="Average Thigh Angle (degree)",
         title ="Mean Gait Cycle over $(length(gait_cycles)) trials",
         linewidth=2)
    display(p)
    return phase, mean_cycle
end




# ─────────────────────────────────────────────────────────────────────────
# 특정 thigh angle 에 대응하는 gait phase (%) 반환 함수
#
# cycle     :: Vector{Float64}   — 한 사이클(θ_filt[peak[i]:peak[i+1]])  
# target    :: Float64           — 찾고 싶은 thigh angle (degree)
# return    :: Float64           — 0…100 범위의 gait phase


# 지금 코드는 특정 target에 대해 두 개의 정답이 있는 상황에서 이를 구별하지 못함.
#─────────────────────────────────────────────────────────────────────────
function phase_for_angle(cycle::Vector{Float64}, target::Float64)
    # 1) phase 벡터 (0…100, cycle 길이에 맞춰)
    phases = collect(LinRange(0, 100, length(cycle)))

    # 2) 선형 보간으로 target 이 끼어드는 구간 찾기
    for i in 1:length(cycle)-1
        y1, y2 = cycle[i], cycle[i+1]
        println("y1 = $(y1), y2 = $(y2), i=$(i)")
        if (y1 - target)*(y2 - target) ≤ 0    # cross or touch
            t1, t2 = phases[i], phases[i+1]
            t = t1 + (target - y1)*(t2 - t1)/(y2 - y1) # linear interpolation
            println("t1 = $(t1), t2 = $(t2), t = $(t)")
            return t
        end
    end

    # 범위 밖이면 에러 또는 clipping
    error("phase_for_angle: target = $target outside [$(minimum(cycle)), $(maximum(cycle))]")
end
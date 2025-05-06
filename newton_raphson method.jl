using ForwardDiff, Plots, Statistics
plotlyjs()

# 로봇 링크 파라미터
# 이미 최적화 한 값을 사용해야 하나???
const r1, r2, r5, r6 = 14, 3, 12, 21
const θ1 = deg2rad(13)

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

function invert_newton_raphson(θ5_target_rad, θ2₀_rad; tol=1e-8, max_iter=1e+3)
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

        Δ  = Fv / Jv
        θ2 -= Δ
        
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

function calc_inverse_data() # 특정 θ5로 θ2를 역으로 계산하는 역방향 함수

    θ2_deg, θ5_deg = calc_forward_data()
    θ5_deg_grad = calc_gradient(θ2_deg, θ5_deg)
    slope_sign = sign.(θ5_deg_grad)
    N = length(θ2_deg)

    # plot(θ2_deg, θ5_deg)
    # plot(θ2_deg, θ5_deg_grad, color="blue")

    # 하나의 theta5에 대해 두 개 이상의 theta2가 대응되는 경우가 있다.
    # 이는 invert_newton_raphson 함수의 입력 파라미터인, th2 init 값을 어떻게 설정하느냐에 달려있다.
    # 원래 처음 생각은 글로벌 미니멈을 기준으로 좌, 우로 그냥 90도, 270를 초기값으로 주려고 생각했다.

    # 하지만, 이 방법이 아니라, 기울기를 고려하는 방법은 어떨까?
    # 일단 애초에 calc_forward_data 의 출력인 th5_deg의 기울기를 구해 놓는 것이다. 
    # 그리고 원함수와 도함수의 그래프를 한번에 띄워서 알아보기 쉽게 하고.
    # 그 후, 이제 특정한 theta5 값을 입력한다면?   
    # 그 theta5 값에 해당하는, 미리 구해논 기울기를 읽어 오고(실시간으로 읽어도 될까)
    # 그 그 기울기가 음수, 양수임에 따라서 만약 음수면 th2 초기값을 90도로, 양수면 270도로?
    
    # 근데 위 두 방법은 같은 거 아닌가?
    # 첫 번째 방법은 th5의 맥스값을 구하려면 모든 값을 비교해야 해서, 두 번째 방법이 계산상 간단한가?

    θ2_est_deg = zeros(Float64, N)
    first_init_deg   = slope_sign[1] > 0 ?  90.0 : 270.0
    θ2_est_deg[1]    = rad2deg(     
                          invert_newton_raphson(
                            deg2rad(θ5_deg[1]),
                            deg2rad(first_init_deg)
                          )
                       )
    for i in 2:N
        θ5_tar_rad    = deg2rad(θ5_deg[i])
        θ2_guess_rad  = deg2rad(θ2_est_deg[i-1])      # ← 이전 스텝 추정값
        θ2_est_deg[i] = rad2deg( invert_newton_raphson(θ5_tar_rad, θ2_guess_rad) )
    end
    err_pct = abs.((θ2_est_deg .- θ2_deg) ./ θ2_deg) * 100
    @show mean(err_pct)
    p = plot(θ2_deg, θ5_deg;       label="forward θ5", xlabel="θ2 (deg)")
    plot!(θ2_deg, θ2_est_deg;  label="inverse θ2_est", lw=2, ls=:dash)
    display(p)
    return θ2_est_deg
end
# using CSV, DataFrames, Plots, Statistics

# FILENAME = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\gaitdata\\imu_data\\1\\imu_trial_1.csv"

# df = CSV.read(FILENAME, DataFrame)   # df[!, "accX1"]처럼 이름으로 접근 가능
# n_samples = nrow(df)
# # data 개수 : 3416 개

# # 가속도·자이로 열 추출
# a_x  = df[!, "accX1"]    # m/s^2
# a_y  = df[!, "accY1"]    # 필요하면 a_z도
# ω_z  = df[!, "gyroZ1"]   # gyroscope Z축 각속도 deg/s
# angle1 = df[!, "angle1"]

# t     = df[!, "time"]
# DT    = mean(diff(t))          # 평균 샘플 간격
# println("샘플 주기 ≈ $(round(DT*1000, digits=2)) ms")

# ALPHA = 0.96        # 0.98~0.995 정도에서 튜닝
# θ      = zeros(Float64, n_samples)
# θ[1]   = atan(a_x[1], a_y[1])   # rad

# for k in 2:n_samples
#     θ_acc  = atan(a_x[k], a_y[k])
#     θ_gyro = θ[k-1] + ω_z[k]*DT        # ω_z가 rad/s라면 그대로, deg/s면 *π/180
#     θ[k]   = ALPHA*θ_gyro + (1-ALPHA)*θ_acc
# end

# θ_deg = θ .* (180/π)

# plot(a_x) 랑 plot(a_y) 했을때 노이즈가 심하다
# theta degree 출력하면 값이 무슨 1000 이상 까지도 나온다?
# 단위 문제, atan 순서 문제, DT는 상수로 바꿔라?, ALPHA 문제?

# 이거 리얼타임아니냐


# julia> df[1, "accX1"]
# 9.82197

# julia> df[1, "accY1"]
# 0.874182

# julia> df[1, "gyroZ1"]
# -2.16463

using CSV, DataFrames, Plots, Statistics
plotlyjs()
# ── 0. 데이터 로드 ───────────────────────────────────────────
FILENAME = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\gaitdata\\imu_data\\1\\imu_trial_1.csv"

df = CSV.read(FILENAME, DataFrame)   # df[!, "accX1"]처럼 이름으로 접근 가능
t        = df.time                              # time [s]
n        = length(t)
DT       = mean(diff(t))                        # 평균 샘플 주기
ALPHA    = 0.96                                 # complementary filter 계수
DEG2RAD  = π/180
RAD2DEG  = 180/π

# ★ LPF 설정
fc   = 3.0                     # Hz (2~4 사이에서 튜닝)
τ    = 1/(2π*fc)               # time constant (≈ 0.08 s at 2 Hz)

# # 그래프 초기화 (GR 백엔드가 가장 가볍습니다)
# plt = plot(xlabel = "time [s]",
#            ylabel = "thigh angle [deg]",
#            ylim   = (-60, 120),
#            legend = false)
# display(plt)

angle_deg = Float64[]    # 누적 각도 저장용
time_buf  = Float64[]

# 버퍼 비우기
empty!(time_buf)          # = deleteat!(time_buf, :)
empty!(angle_deg)

θ_prev = atan(df[1, :accX1], df[1, :accY1])   # 초기 각도(rad)
θ_lpf  = θ_prev                                # ★ LPF 내부 상태(rad)
θ_prev_lpf = θ_lpf

push!(angle_deg, θ_prev * RAD2DEG)
push!(time_buf, df[1, :time])

for k in 2:nrow(df)

    # 1) 한 줄씩 읽기
    ax = df[k, :accX1]
    ay = df[k, :accY1]
    wz = df[k, :gyroZ1] * DEG2RAD       # deg/s → rad/s  (이미 rad/s면 변환 생략)

    # ★ per-sample dt 사용 (LPF/보정 모두에 권장)
    dt = df[k, :time] - df[k-1, :time]


    # 2) 보정 각 계산
    θ_acc  = atan(ax, ay)       # 자세축 확인 후 필요하면 부호 바꾸세요
    θ_gyro = θ_prev + wz*dt
    θ      = ALPHA*θ_gyro + (1-ALPHA)*θ_acc
    
    # ★ 3) 1차 LPF (저지연)
    α_lpf  = dt/(τ + dt)
    θ_lpf  = θ_lpf + α_lpf*(θ - θ_lpf)    # ← 부드러운 각도(rad)

    # ★ 다음 스텝을 위해 이전 상태 업데이트는 '부드럽힌 값'을 사용
    # θ_prev = θ_lpf

    # 3) 버퍼에 쌓기
    push!(time_buf, df[k, :time])
    # push!(angle_deg, θ * RAD2DEG)
    push!(angle_deg, θ_lpf * RAD2DEG)

    # println(angle_deg) # 프린트를 해서 오래걸린 거구나
    # # plot!(plt, time_buf, angle_deg, clear = true)   # 덮어쓰기 방식이 가장 쉽습니다
    # if k % 20 == 0          # 10샘플마다만 재도장
    #     plot!(plt, time_buf, angle_deg, clear=true) # 그래프 업데이트
    #     # Plots.flush()
    #     # gui()
    #     yield()
    # end

    # sleep(DT)   # 실제 센서와 똑같은 속도로 “재생”

    θ_prev = θ
    θ_prev_lpf = θ_lpf

end

# 그래프 초기화 (GR 백엔드가 가장 가볍습니다)
plt = plot(angle_deg,
           title = "real-time thigh angle from imu",
           xlabel = "index",
           ylabel = "thigh angle [deg]",
           ylim   = (-60, 120),
           legend = true)
hline!([0.476826344703436 * RAD2DEG, 0.476826344703436 * RAD2DEG], label = "max value of data")
hline!([-0.49388091910725773 * RAD2DEG, -0.49388091910725773 * RAD2DEG], label = "min value of data")
display(plt)


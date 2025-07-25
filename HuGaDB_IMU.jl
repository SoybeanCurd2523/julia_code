# HuGaDB_IMU.jl

using DelimitedFiles, Statistics, Plots
plotlyjs()

# ---- 상수/공통 파라미터(원하면 여기만 수정) -----------------
FS         = 56.3500           # Hz
DT         = 1 / FS
ACC_SENS   = 32768.0
ACC_RANGE  = 2.0 * 9.80665
GYRO_SENS  = 32768.0
GYRO_RANGE = 2000.0

CUTOFF = 3
ALPHA = 0.96
MIN_CYCLE_TIME = 0.6

subject_number = 11
println("subject_number : ", subject_number)

function main()

    # 1. 설정 값 ----------------------------------------------------

    if subject_number < 10
        FILENAME = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\HuGaDB\\Data\\HuGaDB_v1_walking_0$(subject_number)_00.txt"
    else
        FILENAME = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\HuGaDB\\Data\\HuGaDB_v1_walking_$(subject_number)_00.txt"
    end

    #=
    논문 : Data were collected with accelerometer’s 
    range equal to ±2 g with sensitivity 16.384 LSB/g and gyroscope’s range equal to 
    ±2000◦ /s with sensitivity 16.4 LSB /◦ /s. All sensors are powered from a battery, 
    that helps to minimize electrical grid noise.
    Accelerometer and gyroscope signals were stored in int16 format. EMG sig- 
    nals are stored in uint8. Therefore, accelerometer data can be converted to m/s2 
    by dividing raw data 32768 and multiplying it by 2g. Raw gyroscope data can 
    be converted to ◦ /s by multiplying it by 2000/32768.

    raw data는 int16이므로, -2^15(-32768), 2^15-1(32767) 사이의 값을 가진다.
    가속도계의 측정범위는 +- 2g이고, 자이로스코프 측정 범위는 +-2000◦ /s 이다.
    가속도 [m/s^2] = raw_acc * (2.0 * 9.81) / 32768.0
    각속도 [°/s] = raw_gyro * 2000.0 / 32768.0
    =#


    # 2. 헤더/데이터 로드 함수 ----------------------------------------
    function load_hugadb(filename)
        header = String[]
        open(filename, "r") do io
            for line in eachline(io)
                if !startswith(line, "#") # 로 시작하는 헤더 줄을 건너뛰고, 첫 non-# 줄을 칼럼명(header)으로
                    header = split(chomp(line), '\t')  # 그 뒤의 숫자 데이터를 전부 읽어 반환
                    break
                end
            end
            return header, readdlm(io, '\t')
        end
    end

    # 3. 안전한 칼럼 인덱스 함수 -------------------------------------
    function colidx(colname, header) # header 배열에서 원하는 열(colname)이 몇 번째인지 찾아줌
        idx = findfirst(isequal(colname), header)
        idx === nothing && error("Column '$colname' not found in header")
        return idx
    end

    # 4. 헤더·데이터 읽기 --------------------------------------------
    header, raw = load_hugadb(FILENAME)
    n_samples   = size(raw, 1) #5989

    # 5. 칼럼 인덱스 찾기 ------------------------------------------
    ix_acc_x  = colidx("acc_lt_x", header)
    iz_acc_z  = colidx("acc_lt_z", header)
    iy_acc_y  = colidx("acc_lt_y", header)

    ix_gyro_x = colidx("gyro_lt_x", header)
    iz_gyro_z = colidx("gyro_lt_z", header)
    iy_gyro_y = colidx("gyro_lt_y", header)


    # 6. 물리량 보정 및 complementary filter ------------------------

    # 가속도(raw) → m/s^2 단위, 자이로(raw) → rad/s 단위로 변환
    a_x = raw[:, ix_acc_x]  .* ACC_RANGE ./ ACC_SENS
    a_z = raw[:, iz_acc_z]  .* ACC_RANGE ./ ACC_SENS
    a_y = raw[:, iy_acc_y]  .* ACC_RANGE ./ ACC_SENS

    ω_x = raw[:, ix_gyro_x] .* GYRO_RANGE ./ GYRO_SENS .* (π/180)
    ω_z = raw[:, iz_gyro_z] .* GYRO_RANGE ./ GYRO_SENS .* (π/180)
    ω_y = raw[:, iy_gyro_y] .* GYRO_RANGE ./ GYRO_SENS .* (π/180)

    θ = zeros(Float64, n_samples)
    θ[1] = atan(a_z[1], a_x[1])

    for k in 2:n_samples
        θ_acc  = atan(a_z[k], a_x[k])
        θ_gyro = θ[k-1] + ω_y[k]*DT
        
        # “자이로에는 고역통과 필터(HPF), 가속도에는 저역통과 필터(LPF)를 각각 씌운 뒤 합친다”를 
        # 암묵적(implicit) 으로 구현한 것
        θ[k]   = ALPHA*θ_gyro + (1-ALPHA)*θ_acc
    end

    # 7. 라디안→도 & 저장 -------------------------------------------
    θ_deg = θ .* (180/π)
    # writedlm("thigh_angle_deg.txt", θ_deg)
    # println("완료: thigh_angle_deg.txt 에 $n_samples 개 각도값 저장됨.")

    # 8. 부호 반전 (flexion을 양수로 보기 위해)  <<< 잘못된 듯.
    # θ_calib = -θ_deg
    θ_calib = θ_deg
    # writedlm("thigh_angle_calib.txt", θ_calib)

    N = n_samples                         # 5989개 데이터 중 첫 몇 개만 그래프에 표시할건지
    t = (0:N-1) .* DT                # 0초부터 (N−1)*DT 초까지

    plot(
    t, θ_calib[1:N],
    label  = "α = $(ALPHA)",
    xlabel = "Time (s)",
    ylabel = "Thigh Angle (degree)",
    legend = (0.09, 0.23),
    title  = "Thigh Angle (complementary filtered), First $(N) samples",
    xlabelfontsize  = 14,   # x축 레이블 크기
    ylabelfontsize  = 14,   # y축 레이블 크기
    titlefontsize   = 16,   # 제목 글꼴 크기
    legendfontsize  = 12    # 범례 글꼴 크기
    )

    # --- 9. 이동평균으로 노이즈 억제 (cut-off =  Hz) ---------------------

    window = Int(round(FS / cutoff))     # 윈도우 길이 계산
    θ_filt = similar(θ_calib)
    running = 0.0
    for k in 1:n_samples
        running += θ_calib[k]
        if k > window
            running -= θ_calib[k-window]
            θ_filt[k] = running / window
        else
            θ_filt[k] = running / k
        end
    end

    # 10. 결과를 파일에 저장
    θ_filt_outfile = "result/subject$(subject_number)/thigh_angle_filt_$(cutoff)Hz.txt"
    writedlm(θ_filt_outfile, θ_filt)

    println("완료: 컷오프 빈도 $(cutoff) Hz 로 필터링한 각도값을 $θ_filt_outfile 에 저장했습니다.")
    plot(
    t, θ_filt[1:N],
    label  = "cutoff = $(cutoff) Hz",
    xlabel = "Time (s)",
    ylabel = "Thigh Angle (degree)",
    legend = (0.09, 0.20),
    title  = "Thigh Angle (Moving‐Average filtered), First $(N) samples",
    xlabelfontsize  = 14,   # x축 레이블 크기
    ylabelfontsize  = 14,   # y축 레이블 크기
    titlefontsize   = 16,   # 제목 글꼴 크기
    legendfontsize  = 10    # 범례 글꼴 크기
    )


    # --- 11. 각 gait cycle 분리--------------------------

    # --------------------------------------------------
    # 1) 파라미터 설정

    min_dist = round(Int, min_cycle_time * FS) # → min_dist = 0.6초 × 56.35Hz ≈ 34 샘플

    # 2) 전 구간에 대해 로컬 피크 후보 검출
    candidates = [ i for i in 2:n_samples-1 if
                    θ_filt[i] > θ_filt[i-1] &&
                    θ_filt[i] > θ_filt[i+1]
                ]

    # --------------------------------------------------
    # 3) 최소 간격(min_dist) 필터 적용해서 진짜 피크만 selection
    peaks = Int[] # 각 사이클의 시작 인덱스
    for idx in candidates
        if isempty(peaks) || idx - peaks[end] >= min_dist
            push!(peaks, idx)
        end
    end

    # println("총 $(length(peaks))개의 gait-cycle 대표 피크 인덱스 발생")

    # --------------------------------------------------
    # 4) 피크 사이사이를 하나의 cycle 로 분리
    # 예)
    # 1번 사이클 : 81.68 ~ 61.50
    # 2번 사이클 : 61.50 ~ 49.21
    # 3번 사이클 : 49.21 ~ 49.99
    # 경계값을 공유
    gait_cycles = [
        θ_filt[ peaks[i] : peaks[i+1] ]
        for i in 1:length(peaks)-1
    ]

    # println("총 $(length(gait_cycles))개의 gait cycle 분리 완료")


    gait_cycles_outfile = "result/subject$(subject_number)/gait_cycles_$(length(gait_cycles))ea.txt" # Vector{Vector{Float64}}
    writedlm(gait_cycles_outfile, gait_cycles)

    println("완료: 총 $(length(gait_cycles))개의 분리된 gait cycle들을 $gait_cycles_outfile 에 저장했습니다.")
    # --------------------------------------------------
    # 5) (선택) 각 cycle 별로 다른 색으로 전체 신호 위에 표시

    p = plot(
        # t, θ_filt[1:N], # 가로축을 시간으로
        θ_filt[1:N],      # 가로축을 인덱스로
        color = :black, alpha = 0.3, legend = false,
        # xlabel = "Tims (s)", 
        xlabel = "index",
        ylabel = "Thigh Angle (degree)",
        title  = "분리된 Gait Cycles"
    )
    for (i, c) in enumerate(gait_cycles)
        # 피크 i 에서 다음 피크로 가는 구간의 인덱스
        idx0, idx1 = peaks[i], peaks[i+1]
        plot!(
            p,
            # t[idx0:idx1], θ_filt[idx0:idx1],
            idx0:idx1, θ_filt[idx0:idx1],
            label = "cycle $i",
            linewidth = 2,
        )
    end

    display(p)

    # newton_raphson method.jl 파일에서 복사해옴

    # gait_cycles_length = length(gait_cycles)
    # lines = readlines("result/subject$(subject_number)/gait_cycles_$(gait_cycles_length)ea.txt")

    # # 2) 각 줄을 빈칸(또는 탭)으로 split 하고, Float64 로 parse
    # gait_cycles = [ parse.(Float64, split(line)) for line in lines ]


    # # # 인자
    # # - `gait_cycles` : Vector{Vector{Float64}} (각 element 가 한 사이클)
    # # - `N_phase`     : Int  (재샘플링할 위상 분할 개수, 기본 101)

    # # # 반환
    # # - `phase`      :: Vector{Float64} — 0…100 까지 균일 분할된 위상
    # # - `mean_cycle` :: Vector{Float64} — 위상별 평균 thigh-angle
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
        return
    end

    mean_gait_cycle(gait_cycles)
end

main()


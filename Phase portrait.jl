using DelimitedFiles, Plots
plotlyjs()

theta5_output_file_path = "C:\\Users\\Jehyeon\\Dropbox\\바탕 화면\\GIST\\4-bar linkage\\julia_code\\result\\subject11\\mean_cycle.txt"

output = readdlm(theta5_output_file_path)
output = vec(output)

# 각속도 계산 (단위 시간 간격 기준)
theta5_dot = diff(output)                # θ₅[i+1] - θ₅[i], 길이: 100 각도 변화량
theta5_mid = output[1:end-1]             # 길이: 100. output을 길이를 100으로 맞춘거

# Phase Portrait 그리기
plot(theta5_mid, theta5_dot,
     xlabel="θ₅ (degree)",
     ylabel="δθ₅ (degree)",
     title="Phase Portrait of θ₅",
     legend=false,
     linewidth=2,
     aspect_ratio=:equal)

scatter!([theta5_mid[1]], [theta5_dot[1]], color="red", marker=:circle, label="start")

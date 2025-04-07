using LsqFit, Plots, DelimitedFiles, ForwardDiff

# 데이터 불러오기
file_path = "C:\\Users\\Jehyeon\\Desktop\\GIST\\4-bar linkage\\julia_code\\data\\optimal_theta5_output.txt"
output = vec(readdlm(file_path))  # 101×1 → 1차원 벡터로

# 시간 벡터 (1부터 101까지)
t = collect(1:101)

# 모델 정의: θ(t) = A * sin(ω * t + φ) + B
sine_model(t, p) = @. p[1] * sin(p[2] * t + p[3]) + p[4]

# 초기 추정값 [A, ω, φ, B]
p0 = [1.0, 0.1, 0.0, mean(output)]

# 피팅 수행
fit = curve_fit(sine_model, t, output, p0)

# 결과 계산
fitted_output = sine_model(t, fit.param)

# 시각화
plot(t, output, label="original θ₅", lw=2)
plot!(t, fitted_output, label="sine fit", lw=2, ls=:dash, color=:red)


# # Phase portrait용 값 계산
# θ = fitted_output[1:end-1]              # θ₅(t)
# dθ = diff(fitted_output)                # Δθ₅(t) ≈ θ₅(t+1) - θ₅(t)

# # Phase Portrait 시각화
# plot(θ, dθ,
#      xlabel="θ₅ (rad)",
#      ylabel="Δθ₅ (rad)",
#      title="Phase Portrait of Fitted θ₅",
#      legend=false,
#      linewidth=2,
#      aspect_ratio=:equal)  # x, y축 스케일 동일하게

# # 시작점 표시
# scatter!([θ[1]], [dθ[1]], color="red", label="start", marker=:circle)


# 연속 함수로 모델 만들기
fitted_func(t) = fit.param[1] * sin(fit.param[2] * t + fit.param[3]) + fit.param[4]

t_vals = 1:0.1:101  # 더 촘촘하게 sampling
theta5_vals = [fitted_func(t) for t in t_vals]
theta5_dot_vals = [ForwardDiff.derivative(fitted_func, t) for t in t_vals]

# Plot phase portrait
plot(theta5_vals, theta5_dot_vals, xlabel="θ₅ (rad)", ylabel="dθ₅/dt (rad/frame)", title="Phase Portrait (ForwardDiff)", legend=false, linewidth=2)
scatter!([theta5_vals[1]], [theta5_dot_vals[1]], color="red", label="start")

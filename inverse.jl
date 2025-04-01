using Plots, Dierckx

r1 = 14.0
r2 = 3.0
r5 = 12.0
r6 = 21.0
th1 = deg2rad(13)

function cal_th5(r1, r2, r5, r6, th1, th2)
    rs = sqrt(r1^2 + r2^2 + 2*r1*r2*cos(th1 - th2))
    ths = atan(r1*sin(th1) + r2*sin(th2), r1*cos(th1) + r2*cos(th2))
    input_to_acos = (-r6^2 + r5^2 + rs^2) / (2 * r5 * rs)
    input_to_acos_clamped = clamp(input_to_acos, -1.0, 1.0)  # 안정성 위해
    
    th5 = π/2 - ( ths - acos(input_to_acos_clamped) + π )
    return th5
end

th2_init = deg2rad(248)
th2_values = [th2_init + (i-1)*2*pi/100 for i in 1:101]

# θ_5 = f(θ_2)
th5_values = [cal_th5(r1, r2, r5, r6, th1, th2) for th2 in th2_values]


idx1 = 1:47 # 감소
idx2 = 48:92 # 증가
idx3 = 93:101 # 감소

th5_values_1 = th5_values[idx1]
th2_values_1 = th2_values[idx1]

th5_values_2 = th5_values[idx2]
th2_values_2 = th2_values[idx2]

th5_values_3 = th5_values[idx3]
th2_values_3 = th2_values[idx3]



# θ_2 = f^-1 (θ_5)
th5_values_1_sorted = reverse(th5_values_1) # 보간 입력은 오름차순 이어야함
th2_values_1_sorted = reverse(th2_values_1)

itp_1 = Spline1D(th5_values_1_sorted, th2_values_1_sorted, k=3) # cubic spline


itp_2 = Spline1D(th5_values_2, th2_values_2, k=3) # cubic spline


th5_values_3_sorted = reverse(th5_values_3) # 보간 입력은 오름차순 이어야함
th2_values_3_sorted = reverse(th2_values_3)

itp_3 = Spline1D(th5_values_3_sorted, th2_values_3_sorted, k=3) # cubic spline


# th5_test = th5_values_1_sorted[25]
# th2_real = th2_values_1_sorted[25]
# th2_est = itp_1(th5_test)

# plot(th5_values_1_sorted, th2_values_1_sorted, label="θ₂ vs θ₅ (real)", lw=2, xlabel="θ₅_reverse (rad)", ylabel="θ₂_reverse (rad)", title="Inverse Function Approximation", legend=:bottomright)

# errors = abs.(th2_values_1_sorted .- itp_1.(th5_values_1_sorted))
# plot(th5_values_1_sorted, errors, label="|θ₂_real - θ₂_est|", xlabel="θ₅ (rad)", ylabel="Error (rad)", lw=2, title="Prediction Error")


# plot(th5_values_2, th2_values_2, label="θ₂ vs θ₅ (real)", lw=2, xlabel="θ₅ (rad)", ylabel="θ₂ (rad)", title="Inverse Function Approximation", legend=:bottomright)

# plot(th5_values_3_sorted, th2_values_3_sorted, label="θ₂ vs θ₅ (real)", lw=2, xlabel="θ₅_reverse (rad)", ylabel="θ₂_reverse (rad)", title="Inverse Function Approximation", legend=:bottomright)



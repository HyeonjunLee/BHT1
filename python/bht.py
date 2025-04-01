import numpy as np
import matplotlib.pyplot as plt

# 물리 상수
k_chrome = 93  # 크롬의 열전도율 (W/m·K)
k_steel = 50   # 철의 열전도율 (W/m·K)
h_inner = 1000  # 내부 공기와의 대류 열전달 계수 (W/m²·K)
h_outer = 25    # 외부 공기와의 대류 열전달 계수 (W/m²·K)
T_inner_air_initial = 2500  # Initial internal air temperature (°C)
T_inner_air_final = 1000    # Final internal air temperature (°C)
T_inner_air = lambda t: T_inner_air_initial + (T_inner_air_final - T_inner_air_initial) * (t * dt / (time_steps * dt))
T_outer_air = 25    # 외부 공기 온도 (°C)

# 포신 치수
r_inner = 155 / 2 / 1000  # 포신 내부 반지름 (m)
r_chrome = r_inner + 0.1 / 1000  # 크롬 코팅층 외부 반지름 (m)
r_outer = r_inner + 77.5 / 1000  # 포신 외부 반지름 (m)

# 계산 설정
dr = 0.0001 # 반지름 방향의 간격 (m)
dt = 0.01    # 시간 간격 (s)
time_steps = 2000  # 시간 스텝 수
r = np.arange(r_inner, r_outer + dr, dr)  # 반지름 배열
T = np.ones_like(r) * T_outer_air  # 초기 온도 배열 (외부 공기 온도로 초기화)

# 반지름에 따른 열전도율 설정
k = np.where(r <= r_chrome, k_chrome, k_steel)

# Stability criteria check and automatic adjustment
Fo = np.max(k) * dt / (dr**2)  # Fourier number
if Fo > 0.25:  # Reduce threshold for stability
    dt = 0.25 * (dr**2) / np.max(k)  # Adjust dt to satisfy stability
    print(f"Adjusted time step dt to {dt:.6f} for stability (Fo={Fo:.2f}).")

# Initialize time array for plotting
time_array = np.arange(0, time_steps * dt, dt)
T_inner_air_values = []

# 내부 공기 물리 상수
c_air = 1005  # 공기의 비열 (J/kg·K)
rho_air = 1.225  # 공기의 밀도 (kg/m³)
V_inner = np.pi * r_inner**2 * 1.0  # 내부 공기의 부피 (m³), 길이를 1m로 가정
m_air = rho_air * V_inner  # 내부 공기의 질량 (kg)
A_inner = 2 * np.pi * r_inner * 1.0  # 내부 표면적 (m²), 길이를 1m로 가정

# 내부 공기 온도 초기화
T_inner_air = T_inner_air_initial

# 시간 루프
for t_step in range(time_steps):
    T_new = T.copy()
    T_surface = T[0]  # 포신 내부 표면 온도
    dT_inner_air_dt = -(h_inner * A_inner / (m_air * c_air)) * (T_inner_air - T_surface)
    T_inner_air += dT_inner_air_dt * dt  # 내부 공기 온도 업데이트
    T_inner_air_values.append(T_inner_air)  # 내부 공기 온도 저장

    for i in range(1, len(r) - 1):
        dT_dr = (T[i+1] - T[i-1]) / (2 * dr)
        d2T_dr2 = (T[i+1] - 2 * T[i] + T[i-1]) / (dr**2)
        # Correct radial term to avoid instability
        radial_term = (1 / r[i]) * (T[i+1] - T[i]) / dr if r[i] > 1e-6 else 0
        T_new[i] = T[i] + dt * (k[i] * (d2T_dr2 + radial_term))

    # Clamp temperature values to avoid overflow
    T_new = np.clip(T_new, T_outer_air, T_inner_air_initial)

    # 내부 경계 조건 (대류)
    T_new[0] = T[0] + dt * h_inner * (T_inner_air - T[0]) / (k[0] / dr)

    # 외부 경계 조건 (대류)
    T_new[-1] = T[-1] + dt * h_outer * (T_outer_air - T[-1]) / (k[-1] / dr)

    # Update temperature array
    T = T_new

# 결과 시각화
plt.figure(figsize=(10, 5))

# Plot temperature distribution
plt.subplot(1, 2, 1)
plt.plot(r * 1000, T, label="Temperature Distribution")
plt.xlabel("Radius (mm)")
plt.ylabel("Temperature (°C)")
plt.title("1D Heat Transfer in Gun Barrel")
plt.legend()
plt.grid()

# Plot internal air temperature over time
plt.subplot(1, 2, 2)
plt.plot(time_array, T_inner_air_values, label="Internal Air Temperature", color="red")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (°C)")
plt.title("Internal Air Temperature Over Time")
plt.legend()
plt.grid()

# Save the plots
plt.tight_layout()
plt.switch_backend('Agg')  # Use a non-interactive backend
plt.savefig('/home/hyeonjun/BHT/python/temperature_distribution.png')  # Save the plot as an image
print("Plot saved as 'temperature_distribution.png'.")
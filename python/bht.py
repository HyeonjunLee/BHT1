import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 포신의 기하학적 특성
inner_radius = 0.155 / 2  # 내부 반지름 (m)
thickness = 0.0775  # 두께 (m)
outer_radius = inner_radius + thickness  # 외부 반지름 (m)
barrel_length = 8.092  # 포신의 길이 (m)

# 물리 상수
k_steel = 50.2  # 철의 열전도율 (W/m·K)
rho_steel = 7850  # 철의 밀도 (kg/m^3)
cp_steel = 470  # 철의 비열 (J/kg·K)
T_air_initial = 2500  # 초기 공기 온도 (K)
T_ambient = 300  # 외부 상온 (K)

# 공기의 물성치
k_air = 0.0257  # 공기의 열전도율 (W/m·K)
dynamic_viscosity = 1.81e-5  # 공기의 동점성 계수 (Pa·s)
Pr = 0.71  # Prandtl 수 (공기)

# 특성 길이와 초기 속도
characteristic_length = inner_radius  # 특성 길이 (m)
initial_velocity = 300  # 초기 공기 속도 (m/s)

# 내부 공기 물리 상수
c_air = 1005  # 공기의 비열 (J/kg·K)
rho_air = 1.225  # 공기의 밀도 (kg/m³)
V_inner = np.pi * inner_radius**2 * barrel_length  # 내부 공기의 부피 (m³)
A_inner = 2 * np.pi * inner_radius * barrel_length  # 내부 표면적 (m²)
m_air = rho_air * V_inner  # 초기 내부 공기의 질량 (kg)

# 시뮬레이션 설정
dr = 0.01  # 반지름 방향의 공간 간격 (m)
dt = 0.01  # 시간 간격 (s)
total_time = 5  # 총 시뮬레이션 시간 (s)

# 열확산율 계산
alpha = k_steel / (rho_steel * cp_steel)  # 철의 열확산율 (m²/s)

# Fourier 안정성 조건에 따른 최대 시간 간격 계산
dt_max = 0.5 * dr**2 / alpha

print(f"Fourier 안정성 조건에 따른 최대 시간 간격: dt_max = {dt_max:.6e} s")

# 현재 설정된 시간 간격과 비교
if dt <= dt_max:
    print("현재 시간 간격(dt)은 안정성 조건을 만족합니다.")
else:
    print("현재 시간 간격(dt)은 안정성 조건을 만족하지 않습니다. dt를 줄이세요.")

dt = min(dt, dt_max)  # 안정성 조건을 만족하도록 dt 조정

# 격자 생성
r = np.arange(inner_radius, outer_radius + dr, dr)
n = len(r)
time_steps = int(total_time / dt)

# 초기 온도 설정
T = np.full(n, T_ambient)
T[0] = T_air_initial  # 내부 공기 온도

# 결과 저장
temperature_history = []

# 내부 공기 온도 기록
internal_air_temperature = []

# 대류 열전달 계수의 시간 의존성
def h_air_dynamic(t):
    if t < 0.1:  # 포탄 발사 직후
        velocity = initial_velocity * np.exp(-t / 0.05)  # 지수적 감소
    else:
        velocity = 0  # 발사 후 공기 흐름이 멈춤

    # Reynolds 수 계산
    Re = (rho_air * velocity * characteristic_length) / dynamic_viscosity

    # Nusselt 수 계산 (난류 조건 가정)
    if Re > 4000:  # 난류 조건
        Nu = 0.023 * Re**0.8 * Pr**0.3
    else:  # 층류 조건
        Nu = 3.66  # 층류 조건의 Nusselt 수

    # 대류 열전달 계수 계산
    h_air = (Nu * k_air) / characteristic_length
    return h_air

# 외부 대류 열전달 계수의 시간 의존성
def h_outer_dynamic(t):
    # 외부 대류 열전달 계수의 시간 의존성을 정의합니다.
    # 여기서는 단순히 상수 값을 반환하도록 설정합니다.
    return 10  # W/m²·K

# 상수 정의
ε = 0.7  # 철의 방사율
σ = 5.67e-8  # 스테판-볼츠만 상수 (W/m²·K⁴)

# 시간 루프
for t_step in range(time_steps):
    T_new = T.copy()
    current_time = t_step * dt

    # 동적으로 변화하는 대류 열전달 계수
    h_air = h_air_dynamic(current_time)
    h_outer = h_outer_dynamic(current_time)

    # 내부 공기와 철 사이의 대류 열전달
    velocity = initial_velocity * max(0, 1 - current_time / 0.1)  # 속도 감소 모델
    mass_flow_rate = rho_air * A_inner * velocity  # 질량 흐름률 (kg/s)
    m_air -= mass_flow_rate * dt  # 시간에 따른 질량 감소
    m_air = max(m_air, 1e-3)  # 최소값 설정

    if m_air > 1e-3:
        T_new[0] = T[0] - h_air * (T[0] - T[1]) * dt / (m_air * c_air)
    else:
        T_new[0] = T[0]  # m_air가 너무 작으면 온도 변화 없음

    # 철 내부의 열전달 (6차 중심 차분법)
    for i in range(3, n - 3):
        d2T_dr2 = (2 * T[i - 3] - 27 * T[i - 2] + 270 * T[i - 1] - 490 * T[i] + 270 * T[i + 1] - 27 * T[i + 2] + 2 * T[i + 3]) / (180 * dr**2)
        dT_dr = (-T[i - 3] + 9 * T[i - 2] - 45 * T[i - 1] + 45 * T[i + 1] - 9 * T[i + 2] + T[i + 3]) / (60 * dr)

        if r[i] > 1e-6:
            T_new[i] = T[i] + dt * (k_steel / (rho_steel * cp_steel)) * (d2T_dr2 + dT_dr / r[i])
        else:
            T_new[i] = T[i] + dt * (k_steel / (rho_steel * cp_steel)) * d2T_dr2

    # 외부 표면에서의 열전달 (대류 + 복사)
    q_radiation = ε * σ * (T[-1]**4 - T_ambient**4)  # 복사 열전달
    T_new[-1] = T[-1] + dt * (
        h_outer * (T_ambient - T[-1]) / (rho_steel * cp_steel * dr) +
        q_radiation / (rho_steel * cp_steel * dr)
    )

    # 디버깅: 주요 변수 출력
    print(f"Time: {current_time:.2f}s, T_new[0]: {T_new[0]}, m_air: {m_air}, h_air: {h_air}")

    # 온도 업데이트
    T = T_new.copy()
    temperature_history.append(T.copy())
    internal_air_temperature.append(T[0])

# 결과 시각화
temperature_history = np.array(temperature_history)
time = np.linspace(0, total_time, time_steps)

plt.figure(figsize=(12, 6))

# 첫 번째 그래프: 반지름에 따른 온도 분포
plt.subplot(1, 2, 1)
for i in range(0, len(temperature_history), time_steps // 10):
    plt.plot(r * 1000, temperature_history[i], label=f'Time = {i * dt:.2f}s')  # 반지름을 mm로 변환
plt.xlabel('Radius (mm)')
plt.ylabel('Temperature (K)')
plt.title('Temperature Distribution in Gun Barrel')
plt.legend()
plt.grid()

# 두 번째 그래프: 시간에 따른 내부 공기 온도 변화
plt.subplot(1, 2, 2)
plt.plot(time, internal_air_temperature, color='red', label='Internal Air Temperature')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.title('Internal Air Temperature Over Time')
plt.legend()
plt.grid()

plt.tight_layout()
plt.show()

# 결과 시각화를 위한 애니메이션 생성
fig, ax = plt.subplots(figsize=(8, 6))
line, = ax.plot([], [], lw=2)
ax.set_xlim(r[0] * 1000, r[-1] * 1000)  # 반지름을 mm로 변환
ax.set_ylim(T_ambient - 50, T_air_initial + 50)  # 온도 범위 설정
ax.set_xlabel('Radius (mm)')
ax.set_ylabel('Temperature (K)')
ax.set_title('Temperature Distribution in Gun Barrel (Animation)')
ax.grid()

# 초기화 함수
def init():
    line.set_data([], [])
    ax.set_title('Temperature Distribution in Gun Barrel (Animation)')  # 초기 제목 설정
    return line,

# 업데이트 함수
def update(frame):
    line.set_data(r * 1000, temperature_history[frame])  # 반지름을 mm로 변환
    current_time = frame * dt  # 현재 시간 계산
    ax.set_title(f'Temperature Distribution at Time = {current_time:.2f} s')  # 제목에 시간 표시
    return line,

# 애니메이션 생성
frames = len(temperature_history)  # 총 프레임 수
interval = 5000 / frames  # 각 프레임 간의 시간 간격 (ms) - 총 재생 시간 5초로 설정
ani = FuncAnimation(fig, update, frames=frames, init_func=init, interval=interval, blit=False)

# 애니메이션 저장 또는 표시
plt.show()
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# 포신의 기하학적 특성
inner_radius = 0.155 / 2  # 내부 반지름 (m)
thickness = 0.0775  # 두께 (m)
outer_radius = inner_radius + thickness  # 외부 반지름 (m)

# 물리 상수
k_steel = 50.2  # 철의 열전도율 (W/m·K)
rho_steel = 7850  # 철의 밀도 (kg/m^3)
cp_steel = 470  # 철의 비열 (J/kg·K)
T_air_initial = 2500  # 초기 공기 온도 (K)
T_ambient = 300  # 외부 상온 (K)

# 내부 공기 물리 상수
c_air = 1005  # 공기의 비열 (J/kg·K)
rho_air = 1.225  # 공기의 밀도 (kg/m³)
V_inner = np.pi * inner_radius**2 * 1.0  # 내부 공기의 부피 (m³), 길이를 1m로 가정
m_air = rho_air * V_inner  # 초기 내부 공기의 질량 (kg)

# 시뮬레이션 설정
dr = 0.001  # 반지름 방향의 공간 간격 (m)
dt = 0.01  # 시간 간격 (s)
total_time = 50  # 총 시뮬레이션 시간 (s)

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
        return 500  # 높은 대류 열전달 계수 (W/m²·K)
    else:
        return 100  # 정상 상태로 감소

# 시간 루프
for t_step in range(time_steps):
    T_new = T.copy()
    current_time = t_step * dt

    # 동적으로 변화하는 대류 열전달 계수
    h_air = h_air_dynamic(current_time)

    # 내부 공기와 철 사이의 대류 열전달
    T_new[0] = T[0] - h_air * (T[0] - T[1]) * dt / (m_air * c_air)

    # 철 내부의 열전달 (유한 차분법)
    for i in range(1, n - 1):
        d2T_dr2 = (T[i + 1] - 2 * T[i] + T[i - 1]) / dr**2
        dT_dr = (T[i + 1] - T[i - 1]) / (2 * dr)
        T_new[i] = T[i] + dt * (k_steel / (rho_steel * cp_steel)) * (d2T_dr2 + dT_dr / r[i])

    # 외부 표면에서의 열전달 (대류)
    T_new[-1] = T[-1] + dt * h_air * (T_ambient - T[-1]) / (rho_steel * cp_steel * dr)

    # 내부 공기의 질량 감소 (포탄 발사로 인한 공기 배출)
    if current_time < 0.1:  # 포탄 발사 직후
        m_air *= 0.99  # 질량 감소율 (예: 1% 감소)

    # 온도 업데이트
    T = T_new.copy()
    temperature_history.append(T.copy())

    # 내부 공기 온도 기록
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
    return line,

# 업데이트 함수
def update(frame):
    line.set_data(r * 1000, temperature_history[frame])  # 반지름을 mm로 변환
    current_time = frame * dt  # 현재 시간 계산
    ax.set_title(f'Temperature Distribution at Time = {current_time:.2f} s')  # 제목에 시간 표시
    return line,

# 애니메이션 생성
ani = FuncAnimation(fig, update, frames=len(temperature_history), init_func=init, blit=True)

# 애니메이션 저장 또는 표시
plt.show()
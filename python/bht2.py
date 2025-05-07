import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# -------------------------
# 1. 시뮬레이션 설정 (도메인)
# -------------------------
r_inner = 0.0775     # 내경 반지름 [m]
r_outer = 0.155 + 0.0775  # 외경 반지름 [m]
L = 8.092            # 포신 길이 [m]

nr = 80              # 반경 방향 노드 수
nz = 400             # 축 방향 노드 수
dr = (r_outer - r_inner) / nr
dz = L / nz
dt = 1e-5            # 시간 간격 [s]
nt = 1000            # 시간 단계 수

# -------------------------
# 2. 재료 특성 (AISI 4340)
# -------------------------
rho = 7850           # kg/m^3
cp = 470             # J/kg-K
k = 44.5             # W/m-K
alpha = k / (rho * cp)

# -------------------------
# 3. 외부 대류 조건
# -------------------------
h = 50.0             # W/m^2-K
T_inf = 300.0        # K

# -------------------------
# 4. 초기 온도장
# -------------------------
T = np.ones((nr, nz)) * T_inf

# -------------------------
# 5. 열유속 (가우시안 함수)
# -------------------------
def heat_flux(t):
    q_max = 1e7        # W/m^2
    t0 = 0.001         # s
    sigma = 0.0003     # s
    return q_max * np.exp(-((t - t0)**2) / (2 * sigma**2))

# -------------------------
# 6. 시뮬레이션 루프
# -------------------------
save_interval = 10
T_record = []

for n in range(nt):
    T_new = T.copy()
    t = n * dt
    for i in range(1, nr - 1):
        for j in range(1, nz - 1):
            r = r_inner + i * dr
            Tr = (T[i+1, j] - 2*T[i, j] + T[i-1, j]) / dr**2
            Tr_add = (1/r) * (T[i+1, j] - T[i-1, j]) / (2 * dr) if r != 0 else 0
            Tz = (T[i, j+1] - 2*T[i, j] + T[i, j-1]) / dz**2
            T_new[i, j] = T[i, j] + alpha * dt * (Tr + Tr_add + Tz)

    # 내부면 (r = r_inner): 열유속 경계조건
    q = heat_flux(t)
    T_new[0, :] += 2 * dt * q / (rho * cp * dr)

    # 외부면 (r = r_outer): 대류 경계조건
    T_new[-1, :] = T[-2, :] + 2 * dr * h / k * (T_inf - T[-1, :])

    # 축 방향 양끝: 단열
    T_new[:, 0] = T_new[:, 1]
    T_new[:, -1] = T_new[:, -2]

    T = T_new.copy()

    if n % save_interval == 0:
        T_record.append(T.copy())

# -------------------------
# 7. 애니메이션 생성
# -------------------------
fig, ax = plt.subplots(figsize=(10, 5))
extent = [0, dz * nz, r_inner, r_outer]
img = ax.imshow(T_record[0], extent=extent, origin='lower', aspect='auto',
                cmap='hot', vmin=300, vmax=np.max(T_record))
plt.colorbar(img, ax=ax, label='Temperature [K]')
ax.set_xlabel('Length (z) [m]')
ax.set_ylabel('Radius (r) [m]')
ax.set_title('Time-Evolving Barrel Temperature')

def update(frame):
    img.set_data(T_record[frame])
    ax.set_title(f'Time = {frame * save_interval * dt * 1000:.2f} ms')
    return [img]

anim = animation.FuncAnimation(fig, update, frames=len(T_record), blit=False, interval=150)

# -------------------------
# 8. 출력
# -------------------------
plt.show()

# # 저장하고 싶을 경우 활성화:
anim.save("barrel_temp.gif", writer="pillow", fps=10)

# 시뮬레이션 루프 끝나고 결과 저장
np.save("T_record.npy", np.array(T_record))
print("✅ T_record 저장 완료: T_record.npy")

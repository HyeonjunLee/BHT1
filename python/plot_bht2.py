import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 저장된 결과 불러오기
T_record = np.load("T_record.npy")  # shape: (frame, nr, nz)

# 도메인 재설정 (시뮬레이션과 동일하게)
r_inner = 0.0775
r_outer = 0.2325
nr = T_record.shape[1]
nz = T_record.shape[2]
dz = 8.092 / nz
r_positions = np.linspace(r_inner, r_outer, nr)
z_index = nz // 2  # 중앙 단면 기준

# 애니메이션 준비
fig, ax = plt.subplots()
line, = ax.plot([], [], lw=2)
ax.set_xlim(r_inner, r_outer)
ax.set_ylim(300, np.max(T_record) + 50)
ax.set_xlabel("Radius [m]")
ax.set_ylabel("Temperature [K]")
ax.set_title("Radial Temperature Profile (z ≈ {:.2f} m)".format(z_index * dz))
time_text = ax.text(0.7, 0.9, '', transform=ax.transAxes)

def init():
    line.set_data([], [])
    time_text.set_text('')
    return line, time_text

def update(frame):
    T_profile = T_record[frame][:, z_index]
    line.set_data(r_positions, T_profile)
    time_text.set_text(f"Time = {frame * 0.0001:.2f} s")  # (adjust if needed)
    return line, time_text

ani = animation.FuncAnimation(fig, update, init_func=init,
                              frames=len(T_record), blit=True, interval=150)
plt.show()
ani.save('temperature_profile_animation.mp4', writer='ffmpeg', fps=30)
plt.show()
ani.save('temperature_profile_animation.gif', writer='pillow', fps=30)
plt.show()
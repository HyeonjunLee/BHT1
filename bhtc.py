import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------- Dittus–Boelter Correlation Function ----------------------
def dittus_boelter_h(V, D, T_air, mode='cooling'):
    """
    Compute the convective heat transfer coefficient (h) for air inside a tube 
    using the Dittus–Boelter correlation.
    
    Parameters:
      V : float
         Flow velocity in m/s.
      D : float
         Characteristic diameter in m (for a circular tube, use inner diameter).
      T_air : float
         Air temperature in Kelvin.
      mode : str
         'heating' (n = 0.4) or 'cooling' (n = 0.3). 여기서는 cooling mode 사용.
    
    Returns:
      h : float
         Convective heat transfer coefficient in W/(m²·K).
    """
    # Air properties (예: 평균 온도 주변)
    rho = 1.2             # Density [kg/m³]
    mu = 1.8e-5           # Dynamic viscosity [Pa·s]
    k_air = 0.026         # Thermal conductivity [W/m·K]
    cp_air = 1000         # Specific heat [J/(kg·K)]
    
    # Reynolds number
    Re = (rho * V * D) / mu
    # Prandtl number
    Pr = (cp_air * mu) / k_air
    # Exponent n: cooling mode → 0.3
    n = 0.3 if mode=='cooling' else 0.4
    # Nusselt number from Dittus–Boelter
    Nu = 0.023 * (Re ** 0.8) * (Pr ** n)
    # h 계산: h = Nu * k_air / D
    h = (Nu * k_air) / D
    return h

# ---------------------- Temperature-dependent Material Properties ----------------------
def thermal_conductivity_chrome(T):
    # Linear approximation for chrome (T in °C)
    # At 20°C: ~93 W/m·K, at 1000°C: ~85 W/m·K
    return 93 - 0.008163 * (T - 20)

def thermal_conductivity_steel(T):
    # Linear approximation for steel (T in °C)
    # At 20°C: ~80 W/m·K, at 1000°C: ~70 W/m·K
    return 80 - 0.010204 * (T - 20)

def specific_heat_chrome(T):
    # Linear approximation for chrome specific heat (J/(kg·K))
    # At 20°C: ~450, at 1000°C: ~500
    return 450 + 0.05102 * (T - 20)

def specific_heat_steel(T):
    # Linear approximation for steel specific heat (J/(kg·K))
    # At 20°C: ~450, at 1000°C: ~520
    return 450 + 0.07143 * (T - 20)

# ---------------------- Composite Simulation Function with Dittus–Boelter ----------------------
def transient_heat_conduction_composite_variable_properties(T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm,
                                                              V, rho_air, c_p_air, T_air_initial,
                                                              total_time, dt, N, num_snapshots):
    """
    Transient heat conduction in a composite cylindrical shell (chrome coating + steel substrate)
    with inner surface convection computed using the Dittus–Boelter correlation.
    Material properties (thermal conductivity and specific heat) are temperature-dependent.
    
    Parameters:
      T_outside: Ambient temperature at the outer wall (°C)
      inner_diameter_mm: Inner diameter (mm)
      steel_thickness_mm: Steel thickness (mm)
      d_chrome_mm: Chrome coating thickness (mm) (e.g., 0.1 mm)
      V: Flow velocity of internal air (m/s)
      rho_air, c_p_air: Density and specific heat of the internal air
      T_air_initial: Initial internal air temperature (°C)
      total_time: Total simulation time (s)
      dt: Time step (s)
      N: Number of radial grid points (for composite region)
      num_snapshots: Number of snapshots to record between 0 and total_time
      
    Returns:
      r: Radial coordinate array for the composite (m), from r1 to r2
      T_record: Recorded temperature distribution in the composite (°C), shape = (num_snapshots, N)
      times: Snapshot times (s)
      T_air_record: Internal air temperature evolution (°C) at snapshot times
      r1: Inner radius (m) (end of air region)
      inner_d: Inner diameter in m (used for Dittus–Boelter)
    """
    # Convert dimensions from mm to m
    inner_d = inner_diameter_mm / 1000.0  # inner diameter in m
    steel_thickness = steel_thickness_mm / 1000.0
    d_chrome = d_chrome_mm / 1000.0

    # Define radii
    r1 = inner_d / 2.0                   # Inner radius (chrome inner surface)
    r_chrome_out = r1 + d_chrome         # Chrome outer surface (steel inner surface)
    r2 = r_chrome_out + steel_thickness  # Outer radius (steel outer surface)

    # Create a uniform radial grid from r1 to r2 (composite region)
    r = np.linspace(r1, r2, N)
    dr = r[1] - r[0]

    # Pre-calculate which nodes belong to chrome (True) or steel (False)
    is_chrome = r <= r_chrome_out

    # Stability: use properties at ambient temperature (보수적 추정)
    k_chrome_init = thermal_conductivity_chrome(T_outside)
    cp_chrome_init = specific_heat_chrome(T_outside)
    alpha_chrome = k_chrome_init / (7190.0 * cp_chrome_init)
    
    k_steel_init = thermal_conductivity_steel(T_outside)
    cp_steel_init = specific_heat_steel(T_outside)
    alpha_steel = k_steel_init / (7850.0 * cp_steel_init)
    
    alpha_max = max(alpha_chrome, alpha_steel)
    dt_max = dr**2 / (2 * alpha_max)
    if dt > dt_max:
        print("Warning: dt is too large for stability. dt_max ≈", dt_max)

    num_steps = int(total_time / dt)
    # Snapshot indices를 linspace로 결정
    snapshot_times = np.linspace(0, total_time, num_snapshots)
    snapshot_indices = np.unique(np.array(np.round(snapshot_times / dt), dtype=int))
    
    # 초기 조건: composite는 전체가 T_outside, 내부 공기는 T_air_initial로 시작
    T = np.ones(N) * T_outside
    T_air = T_air_initial

    T_record = []
    T_air_record = []
    times_record = []

    # Inner surface area (per unit length) at r1
    A_inner = 2 * np.pi * r1
    # Node 0 (chrome inner surface)의 유효 질량: V ≈ A_inner*(dr/2), 밀도: 7190 kg/m³
    m_node0 = 7190.0 * A_inner * (dr / 2)

    # Time integration loop
    for n in range(num_steps + 1):
        if n in snapshot_indices:
            T_record.append(T.copy())
            T_air_record.append(T_air)
            times_record.append(n * dt)
        
        T_new = T.copy()
        
        # ----- Node 0: inner boundary (chrome inner surface) with convection -----
        k0 = thermal_conductivity_chrome(T[0])
        cp0 = specific_heat_chrome(T[0])
        # 인접 노드 1의 재료 결정
        if is_chrome[1]:
            k1 = thermal_conductivity_chrome(T[1])
        else:
            k1 = thermal_conductivity_steel(T[1])
        k_eff = 0.5 * (k0 + k1)
        conduction_flux = k_eff * A_inner * (T[1] - T[0]) / dr
        
        # Dittus–Boelter를 사용하여 현재 내부 공기의 유속 V와 inner_diameter를 바탕으로 h 계산
        # (내부 유체는 cooling mode로 가정)
        h_current = dittus_boelter_h(V, inner_d, T_air + 273.15, mode='cooling')
        convective_flux = h_current * A_inner * (T[0] - T_air)
        dT0_dt = (conduction_flux - convective_flux) / (m_node0 * cp0)
        T_new[0] = T[0] + dt * dT0_dt

        # ----- Internal nodes: i = 1 to N-2 -----
        for i in range(1, N - 1):
            if is_chrome[i]:
                k_i = thermal_conductivity_chrome(T[i])
                cp_i = specific_heat_chrome(T[i])
                rho_i = 7190.0
            else:
                k_i = thermal_conductivity_steel(T[i])
                cp_i = specific_heat_steel(T[i])
                rho_i = 7850.0
            if is_chrome[i - 1]:
                k_im = thermal_conductivity_chrome(T[i - 1])
            else:
                k_im = thermal_conductivity_steel(T[i - 1])
            if is_chrome[i + 1]:
                k_ip = thermal_conductivity_chrome(T[i + 1])
            else:
                k_ip = thermal_conductivity_steel(T[i + 1])
            k_im_eff = 0.5 * (k_im + k_i)
            k_ip_eff = 0.5 * (k_i + k_ip)
            r_im = 0.5 * (r[i] + r[i - 1])
            r_ip = 0.5 * (r[i] + r[i + 1])
            flux_term = (r_ip * k_ip_eff * (T[i + 1] - T[i]) - r_im * k_im_eff * (T[i] - T[i - 1])) / (dr ** 2)
            dT_dt = flux_term / (r[i] * rho_i * cp_i)
            T_new[i] = T[i] + dt * dT_dt

        # ----- Outer boundary: Dirichlet condition (steel outer surface) -----
        T_new[-1] = T_outside

        # ----- Update internal air temperature (lumped system) -----
        m_air = rho_air * np.pi * r1 ** 2
        dT_air_dt = - (h_current * A_inner * (T_air - T[0])) / (m_air * c_p_air)
        T_air = T_air + dt * dT_air_dt

        T = T_new.copy()
    
    return r, np.array(T_record), np.array(times_record), np.array(T_air_record), r1, inner_d

# -------------------- Simulation Parameters --------------------
T_outside = 25                # Ambient temperature (°C)
inner_diameter_mm = 155       # Inner diameter (mm)
steel_thickness_mm = 77.5     # Steel thickness (mm)
d_chrome_mm = 0.1             # Chrome coating thickness (mm)
V = 2.0                     # Internal air flow velocity (m/s)

rho_air = 1.2                 # Air density (kg/m³)
c_p_air = 1000                # Air specific heat (J/(kg·K))
T_air_initial = 1000          # Initial internal air temperature (°C)

total_time = 300              # Total simulation time (s)
dt = 0.01                     # Time step (s)
N = 150                       # Number of radial grid points (for composite)
num_snapshots = 50            # Number of snapshots

# Run simulation
r_shell, T_record, times, T_air_record, r1, inner_d = transient_heat_conduction_composite_variable_properties(
    T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm,
    V, rho_air, c_p_air, T_air_initial,
    total_time, dt, N, num_snapshots
)

# -------------------- Animation 1: Composite Domain (Air + Shell) --------------------
fig1, ax1 = plt.subplots(figsize=(8,6))
# Air region: r from 0 to r1 (n_air points)
n_air = 20
r_air = np.linspace(0, r1, n_air) * 1000  # in mm
r_shell_mm = r_shell * 1000              # in mm

# 초기 composite profile: air 영역은 내부 공기 온도, shell 영역은 T_record[0]
T_air_profile = np.full(n_air, T_air_record[0])
composite_T = np.concatenate([T_air_profile, T_record[0]])
composite_r = np.concatenate([r_air, r_shell_mm])
line1, = ax1.plot(composite_r, composite_T, lw=2)
time_text1 = ax1.text(0.05, 0.90, '', transform=ax1.transAxes, fontsize=12)

ax1.set_xlabel('Radius (mm)', fontsize=14)
ax1.set_ylabel('Temperature (°C)', fontsize=14)
ax1.set_title('Composite Domain (Air + Shell)', fontsize=16)
ax1.grid(True)
ax1.set_xlim(0, r_shell_mm[-1])
ax1.set_ylim(min(T_record.min(), T_air_record.min())-10, max(T_record.max(), T_air_record.max())+10)

def update_composite(frame):
    T_air_profile = np.full(n_air, T_air_record[frame])
    T_shell_profile = T_record[frame]
    composite_T = np.concatenate([T_air_profile, T_shell_profile])
    line1.set_ydata(composite_T)
    time_text1.set_text(f'Time = {times[frame]:.1f} s')
    return line1, time_text1

ani_composite = animation.FuncAnimation(fig1, update_composite, frames=len(times), interval=100, blit=True)

# -------------------- Animation 2: Shell-only Domain (Barrel) --------------------
fig2, ax2 = plt.subplots(figsize=(8,6))
line2, = ax2.plot(r_shell_mm, T_record[0], lw=2)
time_text2 = ax2.text(0.05, 0.90, '', transform=ax2.transAxes, fontsize=12)

ax2.set_xlabel('Radius (mm)', fontsize=14)
ax2.set_ylabel('Temperature (°C)', fontsize=14)
ax2.set_title('Shell-only Domain (Barrel)', fontsize=16)
ax2.grid(True)
ax2.set_xlim(r_shell_mm[0], r_shell_mm[-1])
ax2.set_ylim(min(T_record.min())-10, max(T_record.max())+10)

def update_shell(frame):
    line2.set_ydata(T_record[frame])
    time_text2.set_text(f'Time = {times[frame]:.1f} s')
    return line2, time_text2

ani_shell = animation.FuncAnimation(fig2, update_shell, frames=len(times), interval=100, blit=True)

plt.show()

# -------------------- Plot of Internal Air Temperature Evolution --------------------
fig3, ax3 = plt.subplots(figsize=(8,6))
ax3.plot(times, T_air_record, 'r-o')
ax3.set_xlabel('Time (s)', fontsize=14)
ax3.set_ylabel('Internal Air Temperature (°C)', fontsize=14)
ax3.set_title('Internal Air Temperature Evolution', fontsize=16)
ax3.grid(True)
plt.show()


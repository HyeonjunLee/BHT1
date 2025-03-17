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
         Characteristic diameter in m.
      T_air : float
         Air temperature in Kelvin.
      mode : str
         'cooling' (n = 0.3) or 'heating' (n = 0.4); here we assume cooling.
    
    Returns:
      h in W/(m²·K)
    """
    rho = 1.2       # kg/m³
    mu = 1.8e-5     # Pa·s
    k_air = 0.026   # W/m·K
    cp_air = 1000   # J/(kg·K)
    Re = (rho * V * D) / mu
    Pr = (cp_air * mu) / k_air
    n = 0.3 if mode=='cooling' else 0.4
    Nu = 0.023 * (Re ** 0.8) * (Pr ** n)
    h = (Nu * k_air) / D
    return h

# ---------------------- Temperature-dependent Material Properties ----------------------
def thermal_conductivity_chrome(T):
    # T in °C; linear: 93 W/m·K at 20°C, 85 W/m·K at 1000°C
    return 93 - 0.008163 * (T - 20)

def thermal_conductivity_steel(T):
    # T in °C; linear: 80 W/m·K at 20°C, 70 W/m·K at 1000°C
    return 80 - 0.010204 * (T - 20)

def specific_heat_chrome(T):
    # T in °C; linear: 450 J/(kg·K) at 20°C, 500 J/(kg·K) at 1000°C
    return 450 + 0.05102 * (T - 20)

def specific_heat_steel(T):
    # T in °C; linear: 450 J/(kg·K) at 20°C, 520 J/(kg·K) at 1000°C
    return 450 + 0.07143 * (T - 20)

# ---------------------- Helper: Material properties based on radius ----------------------
def get_material_props(T, r_val, r_chrome_out, include_coating):
    """
    Return (k, rho, cp) for a given temperature T and radial location r_val.
    """
    if include_coating and (r_val <= r_chrome_out):
        return thermal_conductivity_chrome(T), 7190.0, specific_heat_chrome(T)
    else:
        return thermal_conductivity_steel(T), 7850.0, specific_heat_steel(T)

# ---------------------- Grid Generation (Non-uniform) ----------------------
def generate_nonuniform_grid(r1, r_chrome_out, r2, N_coat=10, N_steel=140):
    """
    Generate a non-uniform grid from r1 to r2.
    In the coating region [r1, r_chrome_out], use N_coat points (uniformly spaced).
    In the steel region [r_chrome_out, r2], use N_steel points (uniformly spaced).
    Returns the concatenated grid (with r_chrome_out not duplicated).
    """
    r_coat = np.linspace(r1, r_chrome_out, N_coat)
    r_steel = np.linspace(r_chrome_out, r2, N_steel)
    # Avoid duplicate at the interface
    r_full = np.concatenate([r_coat, r_steel[1:]])
    return r_full

# ---------------------- Transient Simulation with Non-uniform Grid ----------------------
def transient_heat_conduction_nonuniform_with_firing(T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm,
                                                     V, rho_air, c_p_air, T_air_initial,
                                                     total_time, dt, N_coat, N_steel,
                                                     firing_events, T_fire, include_coating=True):
    """
    Transient heat conduction simulation using a non-uniform radial grid.
    The inner boundary has convection with the internal air (using Dittus–Boelter for h).
    Firing events reset the internal air temperature to T_fire.
    
    Parameters:
      T_outside: Ambient temperature at outer wall (°C)
      inner_diameter_mm: Inner diameter (mm)
      steel_thickness_mm: Steel thickness (mm)
      d_chrome_mm: Chrome coating thickness (mm); ignored if include_coating is False.
      V: Internal air flow velocity (m/s)
      rho_air, c_p_air: Density and specific heat of the air.
      T_air_initial: Initial internal air temperature (°C)
      total_time: Total simulation time (s)
      dt: Time step (s)
      N_coat: Number of grid points in coating region.
      N_steel: Number of grid points in steel region.
      firing_events: List of firing event times (s).
      T_fire: Internal air temperature after firing (°C)
      include_coating: If True, include chrome coating.
    
    Returns:
      r: Non-uniform radial grid (m) from r1 to r2.
      T_record: Temperature distribution in shell at snapshot times, shape=(num_snapshots, N).
      times: Snapshot times (s)
      T_air_record: Internal air temperature evolution (°C) at snapshot times.
      r1: Inner radius (m) (start of shell).
      inner_d: Inner diameter (m) (used for convection calculation).
      r_chrome_out: Interface between chrome and steel (m).
    """
    # Geometry in m
    inner_d = inner_diameter_mm / 1000.0
    r1 = inner_d / 2.0
    if include_coating:
        d_chrome = d_chrome_mm / 1000.0
    else:
        d_chrome = 0.0
    r_chrome_out = r1 + d_chrome
    steel_thickness = steel_thickness_mm / 1000.0
    r2 = r_chrome_out + steel_thickness

    # Generate non-uniform grid over shell region
    r = generate_nonuniform_grid(r1, r_chrome_out, r2, N_coat, N_steel)
    N = len(r)
    
    # Stability check (using ambient properties; conservative)
    k_chrome_init, _, cp_chrome_init = thermal_conductivity_chrome(T_outside), None, specific_heat_chrome(T_outside)
    alpha_chrome = k_chrome_init / (7190.0 * cp_chrome_init)
    k_steel_init, _, cp_steel_init = thermal_conductivity_steel(T_outside), None, specific_heat_steel(T_outside)
    alpha_steel = k_steel_init / (7850.0 * cp_steel_init)
    alpha_max = max(alpha_chrome, alpha_steel)
    # For non-uniform grid, use smallest Δr
    dr_min = np.min(np.diff(r))
    dt_max = dr_min**2 / (2 * alpha_max)
    if dt > dt_max:
        print("Warning: dt is too large for stability. dt_max ≈", dt_max)

    num_steps = int(total_time / dt)
    num_snapshots = 50
    snapshot_times = np.linspace(0, total_time, num_snapshots)
    snapshot_indices = np.unique(np.array(np.round(snapshot_times/dt), dtype=int))

    # Initial conditions: shell at T_outside; internal air at T_air_initial.
    T = np.ones(N) * T_outside
    T_air = T_air_initial

    T_record = []
    T_air_record = []
    times_record = []

    # Inner boundary: area and control volume for node 0.
    A_inner = 2 * np.pi * r1  # inner surface area per unit length
    # For node 0, control volume extends from r = r1 to r_face = (r[0]+r[1])/2.
    r_face0 = 0.5 * (r[0] + r[1])
    V0 = np.pi * (r_face0**2 - r1**2)
    # Use material properties at node 0 (if coating exists, chrome; else steel)
    k0, rho0, cp0 = get_material_props(T[0], r[0], r_chrome_out, include_coating)

    # Time integration loop (explicit finite volume)
    for n in range(num_steps + 1):
        current_time = n * dt
        # Firing events: if current time is within dt/2 of an event, reset T_air.
        if any(abs(current_time - t_fire) < dt/2 for t_fire in firing_events):
            T_air = T_fire

        if n in snapshot_indices:
            T_record.append(T.copy())
            T_air_record.append(T_air)
            times_record.append(current_time)

        T_new = T.copy()

        # --- Node 0 (inner boundary) update ---
        # Convection flux at inner face (r = r1)
        h_current = dittus_boelter_h(V, inner_d, T_air + 273.15, mode='cooling')
        F_conv = h_current * A_inner * (T_air - T[0])
        # Conduction flux at outer face of control volume (r_face0)
        # Compute effective conductivity between node 0 and node 1:
        k1_eff, _, _ = get_material_props(T[1], r[1], r_chrome_out, include_coating)
        k0_eff, _, _ = get_material_props(T[0], r[0], r_chrome_out, include_coating)
        k_eff = 0.5 * (k0_eff + k1_eff)
        F_cond = - k_eff * (T[1] - T[0]) / (r[1] - r[0]) * (2 * np.pi * r_face0)
        # Energy balance for node 0:
        T_new[0] = T[0] + dt / (rho0 * cp0 * V0) * (F_conv + F_cond)

        # --- Internal nodes: i = 1 to N-2 ---
        for i in range(1, N-1):
            # Determine interface positions
            r_im = 0.5 * (r[i] + r[i-1])
            r_ip = 0.5 * (r[i+1] + r[i])
            V_i = np.pi * (r_ip**2 - r_im**2)
            # Get material properties at node i
            k_i, rho_i, cp_i = get_material_props(T[i], r[i], r_chrome_out, include_coating)
            # For interfaces, use arithmetic average:
            k_im, _, _ = get_material_props(0.5*(T[i]+T[i-1]), 0.5*(r[i]+r[i-1]), r_chrome_out, include_coating)
            k_ip, _, _ = get_material_props(0.5*(T[i]+T[i+1]), 0.5*(r[i]+r[i+1]), r_chrome_out, include_coating)
            F_im = - k_im * (T[i] - T[i-1]) / (r[i] - r[i-1]) * (2 * np.pi * r_im)
            F_ip = - k_ip * (T[i+1] - T[i]) / (r[i+1] - r[i]) * (2 * np.pi * r_ip)
            T_new[i] = T[i] + dt / (rho_i * cp_i * V_i) * (F_im - F_ip)

        # --- Outer boundary (node N-1): Dirichlet condition ---
        T_new[-1] = T_outside

        # Update internal air temperature (lumped system)
        m_air = rho_air * np.pi * r1**2
        dT_air = - (h_current * A_inner * (T_air - T[0])) / (m_air * c_p_air)
        T_air = T_air + dt * dT_air

        T = T_new.copy()
    
    return r, np.array(T_record), np.array(times_record), np.array(T_air_record), r1, inner_d, r_chrome_out

# -------------------- Simulation Parameters --------------------
T_outside = 25                # Ambient temperature (°C)
inner_diameter_mm = 155       # Inner diameter (mm)
steel_thickness_mm = 77.5     # Steel thickness (mm)
d_chrome_mm = 0.1             # Chrome coating thickness (mm)
V = 2.0                       # Internal air flow velocity (m/s)
rho_air = 1.2                 # Air density (kg/m³)
c_p_air = 1000                # Air specific heat (J/(kg·K))
T_air_initial = 2500          # Initial internal air temperature (°C)
T_fire = 2500               # Temperature immediately after firing (°C)
total_time = 300              # Total simulation time (s)
dt = 0.001                  # Time step (s)
# For non-uniform grid: choose N_coat sufficiently high to resolve 0.1 mm coating
N_coat = 10                 # e.g., 10 points in coating region
N_steel = 140               # e.g., 140 points in steel region

# Firing events (in seconds)
firing_events = [60, 120, 180]

# -------------------- Run Simulation --------------------
r_shell, T_record, times, T_air_record, r1_val, inner_d_val, r_chrome_out = \
    transient_heat_conduction_nonuniform_with_firing(T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm,
                                                     V, rho_air, c_p_air, T_air_initial,
                                                     total_time, dt, N_coat, N_steel,
                                                     firing_events, T_fire, include_coating=True)

# -------------------- Animation: Composite Domain (Air + Shell) --------------------
def create_animation(r_shell, T_record, T_air_record, times, r1, title):
    fig, ax = plt.subplots(figsize=(8,6))
    n_air = 20
    r_air = np.linspace(0, r1, n_air) * 1000  # in mm
    r_shell_mm = r_shell * 1000              # in mm
    T_air_profile = np.full(n_air, T_air_record[0])
    composite_T = np.concatenate([T_air_profile, T_record[0]])
    composite_r = np.concatenate([r_air, r_shell_mm])
    line, = ax.plot(composite_r, composite_T, lw=2)
    time_text = ax.text(0.05, 0.90, '', transform=ax.transAxes, fontsize=12)
    ax.set_xlabel('Radius (mm)', fontsize=14)
    ax.set_ylabel('Temperature (°C)', fontsize=14)
    ax.set_title(title, fontsize=16)
    ax.grid(True)
    ax.set_xlim(0, r_shell_mm[-1])
    ax.set_ylim(np.min(T_record)-10, np.max(T_record)+10)
    
    def update(frame):
        T_air_profile = np.full(n_air, T_air_record[frame])
        composite_T = np.concatenate([T_air_profile, T_record[frame]])
        line.set_ydata(composite_T)
        time_text.set_text(f'Time = {times[frame]:.3f} s')
        return line, time_text

    ani = animation.FuncAnimation(fig, update, frames=len(times), interval=50, blit=True)
    plt.show()

create_animation(r_shell, T_record, T_air_record, times, r1_val, 'Composite Domain (Air + Shell) with Non-uniform Grid')

# -------------------- Plot Internal Air Temperature Evolution --------------------
fig2, ax2 = plt.subplots(figsize=(8,6))
ax2.plot(times, T_air_record, 'r-o')
ax2.set_xlabel('Time (s)', fontsize=14)
ax2.set_ylabel('Internal Air Temperature (°C)', fontsize=14)
ax2.set_title('Internal Air Temperature Evolution', fontsize=16)
ax2.grid(True)
plt.show()


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# ---------------------- Temperature-dependent Material Properties ----------------------
def thermal_conductivity_chrome(T):
    # Linear approximation for chrome (T in °C)
    # At 20°C: ~93 W/m·K, at 1000°C: ~85 W/m·K
    return 93 - 0.008163 * (T - 20)

def thermal_conductivity_steel(T):
    # Linear approximation for steel
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

# ---------------------- Composite Simulation Function ----------------------
def transient_heat_conduction_composite_variable_properties(T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm,
                                                              h_in, rho_air, c_p_air, T_air_initial,
                                                              total_time, dt, N, num_snapshots):
    """
    Transient heat conduction in a composite cylindrical shell (chrome coating + steel substrate)
    with convection on the inner (chrome) surface. Material properties (thermal conductivity and 
    specific heat) are temperature-dependent.
    
    Parameters:
      T_outside: Ambient temperature at the outer wall (°C)
      inner_diameter_mm: Inner diameter (mm)
      steel_thickness_mm: Steel thickness (mm)
      d_chrome_mm: Chrome coating thickness (mm) (e.g., 0.1 mm)
      h_in: Convective heat transfer coefficient at inner surface (W/m²·K)
      rho_air, c_p_air: Density and specific heat of the internal air
      T_air_initial: Initial internal air temperature (°C)
      total_time: Total simulation time (s)
      dt: Time step (s)
      N: Number of radial grid points (across the composite)
      num_snapshots: Number of snapshots to record between 0 and total_time
      
    Returns:
      r: Radial coordinate array for the composite (m), from r1 to r2
      T_record: Recorded temperature distribution in the composite (°C), shape = (num_snapshots, N)
      times: Snapshot times (s)
      T_air_record: Internal air temperature evolution (°C) at snapshot times
    """
    # Convert dimensions from mm to m
    inner_diameter = inner_diameter_mm / 1000.0
    steel_thickness = steel_thickness_mm / 1000.0
    d_chrome = d_chrome_mm / 1000.0

    # Define radii
    r1 = inner_diameter / 2.0             # Inner radius (chrome inner surface)
    r_chrome_out = r1 + d_chrome          # Chrome outer surface (steel inner surface)
    r2 = r_chrome_out + steel_thickness   # Outer radius (steel outer surface)

    # Create a uniform radial grid from r1 to r2 (composite region)
    r = np.linspace(r1, r2, N)
    dr = r[1] - r[0]

    # Pre-calculate which nodes belong to chrome (True) or steel (False)
    is_chrome = r <= r_chrome_out

    # Stability: use properties at ambient temperature (conservative estimate)
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
    # Determine snapshot indices using linspace
    snapshot_times = np.linspace(0, total_time, num_snapshots)
    snapshot_indices = np.unique(np.array(np.round(snapshot_times / dt), dtype=int))
    
    # Initial condition: composite starts at ambient temperature T_outside
    T = np.ones(N) * T_outside
    T_air = T_air_initial

    T_record = []
    T_air_record = []
    times_record = []

    # Inner surface area (per unit length) at r1
    A_inner = 2 * np.pi * r1
    # Effective mass for node 0 (inner cell in chrome): volume ≈ A_inner*(dr/2)
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
        if is_chrome[1]:
            k1 = thermal_conductivity_chrome(T[1])
        else:
            k1 = thermal_conductivity_steel(T[1])
        k_eff = 0.5 * (k0 + k1)
        conduction_flux = k_eff * A_inner * (T[1] - T[0]) / dr
        convective_flux = h_in * A_inner * (T[0] - T_air)
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
        dT_air_dt = - (h_in * A_inner * (T_air - T[0])) / (m_air * c_p_air)
        T_air = T_air + dt * dT_air_dt

        T = T_new.copy()
    
    return r, np.array(T_record), np.array(times_record), np.array(T_air_record), r1

# -------------------- Simulation Parameters --------------------
T_outside = 25                # Ambient temperature (°C)
inner_diameter_mm = 155       # Inner diameter (mm)
steel_thickness_mm = 77.5     # Steel thickness (mm)
d_chrome_mm = 0.1             # Chrome coating thickness (mm)
h_in = 25                     # Convective heat transfer coefficient at inner surface (W/m²·K)

rho_air = 1.2                 # Air density (kg/m³)
c_p_air = 1000                # Air specific heat (J/(kg·K))
T_air_initial = 1000          # Initial internal air temperature (°C)

total_time = 15              # Total simulation time (s)
dt = 0.001                     # Time step (s)
N = 200                       # Number of radial grid points (for composite)
num_snapshots = 50            # Number of snapshots

# Run simulation
r_shell, T_record, times, T_air_record, r1 = transient_heat_conduction_composite_variable_properties(
    T_outside, inner_diameter_mm, steel_thickness_mm, d_chrome_mm,
    h_in, rho_air, c_p_air, T_air_initial,
    total_time, dt, N, num_snapshots
)

# -------------------- Animation Including Air and Shell Regions --------------------
fig, ax = plt.subplots(figsize=(8,6))
# Define air region: from r = 0 to r = r1, using n_air points
n_air = 20
r_air = np.linspace(0, r1, n_air) * 1000  # in mm
r_shell_mm = r_shell * 1000  # already in mm

# Prepare initial composite profile:
# For air region: constant T_air (from first snapshot), for shell: T_record[0]
T_air_initial_profile = np.full(n_air, T_air_record[0])
composite_r = np.concatenate([r_air, r_shell_mm])
composite_T = np.concatenate([T_air_initial_profile, T_record[0]])

line, = ax.plot(composite_r, composite_T, lw=2)
time_text = ax.text(0.05, 0.90, '', transform=ax.transAxes, fontsize=12)

ax.set_xlabel('Radius (mm)', fontsize=14)
ax.set_ylabel('Temperature (°C)', fontsize=14)
ax.set_title('Transient Temperature Distribution (Air + Barrel)', fontsize=16)
ax.grid(True)
ax.set_xlim(0, r_shell_mm[-1])
ax.set_ylim(min(T_record.min(), min(T_air_record))-10, max(T_record.max(), max(T_air_record))+10)

def update(frame):
    # For current snapshot: get T_air and shell profile
    T_air_current = T_air_record[frame]
    T_air_profile = np.full(n_air, T_air_current)
    T_shell_profile = T_record[frame]
    composite_T = np.concatenate([T_air_profile, T_shell_profile])
    line.set_ydata(composite_T)
    time_text.set_text(f'Time = {times[frame]:.1f} s')
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=len(times), interval=100, blit=True)
plt.show()

# -------------------- Plot of Internal Air Temperature Evolution --------------------
fig2, ax2 = plt.subplots(figsize=(8,6))
ax2.plot(times, T_air_record, 'r-o')
ax2.set_xlabel('Time (s)', fontsize=14)
ax2.set_ylabel('Internal Air Temperature (°C)', fontsize=14)
ax2.set_title('Internal Air Temperature Evolution', fontsize=16)
ax2.grid(True)
plt.show()


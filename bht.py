import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# Temperature-dependent thermal conductivity (W/m·K)
def thermal_conductivity(T):
    # Linear approximation: at 20°C, k = 80; at 1000°C, k ≈ 60 W/m·K
    return 80 - 0.0204 * (T - 20)

# Temperature-dependent convective heat transfer coefficient (W/m²·K)
def convective_heat_transfer_coefficient(T_surface, T_air):
    # Simple model: base value 25 W/m²·K plus a term proportional to the temperature difference
    return 25 + 0.01 * abs(T_surface - T_air)

def transient_heat_conduction_cylindrical_convective_inner_variable(T_outside, inner_diameter_mm, thickness_mm,
                                                                    rho, c_p, rho_air, c_p_air,
                                                                    T_air_initial,
                                                                    total_time, dt, N):
    """
    Transient heat conduction in a cylindrical shell with variable thermal conductivity and convective heat transfer coefficient.
    The inner wall exchanges energy with internal air via convection and conduction, while the outer wall is fixed at T_outside.
    
    Parameters:
      T_outside: Ambient temperature at the outer wall (°C)
      inner_diameter_mm: Inner wall diameter (mm)
      thickness_mm: Shell thickness (mm)
      rho: Density of the shell material (kg/m³)
      c_p: Specific heat capacity of the shell material (J/(kg·K))
      rho_air: Density of the internal air (kg/m³)
      c_p_air: Specific heat capacity of the internal air (J/(kg·K))
      T_air_initial: Initial temperature of the internal air (°C)
      total_time: Total simulation time (s)
      dt: Time step (s)
      N: Number of grid points in the radial direction
      
    Returns:
      r: Radial coordinate array (m)
      T_record: Recorded temperature distribution in the shell (°C), shape = (num_snapshots, N)
      times: Snapshot times (s)
      T_air_record: Internal air temperature evolution (°C)
    """
    # Convert mm to m
    inner_diameter = inner_diameter_mm / 1000.0
    thickness = thickness_mm / 1000.0

    r1 = inner_diameter / 2.0  # inner radius
    r2 = r1 + thickness        # outer radius

    # Create uniform radial grid from r1 to r2
    r = np.linspace(r1, r2, N)
    dr = r[1] - r[0]

    # Initial conditions
    T = np.ones(N) * T_outside
    T_air = T_air_initial

    # Stability estimate (using approximate maximum k)
    k_max = thermal_conductivity(T_air_initial)
    alpha_max = k_max / (rho * c_p)
    dt_max = dr**2 / (2 * alpha_max)
    if dt > dt_max:
        print("Warning: dt is too large for stability. dt_max ≈", dt_max)

    num_steps = int(total_time / dt)
    
    # Record snapshots at selected times (0, 10, 50, 100, total_time)
    #snapshot_times = [0, 10, 50, 100, total_time]
    snapshot_times = np.linspace(0,total_time,100)
    snapshot_indices = [int(t/dt) for t in snapshot_times]
    T_record = []
    T_air_record = []
    times_record = []

    # Inner surface area (per unit length)
    A_inner = 2 * np.pi * r1
    # Effective mass of node 0 (using half-cell volume for unit length): m_node = rho * A_inner * (dr/2) = rho * π * r1 * dr
    m_node = rho * np.pi * r1 * dr

    for n in range(num_steps + 1):
        if n in snapshot_indices:
            T_record.append(T.copy())
            T_air_record.append(T_air)
            times_record.append(n * dt)
            
        T_new = T.copy()

        # --- Update for node 0 (inner wall) with convection and conduction ---
        # Effective conductivity at inner interface (average of node 0 and 1)
        k_inner = (thermal_conductivity(T[0]) + thermal_conductivity(T[1])) / 2
        conduction_term = k_inner * A_inner * (T[1] - T[0]) / dr
        # Convective heat loss from inner surface
        h_current = convective_heat_transfer_coefficient(T[0], T_air)
        convective_term = h_current * A_inner * (T[0] - T_air)
        dT0_dt = (conduction_term - convective_term) / (m_node * c_p)
        T_new[0] = T[0] + dt * dT0_dt

        # --- Update for internal nodes (i = 1 to N-2) with variable conductivity ---
        for i in range(1, N-1):
            # Compute effective conductivity at interfaces using arithmetic average
            k_ip = (thermal_conductivity(T[i]) + thermal_conductivity(T[i+1])) / 2
            k_im = (thermal_conductivity(T[i]) + thermal_conductivity(T[i-1])) / 2
            # Compute effective radial positions for interfaces
            r_ip = (r[i] + r[i+1]) / 2
            r_im = (r[i] + r[i-1]) / 2
            dT_dr_ip = (T[i+1] - T[i]) / dr
            dT_dr_im = (T[i] - T[i-1]) / dr
            dT_dt = (1 / (rho * c_p)) * (1/(r[i]*dr)) * ( r_ip * k_ip * dT_dr_ip - r_im * k_im * dT_dr_im )
            T_new[i] = T[i] + dt * dT_dt

        # Outer wall (node N-1): fixed at ambient temperature
        T_new[-1] = T_outside

        # --- Update for internal air temperature (lumped system) ---
        m_air = rho_air * np.pi * r1**2  # mass per unit length of air
        dT_air_dt = - (convective_heat_transfer_coefficient(T[0], T_air) * A_inner * (T_air - T[0])) / (m_air * c_p_air)
        T_air = T_air + dt * dT_air_dt

        T = T_new.copy()

    return r, np.array(T_record), np.array(times_record), np.array(T_air_record)

# --------------------- Simulation Parameters ---------------------
T_outside = 25            # Ambient temperature (°C)
inner_diameter_mm = 155   # Inner wall diameter (mm)
thickness_mm = 77.5       # Shell thickness (mm)
rho = 7850                # Density of steel (kg/m³)
c_p = 450                 # Specific heat capacity of steel (J/(kg·K))
rho_air = 1.2             # Density of air (kg/m³)
c_p_air = 1000            # Specific heat capacity of air (J/(kg·K))
T_air_initial = 3000      # Initial internal air temperature (°C)
total_time = 15           # Total simulation time (s)
dt = 0.01                 # Time step (s)
N = 100                   # Number of radial grid points

# Run the simulation
r, T_record, times, T_air_record = transient_heat_conduction_cylindrical_convective_inner_variable(
    T_outside, inner_diameter_mm, thickness_mm,
    rho, c_p, rho_air, c_p_air,
    T_air_initial,
    total_time, dt, N
)

# --------------------- Animation of Temperature Distribution ---------------------
fig, ax = plt.subplots(figsize=(8,6))
r_mm = r * 1000  # Convert radius from m to mm for plotting
line, = ax.plot(r_mm, T_record[0], lw=2)
time_text = ax.text(0.05, 0.90, '', transform=ax.transAxes, fontsize=12)

ax.set_xlabel('Radius (mm)', fontsize=14)
ax.set_ylabel('Temperature (°C)', fontsize=14)
ax.set_title('Transient Temperature Distribution in the Cylinder', fontsize=16)
ax.grid(True)
ax.set_xlim(r_mm[0], r_mm[-1])
ax.set_ylim(np.min(T_record)-10, np.max(T_record)+10)

def update(frame):
    line.set_ydata(T_record[frame])
    time_text.set_text(f'Time = {times[frame]:.1f} s')
    return line, time_text

ani = animation.FuncAnimation(fig, update, frames=len(times), interval=100, blit=True)
plt.show()

# --------------------- Plot of Internal Air Temperature Evolution ---------------------
fig2, ax2 = plt.subplots(figsize=(8,6))
ax2.plot(times, T_air_record, 'r-o')
ax2.set_xlabel('Time (s)', fontsize=14)
ax2.set_ylabel('Internal Air Temperature (°C)', fontsize=14)
ax2.set_title('Internal Air Temperature Evolution', fontsize=16)
ax2.grid(True)
plt.show()

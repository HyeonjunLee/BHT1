import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def transient_heat_conduction_cylindrical_convective_inner(T_outside, inner_diameter_mm, thickness_mm,
                                                            k, rho, c_p,
                                                            h_in, rho_air, c_p_air,
                                                            T_air_initial,
                                                            total_time, dt, N):
    """
    Transient heat conduction in a cylindrical shell with convective cooling on the inner surface.
    The inner wall exchanges energy with high-temperature air via convection while the outer wall is fixed at T_outside.
    
    Parameters:
      T_outside: Ambient temperature at the outer wall (°C)
      inner_diameter_mm: Inner wall diameter (mm)
      thickness_mm: Shell thickness (mm)
      k: Thermal conductivity of the shell material (W/m·K)
      rho: Density of the shell material (kg/m³)
      c_p: Specific heat capacity of the shell material (J/(kg·K))
      h_in: Convective heat transfer coefficient at the inner wall (W/m²·K)
      rho_air: Density of the internal air (kg/m³)
      c_p_air: Specific heat capacity of the air (J/(kg·K))
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
    
    r1 = inner_diameter / 2.0  # Inner radius
    r2 = r1 + thickness        # Outer radius
    
    alpha = k / (rho * c_p)    # Thermal diffusivity
    
    # Create uniform radial grid from r1 to r2
    r = np.linspace(r1, r2, N)
    dr = r[1] - r[0]
    
    # Initial conditions: the shell is at ambient temperature, while the internal air is at T_air_initial
    T = np.ones(N) * T_outside
    T_air = T_air_initial
    
    # Stability condition for explicit scheme: dt <= dr^2/(2*alpha)
    dt_max = dr**2 / (2 * alpha)
    if dt > dt_max:
        print("Warning: dt is too large for stability. dt_max ≈", dt_max)
    
    num_steps = int(total_time / dt)
    
    # Record snapshots at selected times
    #snapshot_times = [0, 10, 50, 100, total_time]
    snapshot_times = np.linspace(0,total_time,200)
    snapshot_indices = [int(t/dt) for t in snapshot_times]
    T_record = []
    T_air_record = []
    times_record = []
    
    # Inner surface area (per unit length)
    A_inner = 2 * np.pi * r1
    
    # Time integration loop
    for n in range(num_steps + 1):
        if n in snapshot_indices:
            T_record.append(T.copy())
            T_air_record.append(T_air)
            times_record.append(n * dt)
        
        T_new = T.copy()
        
        # Update node 0 (inner wall) using energy balance with conduction and convective exchange:
        # m_node * c_p * dT0/dt = k*A_inner*(T[1]-T[0])/dr - h_in*A_inner*(T[0]-T_air)
        # Effective mass of node 0 (using half-cell volume for unit length): m_node = rho * A_inner * (dr/2) = rho * π * r1 * dr
        m_node = rho * np.pi * r1 * dr
        dT0_dt = (k * A_inner * (T[1] - T[0]) / dr - h_in * A_inner * (T[0] - T_air)) / (m_node * c_p)
        T_new[0] = T[0] + dt * dT0_dt
        
        # Update internal nodes (i = 1 to N-2) using the finite-difference approximation in cylindrical coordinates
        for i in range(1, N-1):
            T_new[i] = T[i] + dt * alpha * ( (r[i+1]*(T[i+1]-T[i]) - r[i-1]*(T[i]-T[i-1]) ) / (r[i]*dr**2) )
        
        # Outer wall (node N-1): fixed at ambient temperature
        T_new[-1] = T_outside
        
        # Update the internal air temperature (lumped system assumption)
        # Air mass per unit length: m_air = rho_air * π * r1²
        # Energy balance: m_air*c_p_air*dT_air/dt = - h_in*A_inner*(T_air - T[0])
        m_air = rho_air * np.pi * r1**2
        dT_air_dt = - (h_in * A_inner * (T_air - T[0])) / (m_air * c_p_air)
        T_air = T_air + dt * dT_air_dt
        
        T = T_new.copy()
        
    return r, np.array(T_record), np.array(times_record), np.array(T_air_record)

# --------------------- Simulation Parameters ---------------------
T_outside = 25            # Ambient temperature (°C)
inner_diameter_mm = 155   # Inner wall diameter (mm)
thickness_mm = 77.5       # Shell thickness (mm)
k = 80                    # Thermal conductivity of steel (W/m·K)
rho = 7850                # Density of steel (kg/m³)
c_p = 450                 # Specific heat of steel (J/(kg·K))
h_in = 25                 # Convective heat transfer coefficient at inner wall (W/m²·K)

rho_air = 1.2             # Density of air (kg/m³)
c_p_air = 1000            # Specific heat of air (J/(kg·K))
T_air_initial = 2500      # Initial internal air temperature (°C)

total_time = 10          # Total simulation time (s)
dt = 0.01                 # Time step (s)
N = 100                    # Number of grid points in the radial direction

# Run the simulation
r, T_record, times, T_air_record = transient_heat_conduction_cylindrical_convective_inner(
    T_outside, inner_diameter_mm, thickness_mm,
    k, rho, c_p, h_in, rho_air, c_p_air,
    T_air_initial,
    total_time, dt, N
)

# --------------------- Animation of Temperature Distribution ---------------------
fig, ax = plt.subplots(figsize=(8,6))
# Convert r from m to mm for plotting
r_mm = r * 1000
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

# --------------------- Optional: Plot of Internal Air Temperature ---------------------
fig2, ax2 = plt.subplots(figsize=(8,6))
ax2.plot(times, T_air_record, 'r-o')
ax2.set_xlabel('Time (s)', fontsize=14)
ax2.set_ylabel('Internal Air Temperature (°C)', fontsize=14)
ax2.set_title('Internal Air Temperature Evolution', fontsize=16)
ax2.grid(True)
plt.show()


import numpy as np
import matplotlib.pyplot as plt

def read_output_file(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()

    steps = []
    times = []
    T_air = []
    T_profiles = []

    for line in lines:
        if line.startswith('Step:'):
            steps.append(int(line.split()[-1]))
        elif line.startswith('Time:'):
            times.append(float(line.split()[-1]))
        elif line.startswith('T_air:'):
            T_air.append(float(line.split()[-1]))
        elif line.startswith('T:'):
            T_profile = list(map(float, line.split()[1:]))
            T_profiles.append(T_profile)

    return np.array(times), np.array(T_air), np.array(T_profiles)

def plot_results(times, T_air, T_profiles):
    plt.figure(figsize=(12, 6))

    # Plot internal air temperature over time
    plt.subplot(1, 2, 1)
    plt.plot(times, T_air, label='Internal Air Temperature')
    plt.xlabel('Time (s)')
    plt.ylabel('Temperature (°C)')
    plt.title('Internal Air Temperature Over Time')
    plt.legend()

    # Plot temperature profiles at different times
    plt.subplot(1, 2, 2)
    for i in range(0, len(times), max(1, len(times)//10)):
        plt.plot(T_profiles[i], label=f'Time = {times[i]:.1f} s')
    plt.xlabel('Radial Position Index')
    plt.ylabel('Temperature (°C)')
    plt.title('Temperature Profiles in Shell')
    plt.legend()

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    output_file = 'output.txt'
    times, T_air, T_profiles = read_output_file(output_file)
    plot_results(times, T_air, T_profiles)
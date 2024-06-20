import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

#constants
g = 9.81  # acceleration due to gravity, in m/s^2

#single pendulum simulation first
def simple_pendulum(t, y):
    theta, omega = y
    dydt = [omega, -g / L * np.sin(theta)]
    return dydt

def simulate_simple_pendulum():
    global L
    L = 1.0   # length of the pendulum, in meters
    theta0 = np.pi / 4  # initial angle (in radians)
    omega0 = 0.0  # initial angular velocity

    t_span = (0, 10)
    t_eval = np.linspace(0, 10, 1000)
    y0 = [theta0, omega0]

    sol = solve_ivp(simple_pendulum, t_span, y0, t_eval=t_eval)

    plt.plot(sol.t, sol.y[0], label='Theta (Angular Displacement)')
    plt.plot(sol.t, sol.y[1], label='Omega (Angular Velocity)')
    plt.xlabel('Time (s)')
    plt.ylabel('Value')
    plt.legend()
    plt.title('Simple Pendulum Simulation')
    plt.show()
    
    if __name__ == "__main__":
        simulate_simple_pendulum()
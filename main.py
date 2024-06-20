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
    
# Double Pendulum Simulation
def double_pendulum(t, y):
    theta1, omega1, theta2, omega2 = y
    delta = theta2 - theta1

    denom1 = (2 * m1 + m2 - m2 * np.cos(2 * theta1 - 2 * theta2))
    denom2 = (L2 / L1) * denom1

    dtheta1_dt = omega1
    dtheta2_dt = omega2

    domega1_dt = (m2 * g * np.sin(theta2) * np.cos(delta) - m2 * np.sin(delta) * (L2 * omega2 ** 2 + L1 * omega1 ** 2 * np.cos(delta)) - (m1 + m2) * g * np.sin(theta1)) / (L1 * denom1)

    domega2_dt = ((m1 + m2) * (L1 * omega1 ** 2 * np.sin(delta) - g * np.sin(theta2) + g * np.sin(theta1) * np.cos(delta)) + m2 * L2 * omega2 ** 2 * np.sin(delta) * np.cos(delta)) / (L2 * denom2)

    return [dtheta1_dt, domega1_dt, dtheta2_dt, domega2_dt]
    
    if __name__ == "__main__":
        simulate_simple_pendulum()
        
        simulate_double_pendulum()
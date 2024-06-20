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
    
def simulate_double_pendulum():
    global m1, m2, L1, L2, sol_double
    m1 = 1.0  # mass of the first pendulum
    m2 = 1.0  # mass of the second pendulum
    L1 = 1.0  # length of the first pendulum
    L2 = 1.0  # length of the second pendulum
    theta1_0 = np.pi / 2  # initial angle of the first pendulum
    theta2_0 = np.pi / 2  # initial angle of the second pendulum
    omega1_0 = 0.0  # initial angular velocity of the first pendulum
    omega2_0 = 0.0  # initial angular velocity of the second pendulum

    t_span = (0, 10)
    t_eval = np.linspace(0, 10, 1000)
    y0 = [theta1_0, omega1_0, theta2_0, omega2_0]

    sol_double = solve_ivp(double_pendulum, t_span, y0, t_eval=t_eval)

    plt.plot(sol_double.t, sol_double.y[0], label='Theta1 (Angular Displacement of Pendulum 1)')
    plt.plot(sol_double.t, sol_double.y[1], label='Omega1 (Angular Velocity of Pendulum 1)')
    plt.plot(sol_double.t, sol_double.y[2], label='Theta2 (Angular Displacement of Pendulum 2)')
    plt.plot(sol_double.t, sol_double.y[3], label='Omega2 (Angular Velocity of Pendulum 2)')
    plt.xlabel('Time (s)')
    plt.ylabel('Value')
    plt.legend()
    plt.title('Double Pendulum Simulation')
    plt.show()

    if __name__ == "__main__":
        simulate_simple_pendulum()
        
        simulate_double_pendulum()
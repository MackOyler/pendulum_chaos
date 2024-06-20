import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
import matplotlib.animation as animation

#constants
g = 9.81 #acceleration due to gravity, in m/s^2

#single pendulum simulation first
def simple_pendulu(t, y):
    theta, omega = y
    dydt = [omega, -g / L * np.sin(theta)]
    return dydt
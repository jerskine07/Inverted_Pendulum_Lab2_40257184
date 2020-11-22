""""
John Erskine
40257184
jerskine07@qub.ac.uk
Inverted Pendulum Lab 2
Question 4
"""

import numpy as np
import control as ctrl
import matplotlib.pyplot as plt

### Note: This Question took a while to run for me, please be patient ###

# New Libray
import sympy as sm
from sympy import inverse_laplace_transform as ilt

# Code from Question 1 starting off then laplace is applied

# Define all involved symbolic variables
M, m, ell, g, x1, x2, x3, x4, F = \
    sm.symbols('M, m, ell, g, x1, x2, x3, x4, F')

# Define phi
phi = 4*m*ell*x4**2 * sm.sin(x3) + 4*F - 3*m*g*sm.sin(x3)*sm.cos(x3)
phi /= 4*(M+m) - 3*m*sm.cos(x3)**2

# define psi
psi = -3*(4*m*ell*x4**2 * sm.sin(x3) * sm.cos(x3) + F*sm.cos(x3) - (M + m)*g*sm.sin(x3))
psi /= (4*(M + m) - 3*m*sm.cos(x3)**2)* ell

# differentiation of phi
dphi_F = phi.diff(F)
dphi_x3 = phi.diff(x3)
dphi_x4 = phi.diff(x4)

# differentiation of psi
dpsi_F = psi.diff(F)
dpsi_x3 = psi.diff(x3)
dpsi_x4 = psi.diff(x4)

# substitution of values
def evaluate_at_equilibrium(f):
    return f.subs(([F,0], [x3,0], [x4,0]))

# Simplification using substitued values
dphi_F_eq = evaluate_at_equilibrium(phi.diff(F))
dphi_x3_eq = evaluate_at_equilibrium(phi.diff(x3))
dphi_x4_eq = evaluate_at_equilibrium(phi.diff(x4))

dpsi_F_eq = evaluate_at_equilibrium(psi.diff(F))
dpsi_x3_eq = evaluate_at_equilibrium(psi.diff(x3))
dpsi_x4_eq = evaluate_at_equilibrium(psi.diff(x4))

# Start of original question3 code
# Define symbols s, t, and w for laplace transformations
s, t = sm.symbols('s, t')
w = sm.symbols('w', real=True)

# x2' = aF - bx3
a = dphi_F_eq
b = -dphi_x3_eq

# x4' = -cF + dx3
c = -dpsi_F_eq
d = dphi_x3_eq

# Define the values for constants M, m, ell, and g
M_value = 0.3
m_value = 0.1
ell_value = 0.35
g_value = 9.81

# Define the values for a, b, c, and d
a_value = float(a.subs([(M, M_value), (m, m_value)]))
b_value = float(b.subs([(M, M_value), (m, m_value), (g, g_value)]))
c_value = float(c.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)]))
d_value = float(d.subs([(M, M_value), (m, m_value), (g, g_value), (ell, ell_value)]))

# Define symbols to store the partial derivatives as coefficients
a, b, c, d = sm.symbols('a:d', real=True, positive=True)

# Define the transfer functions G_theta and G_x
G_theta = c / (d - s ** 2)
G_x = (a * d - a * s ** 2 - b * c) / (d * s ** 2 - s ** 4)

# Declare variables for helping to draw the graph
num_points = 1000  # The resolution of the graph
dt = 0.2  # Time t ranges between 0 and 0.2 seconds
t_span = np.linspace(0, dt, num_points)

# Define the input signal for the System
input_signal = np.sin(100 * t_span ** 2)

# Determine the response of the System
G_theta_t_out, G_theta_y_out, G_theta_x_out = ctrl.forced_response(G_theta, t_span, input_signal)
G_x_t_out, G_x_y_out, G_x_x_out = ctrl.forced_response(G_x, t_span, input_signal)

# Plot angle against time using the results from G_theta
plt.plot(G_theta_t_out, G_theta_y_out)
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('Angle (rad)')
plt.savefig('figures\\question4angle.svg', format='svg')  # Save the graph as a .svg file
plt.show()

# Plot x position against time using the results from G_x
plt.plot(G_x_t_out, G_x_y_out)
plt.grid()
plt.xlabel('Time (s)')
plt.ylabel('x position (m)')
plt.savefig('figures\\question4x.svg', format='svg')  # Save the graph as a .svg file
plt.show()

""""
John Erskine
40257184
jerskine07@qub.ac.uk
Inverted Pendulum Lab 2
Question 3
"""
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

# Define F in the s domain for impulse, step, and frequency response
F_s_impulse = 1
F_s_step = 1 / s
F_s_frequency = w / (s ** 2 + w ** 2)

# Define the impulse, step, and frequency responses for X3_s and X1_s (s domain)
X3_s_impulse = G_theta * F_s_impulse
X3_s_step = G_theta * F_s_step
X3_s_frequency = G_theta * F_s_frequency
X1_s_impulse = G_x * F_s_impulse
X1_s_step = G_x * F_s_step
X1_s_frequency = G_x * F_s_frequency

# Define the impulse, step, and frequency responses for x3_t and x1_t (t domain)
x3_t_impulse = ilt(X3_s_impulse, s, t).simplify()
x3_t_step = ilt(X3_s_step, s, t).simplify()
x3_t_frequency = ilt(X3_s_frequency, s, t, w).simplify()
x1_t_impulse = ilt(X1_s_impulse, s, t).simplify()
x1_t_step = ilt(X1_s_step, s, t).simplify()
x1_t_frequency = ilt(X1_s_frequency, s, t, w).simplify()

# Pretty print all of the calculated responses in the t domain
print('\n x3 t impulse')
sm.pprint(x3_t_impulse)
print('\n x3 t step')
sm.pprint(x3_t_step)
print('\n x3 t freq')
sm.pprint(x3_t_frequency)
print('\n x1 t impulse')
sm.pprint(x1_t_impulse)
print('\n x1 t step')
sm.pprint(x1_t_step)
print('\n x3 t freq')
sm.pprint(x1_t_frequency)

# Print the equations again in LaTeX format
print('\n Latex form: x3 t impulse')
print(sm.latex(x3_t_impulse))
print('\n Latex form: x3 t step')
print(sm.latex(x3_t_step))
print('\n Latex form: x3 t freq')
print(sm.latex(x3_t_frequency))
print('\n Latex form: x1 t impulse')
print(sm.latex(x1_t_impulse))
print('\n Latex form: x1 t step')
print(sm.latex(x1_t_step))
print('\n Latex form: x1 t freq')
print(sm.latex(x1_t_frequency))


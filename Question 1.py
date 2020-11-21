""""
John Erskine
40257184
jerskine07@qub.ac.uk
Inverted Pendulum Lab 2
Question 1
"""
#import libraries
import sympy as sm

# Define all involved symbolic variables
M, m, ell, g, x1, x2, x3, x4, F = \
    sm.symbols('M, m, ell, g, x1, x2, x3, x4, F')

#Define phi
phi = 4*m*ell*x4**2 * sm.sin(x3) + 4*F - 3*m*g*sm.sin(x3)*sm.cos(x3)
phi /= 4*(M+m) - 3*m*sm.cos(x3)**2

#define psi
psi = -3*(4*m*ell*x4**2 * sm.sin(x3) * sm.cos(x3) + F*sm.cos(x3) - (M + m)*g*sm.sin(x3))
psi /= (4*(M + m) - 3*m*sm.cos(x3)**2)* ell

#differentiation of phi
dphi_F = phi.diff(F)
dphi_x3 = phi.diff(x3)
dphi_x4 = phi.diff(x4)

#differentiation of psi
dpsi_F = psi.diff(F)
dpsi_x3 = psi.diff(x3)
dpsi_x4 = psi.diff(x4)

#substitution of values
def evaluate_at_equilibrium(f):
    return f.subs(([F,0], [x3,0], [x4,0]))

#Simplification using substitued values
dphi_F_eq = evaluate_at_equilibrium(phi.diff(F))
dphi_x3_eq = evaluate_at_equilibrium(phi.diff(x3))
dphi_x4_eq = evaluate_at_equilibrium(phi.diff(x4))

dpsi_F_eq = evaluate_at_equilibrium(psi.diff(F))
dpsi_x3_eq = evaluate_at_equilibrium(psi.diff(x3))
dpsi_x4_eq = evaluate_at_equilibrium(psi.diff(x4))

#Printing of results
print('Phi Partial differentiation 1 =')
sm.pprint(dphi_F_eq)
print('Phi Partial differentiation 2 =')
sm.pprint(dphi_x3_eq)
print('Phi Partial differentiation 3 =')
sm.pprint(dphi_x4_eq)

print('Psi Partial differentiation 1 =')
sm.pprint(dpsi_F_eq)
print('Psi Partial differentiation 2 =')
sm.pprint(dpsi_x3_eq)
print('Psi Partial differentiation 3 =')
sm.pprint(dpsi_x4_eq)

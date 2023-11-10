import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np

L = 50

##############################################################################################################

##############################################################################################################

A = np.zeros((16,16))
B = np.zeros((16))

# field equations 1
s = L
A[0][4], A[0][0], A[0][1], A[0][2], A[0][3], B[0] = 1, -1, -s, -(s**2)/2, -(s**3)/6, ((s-L/2)**3)/6
A[1][5], A[1][1], A[1][2], A[1][3], B[1] = 1, -1, -s, -(s**2)/2, ((s-L/2)**2)/2
A[2][6], A[2][2], A[2][3], B[2] = 1, -1, -s, (s-L/2)
A[3][7], A[3][3], B[3] = 1, -1, 1

# edge equations left
A[4][0] = 1
A[5][2] = 1

# transition equations 1-2
A[6][8], A[6][4] = 1, -1
A[7][9], A[7][5] = 1, -1
A[8][10], A[8][6] = 1, -1
A[9][4] = 1

# field equations 2
s = 2*L
A[10][12], A[10][8], A[10][9], A[10][10], A[10][11] = 1, -1, -(s-L), -((s-L)**2)/2, -((s-L)**3)/6
A[11][13], A[11][9], A[11][10], A[11][11] = 1, -1, -(s-L), -((s-L)**2)/2
A[12][14], A[12][10], A[12][11] = 1, -1, -(s-L)
A[13][15], A[13][11] = 1, -1

# edge equations right
A[14][12] = 1
A[15][14] = 1

A_inv = np.linalg.inv(A)
X = np.dot(A_inv,B)
S_tmm = np.arange(0,2*L+1,1)
W_tmm = []

for s in S_tmm:
    if s <= L/2:
        W_tmm.append(X[0] + X[1]*s + X[2]*(s**2)/2 + X[3]*(s**3)/6)
    if L/2 < s <= L:
        W_tmm.append(X[0] + X[1]*s + X[2]*(s**2)/2 + X[3]*(s**3)/6 + ((s-L/2)**3)/6)
    if L < s:
        W_tmm.append(X[8] + X[9]*(s-L) + X[10]*((s-L)**2)/2 + X[11]*((s-L)**3)/6)

plt.plot(S_tmm, W_tmm, marker="o", color="orange")
plt.grid()
plt.show()

##############################################################################################################

##############################################################################################################

A1 = np.array([
    [1, 1, 1, 1], 
    [np.exp(-3), np.exp(-2), np.exp(-1), 1],
    [9, 4, 1, 0],
    [9*np.exp(-3), 4*np.exp(-2), np.exp(-1), 0],
])

B1 = np.array([-1, -np.exp(-4), -16, -16*np.exp(-4)])

A1_inv = np.linalg.inv(A1)
a1, b1, c1, d1 = np.dot(A1_inv,B1)

A2 = np.array([
    [1, 1, 1, 1], 
    [np.exp(-2), np.exp(-1), 1, np.exp(1)],
    [4, 1, 0, 1],
    [4*np.exp(-2), np.exp(-1), 0, np.exp(1)],
])

B2 = np.array([-1, -np.exp(-3), -9, -9*np.exp(-3)])

A2_inv = np.linalg.inv(A2)
a2, b2, c2, d2 = np.dot(A2_inv,B2)

def phi_1(x):
    return np.exp(-(4*x)/(2*L)) + a1*np.exp(-(3*x)/(2*L)) + b1*np.exp(-(2*x)/(2*L)) + c1*np.exp(-x/(2*L)) + d1

def d4_phi_1(x):
    return (256*np.exp(-(4*x)/(2*L)) + 81*a1*np.exp(-(3*x)/(2*L)) + 16*b1*np.exp(-(2*x)/(2*L)) + c1*np.exp(-x/(2*L))) / (16*L**4)


def phi_2(x):
    return np.exp(-(3*x)/(2*L)) + a2*np.exp(-(2*x)/(2*L)) + b2*np.exp(-x/(2*L)) + c2 + d2*np.exp(x/(2*L))

def d4_phi_2(x):
    return (81*np.exp(-(3*x)/(2*L)) + 16*a2*np.exp(-(2*x)/(2*L)) + b2*np.exp(-x/(2*L)) + d2*np.exp(x/(2*L))) / (16*L**4)

def integrand_11(x):
    return (d4_phi_1(x)) * (phi_1(x))

def integrand_21(x):
    return (d4_phi_2(x)) * (phi_1(x))

def integrand_12(x):
    return (d4_phi_1(x)) * (phi_2(x))

def integrand_22(x):
    return (d4_phi_2(x)) * (phi_2(x))

I11 = quad(integrand_11, 0, 2*L)[0]
I21 = quad(integrand_21, 0, 2*L)[0]
I12 = quad(integrand_12, 0, 2*L)[0]
I22 = quad(integrand_22, 0, 2*L)[0]

U = np.array([
    [I11, I21, -phi_1(L)], 
    [I12, I22, -phi_2(L)],
    [phi_1(L), phi_2(L), 0]
])

V = np.array([phi_1(L/2), phi_2(L/2), 0])

U_inv = np.linalg.inv(U)
alpha1, alpha2, Z = np.dot(U_inv,V)

S_wrm = np.arange(0,2*L+1,1)
W_wrm = np.array([alpha1*phi_1(s) + alpha2*phi_2(s) for s in S_wrm])

plt.plot(S_wrm, W_wrm, marker="o", color="orange")
plt.grid()
plt.show()
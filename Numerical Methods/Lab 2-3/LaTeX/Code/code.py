import matplotlib.pyplot as plt
from scipy.integrate import quad
import numpy as np
import copy

L = 10

def Krylov(index,varphi,s):    
    if index == 1:
        return (1.0/2)*(np.cosh(varphi*s) + np.cos(varphi*s))
    if index == 2:
        return (1.0/(2*varphi))*(np.sinh(varphi*s) + np.sin(varphi*s))
    if index == 3:
        return (1.0/(2*np.power(varphi,2)))*(np.cosh(varphi*s) - np.cos(varphi*s))
    if index == 4:
        return (1.0/(2*np.power(varphi,3)))*(np.sinh(varphi*s) - np.sin(varphi*s))

step = 1e-4
w_candidates = np.arange(step, 1+step/2, step)

M2L = np.zeros(len(w_candidates))
theta0 = 1.0

w_zero_M2L = []
for i in range(len(w_candidates)):
    varphi = np.sqrt(w_candidates[i])
    X = homogeneous_equations(theta0,varphi)
    M2L[i] = X[14]

    if M2L[i] == 0:
        w_zero_M2L.append(w_candidates[i])

W_tmm = []
for w_i in w_zero_M2L:
    varphi = np.sqrt(w_i)
    X = homogeneous_equations(theta0,varphi)

    step = 0.01
    S = np.arange(0, 2*L+step/2, step)
    W_i = []

    for s in S:
        if 0 <= s <= L:
            W_i.append(X[1]*Krylov(2,varphi,s) + X[3]*Krylov(4,varphi,s))
        if L < s <= 2*L:
            W_i.append(X[9]*Krylov(2,varphi,s-L) + X[10]*Krylov(3,varphi,s-L) + X[11]*Krylov(4,varphi,s-L))
    
    W_tmm.append(copy.deepcopy(W_i))

step = 1e-4
w_candidates = np.arange(0, 1+step/2, step)
determinant = np.zeros(len(w_candidates))

w_zero_determinant = []
for i in range(len(w_candidates)):
    I = calculate_integrads_2(w_candidates[i])

    U = np.array([
        [I[1][1], I[2][1], -phi_1(L)], 
        [I[1][2], I[2][2], -phi_2(L)],
        [phi_1(L), phi_2(L), 0]
    ])

    determinant[i] = np.linalg.det(U)
    if determinant[i] == 0:
        w_zero_determinant.append(w_candidates[i])

I = calculate_integrads_2(w_zero_determinant[0])
alpha1 = 1000

U = np.array([
    [I[2][1], -phi_1(L)], 
    # [I[2][2], -phi_2(L)],
    [phi_2(L), 0]
])

V = np.array([-alpha1*I[1][1], -alpha1*phi_1(L)])

U_inv = np.linalg.inv(U)
alpha2, Z = np.dot(U_inv,V)

eta = 0.8*0.0987
w = 0.0987
varphi = np.sqrt(w)

A, B = define_inhomogeneous_matrices(varphi)
A_inv = np.linalg.inv(A)
X = np.dot(A_inv,B)

step = 0.2
S = np.arange(0, 2*L+step/2, step)

W_tmm_classic_t = []
for t in [0]:
    W_t = []
    for s in S:
        if s <= L/2:
            W_t.append((X[1]*Krylov(2,varphi,s) + X[3]*Krylov(4,varphi,s)) * np.cos(eta*t))
        if L/2 < s <= L:
            W_t.append((X[1]*Krylov(2,varphi,s) + X[3]*Krylov(4,varphi,s) - Krylov(4,varphi,s-L/2)) * np.cos(eta*t))
        if s > L:
            W_t.append((X[9]*Krylov(2,varphi,s-L) + X[10]*Krylov(3,varphi,s-L) + X[11]*Krylov(4,varphi,s-L)) * np.cos(eta*t))
    
    plt.plot(S, W_t, marker="o", color="blue")
    W_tmm_classic_t.append(copy.deepcopy(W_t))

plt.grid()
plt.show()

eta = 0.8*0.0987

scale = 1200

step = 0.2
S = np.arange(0, 2*L+step/2, step)

W_tmm_eigenvectors_t = []
for t in [40]:
    W_t = np.zeros(len(S))
    for i in range(len(S)):
        if 0 <= s <= L:
            for w_i in w_zero_M2L[0:5]:
                varphi = np.sqrt(w_i)
                W_t[i] += scale*F(S[i],1,varphi)*T(t,1,w_i,eta)
        if L < s <= 2*L:
            for w_i in w_zero_M2L[0:5]:
                varphi = np.sqrt(w_i)
                W_t[i] += scale*F(S[i]-L,1,varphi)*T(t,1,w_i,eta)
        
    plt.plot(S, W_t, marker="o", color="blue")
    W_tmm_eigenvectors_t.append(copy.deepcopy(W_t))

plt.grid()
plt.show()

eta = 0.8*0.0987
w = 0.0987
varphi = np.sqrt(w)

I = calculate_integrads(w)

U = np.array([
    [I[1][1], I[2][1], I[3][1], I[4][1], I[5][1], -phi_1(L)], 
    [I[1][2], I[2][2], I[3][2], I[4][2], I[5][2], -phi_2(L)],
    [I[1][3], I[2][3], I[3][3], I[4][3], I[5][3], -phi_3(L)],
    [I[1][4], I[2][4], I[3][4], I[4][4], I[5][4], -phi_4(L)],
    [I[1][5], I[2][5], I[3][5], I[4][5], I[5][5], -phi_5(L)],
    [phi_1(L), phi_2(L), phi_3(L), phi_4(L), phi_5(L), 0],
])

W_wrm_t = []
for t in [0]:
    V = np.array([
        phi_1(L/2)*np.cos(eta*t) / np.cos(w*t), 
        phi_2(L/2)*np.cos(eta*t) / np.cos(w*t),
        phi_3(L/2)*np.cos(eta*t) / np.cos(w*t),
        phi_4(L/2)*np.cos(eta*t) / np.cos(w*t),
        phi_5(L/2)*np.cos(eta*t) / np.cos(w*t),
        0
    ])

    U_inv = np.linalg.inv(U)
    alpha1, alpha2, alpha3, alpha4, alpha5, Z = np.dot(U_inv,V)

    step = 0.2
    S = np.arange(0, 2*L+step/2, step)
    W_t = [alpha1*phi_1(s) + alpha2*phi_2(s) + alpha3*phi_3(s) + alpha4*phi_4(s) + alpha5*phi_5(s) for s in S]
        
    plt.plot(S, W_t, marker="o", color="blue")
    W_wrm_t.append(copy.deepcopy(W_t))

plt.grid()
plt.show()
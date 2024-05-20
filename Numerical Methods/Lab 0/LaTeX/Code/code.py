import matplotlib.pyplot as plt
import numpy as np

##############################################################################################################
##############################################################################################################

def transition_matrix_method(D,K,C,triangle):
    A = np.zeros((D,D))
    B = np.zeros((D))

    if triangle == "ABC":
        L, alpha = 5, 1

        # field equations AB
        A[0][2], A[0][0], A[0][1], B[0] = 1, -1, 2*alpha*abs(C[1])*L, -abs(C[1])*C[1]*L
        A[1][3], A[1][1], B[1] = 1, -1, C[1] - C[3]

        # field equations BC
        A[2][6], A[2][4], A[2][5], B[2] = 1, -1, 2*alpha*abs(C[5])*L, -abs(C[5])*C[5]*L
        A[3][7], A[3][5], B[3] = 1, -1, C[5] - C[7]

        # field equations CA
        A[4][10], A[4][8], A[4][9], B[4] = 1, -1, 2*alpha*abs(C[9])*L, -abs(C[9])*C[9]*L
        A[5][11], A[5][9], B[5] = 1, -1, C[9] - C[11]

        # edge equations
        A[6][0], B[6] = 1, 20*20
        A[7][4], B[7] = 1, 6*6
        A[8][8], B[8] = 1, 10*10

        # transition equations A
        A[9][10], A[9][0] = 1, -1

        # transition equations B
        A[10][4], A[10][2] = 1, -1

        # transition equations C
        A[11][8], A[11][6] = 1, -1



    if triangle == "ABCD":
        # field equations
        for i,j in zip(range(0,2*K-1,2),range(0,D,4)):
            if i < 12:
                alpha, L = 1, 2.5
                A[i][2+j], A[i][0+j], A[i][1+j], B[i] = 1, -1, 2*alpha*abs(C[1+j])*L, -abs(C[1+j])*C[1+j]*L
                A[i+1][3+j], A[i+1][1+j], B[i+1] = 1, -1, C[1+j] - C[3+j]
            else:
                alpha, L = 1/3, 5*np.sqrt(3)/6
                A[i][2+j], A[i][0+j], A[i][1+j], B[i] = 1, -1, 2*alpha*abs(C[1+j])*L, -abs(C[1+j])*C[1+j]*L
                A[i+1][3+j], A[i+1][1+j], B[i+1] = 1, -1, C[1+j] - C[3+j]

        # edge equations
        A[18][0], B[18] = 1, 20*20
        A[19][8], B[19] = 1, 6*6
        A[20][16], B[20] = 1, 10*10

        # transition equations A
        A[21][22], A[21][0] = 1, -1

        # transition equations B
        A[22][8], A[22][6] = 1, -1

        # transition equations C
        A[23][16], A[23][14] = 1, -1

        # node equations D0
        A[24][32], A[24][28] = 1, -1
        A[25][32], A[25][26] = 1, -1
        A[26][27], A[26][33], A[26][29], B[26] = 1, -1, -1, C[29] + C[33] - C[27]

        # node equations D1
        A[27][20], A[27][18] = 1, -1
        A[28][20], A[28][34] = 1, -1
        A[29][35], A[29][19], A[29][21], B[29] = 1, 1, -1, C[21] - C[19] - C[35]

        # node equations D2
        A[30][12], A[30][10] = 1, -1
        A[31][12], A[31][30] = 1, -1
        A[32][11], A[32][31], A[32][13], B[32] = 1, 1, -1, C[13] - C[31] - C[11]

        # node equations D2
        A[33][24], A[33][2] = 1, -1
        A[34][24], A[34][4] = 1, -1
        A[35][3], A[35][25], A[35][5], B[35] = 1, -1, -1, C[5] + C[25] - C[3]

    A_inv = np.linalg.inv(A)
    X = np.dot(A_inv,B)
    return X

##############################################################################################################
##############################################################################################################

triangle = "ABC"
iterations = 1000
eta = 0.01

N = 2
K = 3
D = 2*N*K

C = np.zeros((D))
for i in range(1,len(C),2):
    C[i] = 5

Qs = [[C[i] for i in range(len(C)) if i % (2*N) == 1]]

for cycle in range(iterations):
    X = transition_matrix_method(D,K,C,triangle)
    Qs.append([X[i] for i in range(len(X)) if i % (2*N) == 1])
    for i in range(1,len(C),2):
        C[i] = C[i] + eta*X[i]

for j in range(len(Qs[0])):
    plt.plot([Qs[i][j] for i in range(len(Qs))], marker="o", color="blue")
    plt.show()

##############################################################################################################
##############################################################################################################

triangle = "ABCD"
iterations = 1000
eta = 0.01

N = 2
K = 9
D = 2*N*K

C = np.zeros((D))
for i in range(1,len(C),2):
    C[i] = 5

Qs = [[C[i] for i in range(len(C)) if i % (2*N) == 1]]

for cycle in range(iterations):
    X = transition_matrix_method(D,K,C,triangle)
    Qs.append([X[i] for i in range(len(X)) if i % (2*N) == 1])
    for i in range(1,len(C),2):
        C[i] = C[i] + eta*X[i]

for j in range(len(Qs[0])):
    plt.plot([Qs[i][j] for i in range(len(Qs))], marker="o", color="blue")
    plt.show()
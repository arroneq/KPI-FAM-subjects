import matplotlib.pyplot as plt
import numpy as np

N = 13
b = 1667 # insurance premium at the initial time (2022)
u1 = 576 # initial (2022) capital of the company

years = np.arange(2014,2023,1)
annual_payout_levels = np.array([20.33, 29.12, 20.1, 20.42, 25.74, 30.68, 36.29, 41.93, 43.17])/100	

plt.title("Annual payout levels")
plt.ylabel("Payout level")
plt.xlabel("Year")
plt.plot(years, annual_payout_levels, color="blue")
plt.plot(years, annual_payout_levels, "o", color="blue")
plt.show()

def generate_matrix_P(M, a, b, annual_payout_levels):
    P = np.zeros((M+1,M+1))

    P[0][0] = 1.0
    for i in range(1,len(P)):
        for j in range(len(P[0])):
            if j == M:
                P[i][j] = len(np.where(
                    i + round(b*(1-a)) - (b*annual_payout_levels).round() >= M
                    )[0]) / len(annual_payout_levels)
            elif j == 0:
                P[i][j] = len(np.where(
                    i + round(b*(1-a)) - (b*annual_payout_levels).round() <= 0
                    )[0]) / len(annual_payout_levels)
            else:
                P[i][j] = len(np.where(
                    i + round(b*(1-a)) - (b*annual_payout_levels).round() == j
                    )[0]) / len(annual_payout_levels)
    
    return P

def P_to_power_T(P, T):
    return np.linalg.matrix_power(P, T)

a = 0.7 # share of premiums spent on servicing insurance contracts
M = 2452

P = generate_matrix_P(M, a, b, annual_payout_levels)

for T in range(10, 101, 10):
    P_to_power = P_to_power_T(P, T)
    print("Default probability at time T =", T, round(P_to_power[u1][0],2))
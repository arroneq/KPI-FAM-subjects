import matplotlib.pyplot as plt
import numpy as np
import copy

N = 13
x_0 = 1667017 # insurance premium at the initial time (2022)
u_0 = 576078 # initial (2022) capital of the company

years = np.arange(2014,2023,1)
annual_payout_levels = np.array([20.33, 29.12, 20.1, 20.42, 25.74, 30.68, 36.29, 41.93, 43.17])		

plt.title("Annual payout levels")
plt.ylabel("Payout level")
plt.xlabel("Year")
plt.plot(years, annual_payout_levels, color="blue")
plt.plot(years, annual_payout_levels, "o", color="blue")
plt.show()

def x(t):
    # insurance premium
    return x_0 + t/N

def evolution_function(u_previous, v, t):
    return u_previous + (1-v)*x(t) - x(t)*np.random.choice(annual_payout_levels)/100

v = 0.7 # share of premiums spent on servicing insurance contracts
T = 31 # years of predictions
L = u_0/N
u_list = []
iterations = 1000

for i in range(iterations):
    u = np.zeros(T)
    for t in range(T):
        if t == 0:
            u[t] = u_0
        else:
            u[t] = evolution_function(u[t-1], v, t-1)
    plt.plot(np.arange(2022,2022+T,1), u)
    u_list.append(copy.deepcopy(u))

plt.show()

u_mean = np.zeros(T)
prob = np.zeros(T)

for t in range(T):    
    u_mean[t] = np.mean([u_list[i][t] for i in range(iterations)])
    prob[t] = [len(np.where(u_list[i][:t] <= L)[0]) > 0 for i in range(iterations)].count(True)/iterations

plt.plot(u_mean, prob, "o")
plt.show()
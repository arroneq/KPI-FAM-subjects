import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binom

n = 3

L = 10**6 * n**2
M = 10**3 * n
N = 10 * n
K = 10**3 * n
U = 10**4

alpha = 1/(n+1)
epsilon = 1/(n+1)
q = 1 - 1/(n+1)

pi_star = K/N

p = 1 - pow(1-M/L,365)

x = np.arange(0, N+1, 1)
cdf = binom.cdf(x, N, p*alpha)

plt.plot(x, cdf, marker="o", linestyle="-", color="b")
plt.title(f"Binomial CDF (n={N}, p={p*alpha})")
plt.xlabel("Number of Successes")
plt.ylabel("Cumulative Probability")
plt.grid()
plt.show()

p_bankruptcy = binom.cdf((U + N*pi_star)/K, N, p*alpha)
print("Probability of bankruptcy:", round(1 - p_bankruptcy,4))

pi = (K * binom.ppf(1-epsilon, N, p*alpha) - U)/N
print("Optimal cost of the policy (at least):", pi)

quantile = binom.ppf(q, N, p*alpha)
print("Quantile:", quantile)

x_pi = [pi for pi in range(0,K+1)]
y_bankruptcy = [1-binom.cdf((U + N*pi)/K, N, p*alpha) for pi in range(0,K+1)]

plt.plot(x_pi, y_bankruptcy, linestyle="-", color="b")
plt.title("Bankruptcy probability as a function of policy price")
plt.xlabel("Policy Price")
plt.ylabel("Bankruptcy Probability")
plt.grid()
plt.show()
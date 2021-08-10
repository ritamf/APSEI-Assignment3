import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

import os

from SSA import SSA
from SSAModel import SSAModel

# Selection of N (1000) people from the portuguese people with access to SearX 
N = 10**3

# Initial number of infected and recovered individuals, I0 and R0.
I0, R0 = 1, 0

# Everyone else, S0, is susceptible to infection initially.
S0 = N - I0 - R0

# Initial conditions vector
y0 = S0, I0, R0

# A grid of time points from year 2014 to 2050
YEAR_START = 2014
YEAR_END = 2050
t = np.linspace(YEAR_START, YEAR_END, 12 * (YEAR_END-YEAR_START) )

beta = 41.25/(YEAR_END-YEAR_START)
gamma = 22.625/(YEAR_END-YEAR_START)

# initial species counts and sojourn times
initital_conditions = {
    "s": [S0],
    "i": [I0],
    "r": [R0],
    "time": [YEAR_START],
}

# propensity functions
propensities = {
    0: lambda d: beta * d["s"][-1] * d["i"][-1] / N,
    1: lambda d: gamma * d["i"][-1],
}

# change in species for each propensity
stoichiometry = {
    0: {"s": -1, "i": 1, "r": 0},
    1: {"s": 0, "i": -1, "r": 1},
}

# instantiate the epidemic SSA model container
epidemic = SSAModel(
    initital_conditions,
    propensities,
    stoichiometry
)

# instantiate the SSA container with model
epidemic_generator = SSA(epidemic)

# make a nice, big figure 
fig = plt.figure(facecolor='w',figsize=(10,5))

# make a subplot for the susceptible, infected and recovered individuals
axes_s = plt.subplot(311)
axes_s.set_title(f"SIR: SearX's popularity for N={N}, beta={beta:.3f} and gamma={gamma:.3f}\n")
axes_s.set_ylabel("S: never been \n SearX users")

axes_i = plt.subplot(312)
axes_i.set_ylabel("I: are SearX \n users")

axes_r = plt.subplot(313)
axes_r.set_ylabel("R: are no longer \n SearX users")

axes_r.set_xlabel("t: Year")

# simulate and plot trajectories
trajectories = 0
for trajectory in epidemic_generator.direct():
    axes_s.plot(trajectory["time"], trajectory["s"], color="b", ls="--")
    axes_i.plot(trajectory["time"], trajectory["i"], color="r", ls="--")
    axes_r.plot(trajectory["time"], trajectory["r"], color="g", ls="--")
    trajectories += 1
    if trajectories == 10:
        break

def deriv(n_SIR, t, beta, gamma):
    dS_dt = -beta * n_SIR[0] * n_SIR[1] / N
    dI_dt = ((beta * n_SIR[0] / N) - gamma) * n_SIR[1]
    dR_dt = gamma * n_SIR[1]
    return dS_dt, dI_dt, dR_dt

solution = odeint(deriv, y0, t, args=(beta, gamma))
solution = [[row[i] for row in solution] for i in range(3)]

# plot numerical solution
axes_s.plot(t, solution[0], color="black", lw=3)
axes_i.plot(t, solution[1], color="black", lw=3)
axes_r.plot(t, solution[2], color="black", lw=3)


# save image
img_name="plts3_SIR_gillespie_algo"

os.system(f"rm {img_name}.png")

plt.savefig(f"{img_name}")
plt.close()
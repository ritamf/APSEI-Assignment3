import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

import os

from SSA import SSA
from SSAModel import SSAModel

# Selection of N (1000) people from the portuguese people with access to SearX 
N = 10**3

# Initial conditions
I0 = 1
S0 = N - I0

# Initial conditions vector
y0 = S0, I0

# A grid of time points from year 2014 to 2050
YEAR_START = 2014
YEAR_END = 2050
t = np.linspace(YEAR_START, YEAR_END, 12 * (YEAR_END-YEAR_START) )

beta = 41.25/(YEAR_END-YEAR_START)

# initial species counts and sojourn times
initital_conditions = {
    "s": [S0],
    "i": [I0],
    "time": [YEAR_START],
}

# propensity functions
propensities = {
    0: lambda d: beta * d["s"][-1] * d["i"][-1] / N,
}

# change in species for each propensity
stoichiometry = {
    0: {"s": -1, "i": 1},
    1: {"s": 1, "i": -1},
}

# instantiate the epidemic SSA model container
epidemic = SSAModel(
    initital_conditions,
    propensities,
    stoichiometry
)

def derivSI(n_SI, t, beta):
    dS_dt = -beta * n_SI[0] * n_SI[1] / N
    dI_dt = beta * n_SI[0] * n_SI[1] / N
    return dS_dt, dI_dt

# instantiate the SSA container with model
epidemic_generator = SSA(epidemic)

# make a nice, big figure 
fig = plt.figure(facecolor='w',figsize=(10,5))

# make a subplot for the susceptible, infected and recovered individuals
axes_s =fig.add_subplot(311, facecolor='#dddddd', axisbelow=True)
axes_s.grid(b=True, which='major', c='w', lw=2, ls='-')
for spine in ('right', 'top'):
    axes_s.spines[spine].set_visible(False)

axes_s.set_title(f"SIR: SearX's popularity for N={N} and beta={beta:.3f}\n")
axes_s.set_ylabel("S: never been \n SearX users")

axes_i = fig.add_subplot(312, facecolor='#dddddd', axisbelow=True)
axes_i.grid(b=True, which='major', c='w', lw=2, ls='-')
for spine in ('right', 'top'):
    axes_i.spines[spine].set_visible(False)
axes_i.set_ylabel("I: are SearX \n users")

axes_i.set_xlabel("t: Year")

# simulate and plot trajectories
trajectories = 0
for trajectory in epidemic_generator.direct():
    axes_s.plot(trajectory["time"], trajectory["s"], color="b", ls="--")
    axes_i.plot(trajectory["time"], trajectory["i"], color="r", ls="--")
    trajectories += 1
    if trajectories == 10:
        break

solution = odeint(derivSI, y0, t, args=(beta,))
solution = [[row[i] for row in solution] for i in range(2)]

# plot numerical solution
axes_s.plot(t, solution[0], color="black", lw=3)
axes_i.plot(t, solution[1], color="black", lw=3)

# save image
img_name="plts1_SI_gillespie_algo"

os.system(f"rm {img_name}.png")

plt.savefig(f"{img_name}")
plt.close()
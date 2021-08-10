import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import os

# Total population, N.
N = 10166984 * 0.76

# Initial number of infected and recovered individuals, E0, I0 and R0.
E0, I0, R0 = 0, 1, 0

# Everyone else, S0, is susceptible to infection initially.
S0 = N - E0 - I0 - R0

# A grid of time points from year 2014 to 2050
YEAR_START = 2014
YEAR_END = 2050
t = np.linspace(YEAR_START, YEAR_END, 12 * (YEAR_END-YEAR_START) )

# The SEIR model differential equations.
def deriv(y, t, N, alpha, beta, gamma, delta):
    S, E, I, R = y
    dSdt = -alpha * S * I / N
    dEdt = alpha * S * I / N - beta * E
    dIdt = beta * E - gamma * I + delta * R
    dRdt = gamma * I - delta * R
    return dSdt, dEdt, dIdt, dRdt

# Initial conditions vector
y0 = S0, E0, I0, R0

img_name="plts6_SEIRI"

# only the images generated from the last execution of this script are shown in the directory
os.system(f"rm {img_name}*.png")

num_figures = 3
for i in range(num_figures):

    # exposition rate, alpha; contact rate, beta; mean recovery rate, gamma, (in 1/days); retransmission rate, delta
    alpha, beta, gamma, delta = (47.9+i*4.75)/(YEAR_END-YEAR_START), (41.25+i*4.75)/(YEAR_END-YEAR_START), (22.625-i*4.75)/(YEAR_END-YEAR_START), (36.5+i*4.75)/(YEAR_END-YEAR_START)

    # Integrate the SEIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, alpha, beta, gamma, delta))
    S, E, I, R = ret.T

    # Plot the data on three separate curves for S(t), E(t), I(t) and R(t)
    fig = plt.figure(facecolor='w',figsize=(10,5))
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(t, S, 'b', lw=2, label='S: never been SearX users and are not familiar with it')
    ax.plot(t, E, 'y', lw=2, label='E: never been SearX users but are familiar with it')
    ax.plot(t, I, 'r', lw=2, label='I: are SearX users')
    ax.plot(t, R, 'g', lw=2, label='R: are no longer SearX users')

    ax.set_title(f"SEIRI: SearX's popularity for N={int(N)}, alpha={alpha:.3f}, beta={beta:.3f}, gamma={gamma:.3f} and delta={delta:.3f}\n")
    ax.set_xlabel('t: Year')
    ax.set_ylabel('N: Number of people with access to SearX')

    ylim = N * 1.2
    ax.set_ylim(0,ylim)
    ax.yaxis.set_tick_params(length=0)
    ax.xaxis.set_tick_params(length=0)
    ax.grid(b=True, which='major', c='w', lw=2, ls='-')
    legend = ax.legend()
    legend.get_frame().set_alpha(0.5)
    for spine in ('right', 'top'):
        ax.spines[spine].set_visible(False)

    plt.savefig(f"{img_name}{i}")
    plt.close()

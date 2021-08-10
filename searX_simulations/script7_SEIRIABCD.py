import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

import os

# Total population, N.
N = 10166984

# Initial number of infected and recovered individuals, E0, I0 and R0.
E0, I0, R0, A0, B0, C0, D0 = 0, 1, 0, 10.2863**7 * (1-0.76), 0, 0, 0

# Everyone else, S0, is susceptible to infection initially.
S0 = N - E0 - I0 - R0 - A0 - B0 - C0 - D0

# A grid of time points from year 2014 to 2050
YEAR_START = 2014
YEAR_END = 2050
t = np.linspace(YEAR_START, YEAR_END, 12 * (YEAR_END-YEAR_START) )

# The SEIR model differential equations.
def deriv(y, t, N, alpha, beta, gamma, delta, epsilon):
    S, E, I, R, A, B, C, D = y

    dSdt = -alpha * S * I / N - epsilon*S + epsilon*A
    dEdt = alpha * S * I / N - beta * E - epsilon*E + epsilon*B
    dIdt = beta * E - gamma * I + delta * R - epsilon*I + epsilon*C
    dRdt = gamma * I - delta * R - epsilon*R + + epsilon*D

    dAdt = epsilon*S - epsilon*A
    dBdt = epsilon*E - epsilon*B
    dCdt = epsilon*I - epsilon*C
    dDdt = epsilon*R - epsilon*D

    return dSdt, dEdt, dIdt, dRdt, dAdt, dBdt, dCdt, dDdt

# Initial conditions vector
y0 = S0, E0, I0, R0, A0, B0, C0, D0

img_name="plts7_SEIRIABCD"

# only the images generated from the last execution of this script are shown in the directory
os.system(f"rm {img_name}*.png")

num_figures = 3
for i in range(num_figures):

    # exposition rate, alpha; contact rate, beta; mean recovery rate, gamma, (in 1/days); retransmission rate, delta
    alpha, beta, gamma, delta, epsilon = (47.9+i*4.75)/(YEAR_END-YEAR_START), (41.25+i*4.75)/(YEAR_END-YEAR_START), (22.625-i*4.75)/(YEAR_END-YEAR_START), (36.5+i*4.75)/(YEAR_END-YEAR_START), (0+i*0.25)/(YEAR_END-YEAR_START)

    # Integrate the SEIR equations over the time grid, t.
    ret = odeint(deriv, y0, t, args=(N, alpha, beta, gamma, delta, epsilon))
    S, E, I, R, A, B, C, D = ret.T

    # Plot the data on three separate curves for S(t), E(t), I(t) and R(t)
    fig = plt.figure(facecolor='w',figsize=(10,5))
    ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
    ax.plot(t, S, 'b', lw=2, label='S: never been SearX users and are unfamiliar with it')
    ax.plot(t, E, 'y', lw=2, label='E: never been SearX users but are familiar with it')
    ax.plot(t, I, 'r', lw=2, label='I: are SearX users')
    ax.plot(t, R, 'g', lw=2, label='R: are no longer SearX users')

    ax.plot(t, S, 'b--', lw=4, label='A: S group but without access to searX')
    ax.plot(t, E, 'y--', lw=4, label='B: E group but without access to searX')
    ax.plot(t, I, 'r--', lw=4, label='C: I group but without access to searX')
    ax.plot(t, R, 'g--', lw=4, label='D: R group but without access to searX')

    ax.set_title(f"SEIRIABCD: SearX's popularity for alpha={alpha:.3f}, beta={beta:.3f}, gamma={gamma:.3f}, delta={delta:.3f} and epsilon={epsilon:.3f}\n")
    ax.set_xlabel('t: Year')
    ax.set_ylabel('N: Number of people')

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

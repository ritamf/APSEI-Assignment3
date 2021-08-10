import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

f = open("data.txt", "r") 

lines = f.readlines()

N = 7726907 * 0.76

t, S, E, I, R, u = list(), list(), list(), list(), list(), list() 

for i in range(len(lines)):
    lines[i] = lines[i].split(" ")
    
    lines[i] = [float(elem.strip()) for elem in lines[i] if elem != ""]
            
    t += [ lines[i][0] + 2014 ]
    S += [ lines[i][1] * N ]
    E += [ lines[i][2] * N ]
    I += [ lines[i][3] * N ]
    R += [ lines[i][4] * N ]
    u += [ lines[i][5] * N ]

    print(lines[i])

# Plot the data on three separate curves for S(t), E(t), I(t) and R(t)
fig = plt.figure(facecolor='w',figsize=(10,5))
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', lw=2, label='S: never been SearX users and are unfamiliar with it')
ax.plot(t, E, 'y', lw=2, label='E: never been SearX users but are familiar with it')
ax.plot(t, I, 'r', lw=2, label='I: are SearX users')
ax.plot(t, R, 'g', lw=2, label='R: are no longer SearX users')
ax.plot(t[1:-1], [N - ui for ui in u[1:-1] ], 'm--', lw=2, label='u: optimal control')

ax.set_title(f"SEIR: SearX's popularity \n")
ax.set_xlabel('t: Year')
ax.set_ylabel('N: Number of people with access to SearX')

ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('right', 'top'):
    ax.spines[spine].set_visible(False)

plt.savefig(f"plts5_SEIR_optimal_control")
plt.close()

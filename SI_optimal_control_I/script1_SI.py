import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

f = open("data.txt", "r") 

lines = f.readlines()

N = 10166984 * 0.76

t, S, I, u = [], [], [], []

for i in range(len(lines)):
    lines[i] = lines[i].split(" ")
    
    lines[i] = [float(elem.strip()) for elem in lines[i] if elem != ""]
            
    t += [ lines[i][0] + 2014 ]
    S += [ lines[i][1] * N ]
    I += [ lines[i][2] * N ]
    u += [ lines[i][3] * N ]

    print(lines[i])

# Plot the data on three separate curves for S(t), I(t) and u(t)
fig = plt.figure(facecolor='w',figsize=(10,5))
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(t, S, 'b', lw=2, label='S: people who have not used SearX')
ax.plot(t, I, 'r', lw=2, label='I: people who have used SearX')
ax.plot(t[1:-1], u[1:-1], 'm--', lw=2, label='u: optimal control')

ax.set_title(f"SI: SearX's popularity for N={int(N)} and beta={1.41:.3f}\n")
ax.set_xlabel('t: Year')
ax.set_ylabel('N: Number of people with access to SearX')

ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('right', 'top'):
    ax.spines[spine].set_visible(False)

plt.savefig(f"plts1_SI_optimal_control")
plt.close()

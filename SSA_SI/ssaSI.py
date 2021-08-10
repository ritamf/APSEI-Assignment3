import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

from SI_Simulation import SI_Simulation # Import the model.

G = nx.erdos_renyi_graph(1000,0.05) # Use Networkx to generate a random graph.
A = nx.to_numpy_array(G)
model = SI_Simulation(A, lam=0.8, gam=0.001, i0=0.0, prop="TAKE") # Setup the simulation with given parameters.
model.RunToConvergence() # Run the simulation.
model.IntegrateSolution() # Numerically integrate the mean field equations.

# Simulation.
fig = plt.figure(facecolor='w',figsize=(10,5))
ax = fig.add_subplot(111, facecolor='#dddddd', axisbelow=True)
ax.plot(model.times,model.S, 'b--', label='S: have not used SearX (sim)')
ax.plot(model.times,model.I, 'r--', label='I: have used SearX (sim)')

# Numerical Integration.
t = np.linspace(0,model.times[-1],100)

ax.plot(t, model.solution[:,0] * model.N,'b-', label='S: have not used SearX (approx)')
ax.plot(t, model.solution[:,1] * model.N,'r-', label='I: have used SearX (approx)')

ax.set_title(f"SI: SearX's popularity for N={model.N}\n")
ax.set_xlabel('t: Time')
ax.set_ylabel('N: Number of people with internet access')

ax.legend(bbox_to_anchor=(1.0,0.7))

ax.yaxis.set_tick_params(length=0)
ax.xaxis.set_tick_params(length=0)
ax.grid(b=True, which='major', c='w', lw=2, ls='-')
legend = ax.legend()
legend.get_frame().set_alpha(0.5)
for spine in ('right', 'top'):
    ax.spines[spine].set_visible(False)

plt.savefig(f"si")
plt.close()
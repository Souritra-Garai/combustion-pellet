from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

from scipy.optimize import minimize_scalar

particle_radius = np.array([50, 40, 30, 25, 20, 10])

flame_speeds = np.array([2.8, 3.5, 4.7, 5.6, 7.0, 13.7])

def fun(x) :

	return np.sum((flame_speeds - x / particle_radius)**2)

solution = minimize_scalar(fun)

x = solution.x

plt.style.use(['science', 'high-vis'])

plt.plot(particle_radius, flame_speeds, 'o', label='Simulation Results')

r = np.linspace(8, 60)
plt.plot(r, x / r, label='$y = ' + str(round(x)) + '/x $')

plt.xlabel('Particle Radius, $r_P$ (micron)')
plt.ylabel('Combustion Velocity (mm / s)')

plt.legend()

# plt.show()
plt.savefig('Particle Size Study.png', dpi=600, transparent=True)
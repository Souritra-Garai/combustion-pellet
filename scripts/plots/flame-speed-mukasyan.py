from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'high-vis'])

phi = np.array([0.24, 0.48,	0.52, 0.55, 0.59, 0.62, 0.64, 0.65, 0.68, 0.70, 0.75, 0.80, 0.90, 0.99])

flame_speeds_experiment = np.array([
	[5.2,	2.1,	2.7,	3.7,	4.2,	5.2,	5.5,	5.8,	32.4,	43.9],
	[4.9,	2.0,	2.6,	3.5,	4.0,	5.0,	5.2,	5.5,	30.8,	41.7],
	[5.5,	2.2,	2.8,	3.9,	4.4,	5.5,	5.8,	6.1,	33.9,	46.0]
])

errors = np.absolute(flame_speeds_experiment[1:3] - flame_speeds_experiment[0])

flame_speeds_du			= np.array([6.8,	8.9,	10.0,	11.2,	13.6,	16.9,	20.6,	23.3,	35.4,	44.7,	63.3,	76.8,	95.3,	106.8])
flame_speeds_alawieh	= np.array([1.0,	2.3,	2.5,	2.8,	3.4,	4.1,	5.0,	5.6,	8.2,	10.3,	14.5,	17.5,	21.7,	24.4])

plt.errorbar(
	phi[:8],
	flame_speeds_experiment[0, :8],
	errors[:, :8], 
	fmt='o',
	capsize=5,
	ecolor='black',
	label='Experimental Data')

plt.plot(phi[:8], flame_speeds_du[:8], label='Simulation (Du et al)')
plt.plot(phi[:8], flame_speeds_alawieh[:8], label='Simulation (Alawieh et al)')

plt.ylabel('Combustion Velocity (mm / s)')
plt.xlabel('Particle Volume Fractions, $\phi_P$')

plt.legend()

# plt.show()
plt.savefig('Mukasyan Zoom 1.png', dpi=600, transparent=True)
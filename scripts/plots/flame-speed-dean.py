from cProfile import label
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'high-vis'])

phi = np.array([0.55, 0.68, 0.72, 0.80, 0.83, 0.99])

flame_speeds_experiment = np.array([
	[5.0,	21.8,	50.7,	72.7,	81.8,	118.9],
	[2.9,	13.8,	32.3,	46.5,	52.1,	75.9],
	[7.0,	29.6,	69.1,	98.9,	111.3,	161.9]
])

errors = np.absolute(flame_speeds_experiment[1:3] - flame_speeds_experiment[0])

flame_speeds_du			= np.array([12.9,	39.3,	58.4,	84.4,	91.4,	116.8])
flame_speeds_alawieh	= np.array([3.0,	8.7,	12.5,	18.8,	20.3,	24.2])

flame_speeds_sundaram	= np.array([
	[9.2,	38.2,	50.5,	79.0,	88.2,	115.9],
	[6.2,	24.1,	31.9,	50.0,	55.8,	72.7]
])

plt.errorbar(
	phi,
	flame_speeds_experiment[0],
	errors.reshape((2,-1)), 
	fmt='o',
	capsize=5,
	ecolor='black',
	label='Experimental Data')

plt.plot(phi, flame_speeds_du, label='Simulation (Du et al)')
plt.plot(phi, flame_speeds_alawieh, label='Simulation (Alawieh et al)')

# plt.plot(phi, flame_speeds_sundaram[0], label='c=15')
# plt.plot(phi, flame_speeds_sundaram[1], label='c=6')

plt.ylabel('Combustion Velocity (mm / s)')
plt.xlabel('Particle Volume Fractions, $\phi_P$')

plt.legend()

# plt.show()
plt.savefig('Dean.png', dpi=600, transparent=True)
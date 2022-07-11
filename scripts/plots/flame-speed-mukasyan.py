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

flame_speeds_alawieh	= np.array([
	[1.0,	1.7,	1.7,	1.7,	1.8,	1.8,	1.9,	1.9,	2.0,	2.0,	2.2,	2.4,	3.3,	9.5],
	[1.0,	2.3,	2.5,	2.8,	3.4,	4.1,	5.0,	5.6,	8.2,	10.3,	14.5,	17.5,	21.7,	24.4],
	[1.0,	16.7,	18.1,	19.0,	19.9,	20.4,	20.8,	21.0,	21.6,	21.9,	22.4,	23.0,	24.0,	24.4]
])

plt.errorbar(
	phi[:10],
	flame_speeds_experiment[0, :],
	errors[:, :], 
	fmt='o',
	capsize=5,
	ecolor='black',
	label='Experimental Data')

plt.plot(phi, flame_speeds_du[:], label='Simulation (Du et al)')
line, = plt.plot(phi, flame_speeds_alawieh[1, :], label='Simulation (Alawieh et al)')

plt.plot(phi, flame_speeds_alawieh[0, :], '-', c=line.get_color())
plt.plot(phi, flame_speeds_alawieh[2, :], '-', c=line.get_color())
plt.fill_between(phi, flame_speeds_alawieh[0, :], flame_speeds_alawieh[2, :], alpha=0.2, color=line.get_color(), lw=0)

plt.ylabel('Combustion Velocity (mm / s)')
plt.xlabel('Particle Volume Fractions, $\phi_P$')

plt.legend()

# plt.show()
plt.savefig('Mukasyan.png', dpi=600, transparent=True)
import numpy as np
import matplotlib.pyplot as plt

MAX_ITER = 1000

def getEMT(k1, k2, v1, v2, phi12, phi22) :

	B = (v1 * phi12 * (2 * k1 - k2) + v2 * phi22 * (2 * k2 - k1)) / (v1 * phi12 + v2 * phi22)

	D = np.sqrt(B**2 + 8 * k1 * k2)

	return (D + B) / 4

def getME2(k1, k2, v1, v2, phi11, phi21) :

	f = 3 * k2 / (2 * k2 + k1)

	return (k2 * v2 * phi21 + k1 * v1 * phi11 * f) / (v2 * phi21 + v1 * phi11 * f)

def getCompFactors(v1, v2, phi11, e1) :

	phi12 = 1 - phi11

	phi21 = (e1 - v1 * phi11) / v2

	phi22 = (1 - e1 - v1 * phi12) / v2

	return phi12, phi21, phi22

def getEqn(v1, v2, phi11, k1, k2, e1) :

	phi12, phi21, phi22 = getCompFactors(v1, v2, phi11, e1)

	return getEMT(k1, k2, v1, v2, phi12, phi22) - getME2(k1, k2, v1, v2, phi11, phi21)

def getMEB(k1, k2, v1, e1) :

	v2 = 1 - v1

	phi_11_u = (2.0 * k2 + k1) / (2.0 * v1 * (k1 - k2)) - 0.0000001
	phi_11_l = 0.0

	k_e_u = getEqn(v1, v2, phi_11_u, k1, k2, e1)
	k_e_l = getEqn(v1, v2, phi_11_l, k1, k2, e1)

	__iter = 0

	while abs(phi_11_u - phi_11_l) > 0.0000001 and __iter < MAX_ITER :

		phi_11_m = 0.5 * (phi_11_u + phi_11_l)
		
		k_e_m = getEqn(v1, v2, phi_11_m, k1, k2, e1)

		if k_e_m == 0 :
		
			break

		elif k_e_l * k_e_m < 0.0 :

			phi_11_u = phi_11_m
			k_e_u = k_e_m

		else :

			phi_11_l = phi_11_m
			k_e_l = k_e_m

		__iter += 1

	phi_11_m = 0.5 * (phi_11_u + phi_11_l)

	return getME2(k1, k2, v1, v2, phi_11_m, (e1 - v1 * phi_11_m) / v2), phi_11_m

lambda_P = 163.002		# W / m-K
lambda_F = 0.0176054	# W / m-K

v_1 = np.linspace(0.1, 1, 1000)

phi_11 = np.zeros_like(v_1)
lambda_m = np.zeros_like(v_1)

for i in np.arange(v_1.shape[0]) :

	lambda_m[i], phi_11[i] = getMEB(lambda_P, lambda_F, v_1[i], 0.001)

fig1 = plt.figure()

# fig1.suptitle('Thermal Conductivity for ME2 + EMT Combined Structure')

ax1 = fig1.add_subplot()
ax2 = ax1.twinx()

plot_phi, = ax1.plot(v_1, phi_11)

ax1.grid(which='major', color='grey')
ax1.minorticks_on()
ax1.grid(which='minor', color='grey', ls='--')

ax1.set_xlabel('Particle Volume Fractions')
ax1.set_ylabel('Volume Fraction of Particles in ME2')

plot_lambda, = ax2.plot(v_1, lambda_m, color='orange')

# ax2.grid(which='major', color='grey')
# ax2.minorticks_on()
# ax2.grid(which='minor', color='grey', ls='--')

# ax2.set_xlabel('Particle Volume Fractions')
ax2.set_ylabel('Effective Thermal Conductivity (W / m-K)')

ax1.set_title('Thermal Conductivity for ME2 + EMT Combined Structure')

ax1.legend([plot_phi, plot_lambda], ['Volume Fractions in ME2', 'Thermal Conductivity'], loc='upper left')

plt.show()
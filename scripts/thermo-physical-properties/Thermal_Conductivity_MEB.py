import numpy as np
import matplotlib.pyplot as plt

MAX_ITER = 1000

def getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11) :

	f_1 = (1.0 - 2.0 * v_1 * phi_11) / 2.0
	f_2 = v_1 * phi_11 * 3.0 * k_2 / (2.0 * k_2 + k_1)

	return (f_1 * k_2 + f_2 * k_1) / (f_1 + f_2)

def getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11) :
	
	D = (2.0 * k_1 - k_2) * v_1 * (1.0 - phi_11) +  (2.0 * k_2 - k_1) * (2.0 * v_2 + 2.0 * v_1 * phi_11 - 1.0) / 2.0

	return 0.5 * (D + np.sqrt(D * D + 2.0 * k_1 * k_2))

def getThermalConductivityMEB(v_2, k_2, k_1) :
	
	v_1 = 1.0 - v_2

	phi_11_u = (2.0 * k_2 + k_1) / (2.0 * v_1 * (k_1 - k_2)) - 0.0000001
	phi_11_l = 0.0
	
	k_e_u = getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_u) - getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_u)

	k_e_l = getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_l) - getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_l)

	__iter = 0

	while abs(phi_11_u - phi_11_l) > 0.00001 and __iter < MAX_ITER :

		phi_11_m = 0.5 * (phi_11_u + phi_11_l)
		
		k_e_m = getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_m) - getThermalConductivityME2(v_1, k_1, v_2, k_2, phi_11_m)

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

	return getThermalConductivityEMT(v_1, k_1, v_2, k_2, phi_11_m), phi_11_m

lambda_P = 163.002		# W / m-K
lambda_F = 0.0176054	# W / m-K

v_1 = np.linspace(0.1, 1, 1000)

phi_11 = np.zeros_like(v_1)
lambda_m = np.zeros_like(v_1)

for i in np.arange(v_1.shape[0]) :

	lambda_m[i], phi_11[i] = getThermalConductivityMEB(1.0 - v_1[i], lambda_F, lambda_P)

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

fig1.set_size_inches(10, 4.5)
fig1.set_dpi(300)

plt.savefig('Thermal_Conductivity_MEB.png', dpi=600)

#####################################################################################################

fig2 = plt.figure()
ax = fig2.add_subplot()

phi_11 = np.linspace(0, 1, 2000)

for v_1 in [0.8] :

	ax.plot(phi_11, getThermalConductivityME2(v_1, lambda_P, 1.0 - v_1, lambda_F, phi_11), label=r'ME2 $v_1$ = ' + str(round(v_1, 2)))
	ax.plot(phi_11, getThermalConductivityEMT(v_1, lambda_P, 1.0 - v_1, lambda_F, phi_11), label=r'EMT $v_1$ = ' + str(round(v_1, 2)))

ax.grid(which='major', color='grey')
ax.minorticks_on()
ax.grid(which='minor', color='grey', ls='--')

ax.set_xlabel('Volume Fraction of Particles in ME2')
ax.set_ylabel('Effective Thermal Conductivity (W / m-K)')

ax.legend()

plt.show()
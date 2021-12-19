import numpy as np
import matplotlib.pyplot as plt

# Universal Gas Constant
R = 8.314 # J / mol. - K

class DiffusivityParameters :

	def __init__(self, pre_exponential_factor, activation_energy, name='') -> None:

		self.name = name

		self.__pre_exponential_factor = pre_exponential_factor
		self.__activation_energy = activation_energy

		pass

	def getDiffusivity(self, temperature) :

		return self.__pre_exponential_factor * np.exp(- self.__activation_energy / (R * temperature))

	def plot(self, ax, temperature) :

		ax.plot(temperature, self.getDiffusivity(temperature), label=self.name)
		pass

Du = DiffusivityParameters(9.54E-8, 26E3, 'Du')
Alawieh = DiffusivityParameters(2.56E-6, 102.191E3, 'Alawieh')

temperature_range = np.linspace(900, 3500, 1000)

fig, ax = plt.subplots()

ax.set_xlabel('Temperature (K)')

ax.set_ylabel(r'Diffusivity ($m^2 / s$)')

Du.plot(ax, temperature_range)
Alawieh.plot(ax, temperature_range)

ax.legend()

ax.grid(which='major', color='grey')
ax.minorticks_on()
ax.grid(which='minor', color='grey', ls='--')

ax.set_title('Variation of Diffusivity with Temperature for different parameters')

# fig.set_size_inches(10, 4.5)
# fig.set_dpi(300)

# plt.savefig('Thermodynamic_Properties.png', dpi=600)
plt.show()
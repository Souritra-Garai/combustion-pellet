import numpy as np
import matplotlib.pyplot as plt

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder

solution_folder = getlatestfolder()

data = np.genfromtxt(solution_folder + '/thermal_conductivity_pellet.csv', delimiter=',', skip_header=True)

fig, ax = plt.subplots()

ax.set_title('Variation of Thermal Conductivity of Pellet with Particle Volume Fractions')
ax.set_xlabel('Particle Volume Fraction')
ax.set_ylabel('Thermal Conductivity (W / m-K)')

ax.plot(data[:, 0], data[:, 1], label='ME1')
ax.plot(data[:, 0], data[:, 2], label='EMT')
ax.plot(data[:, 0], data[:, 3], label='CC + EMT')
ax.plot(data[:, 0], data[:, 4], label='ME2 + EMT')
ax.plot(data[:, 0], data[:, 5], label='CC')
ax.plot(data[:, 0], data[:, 6], label='ME2')

ax.grid(which='major', color='grey')
ax.minorticks_on()
ax.grid(which='minor', color='grey', ls='--')

ax.legend()

# fig.set_size_inches(10, 4.5)
# fig.set_dpi(300)

# plt.savefig(solution_folder + '/Thermal_Conductivity_Models.png', dpi=600)

plt.show()
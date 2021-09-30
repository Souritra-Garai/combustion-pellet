import numpy as np
import matplotlib.pyplot as plt

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder

solution_folder = getlatestfolder()

data = np.genfromtxt(solution_folder + '/Adiabatic_Combustion_Temperature.csv', delimiter=',', skip_header=True)

fig, ax = plt.subplots()

ax.set_title('Variation of Adiabatic Combustion Temperature with Initial Temperature')
ax.set_xlabel('Initial Temperature (K)')
ax.set_ylabel('Adiabatic Combustion Temperature (K)')

ax.set_ylim([1800, 3400])

ax.plot(data[:, 0], data[:, 1])

ax.grid(which='major', color='grey')
ax.minorticks_on()
ax.grid(which='minor', color='grey', ls='--')

plt.show()
import numpy as np
import matplotlib.pyplot as plt

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder

solution_folder = getlatestfolder()

data_Al     = np.genfromtxt(solution_folder + '/Thermo_Physical_Properties_Al.csv', delimiter=',', skip_header=True)
data_Ni     = np.genfromtxt(solution_folder + '/Thermo_Physical_Properties_Ni.csv', delimiter=',', skip_header=True)
data_NiAl   = np.genfromtxt(solution_folder + '/Thermo_Physical_Properties_NiAl.csv', delimiter=',', skip_header=True)

fig, axes_left = plt.subplots()

plt.figure(fig.number, figsize=(16, 9), dpi=600)

axes_left.set_xlabel('Temperature (K)')

axes_right = axes_left.twinx()

axes_left.set_ylabel('Specific Internal Energy (kJ / mol)')
axes_right.set_ylabel('Specific Heat Capacity (J / mol - K)')

for data, material in zip([data_Al, data_Ni, data_NiAl], ['Al', 'Ni', 'NiAl']) : 

    T = data[:, 0]

    c = data[:, 1]
    h = data[:, 2]

    axes_right.plot(T, c, label=material+' Heat Capacity')
    axes_left.plot(T, h / 1000, '--', label=material + ' Internal Energy')

axes_right.set_ylim([0, 100])

axes_left.legend(loc='upper left')
axes_right.legend(loc='lower right')

axes_left.grid(which='major', color='grey')
axes_left.minorticks_on()
axes_left.grid(which='minor', color='grey', ls='--')

axes_left.set_title('Variation of Thermodynamic Quantities with Temperature')

# axes_left.set_aspect('auto')
# plt.savefig('Thermodynamic_Properties.png', dpi=300)
plt.show()
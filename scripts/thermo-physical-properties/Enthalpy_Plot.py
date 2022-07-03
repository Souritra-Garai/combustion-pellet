import numpy as np
import matplotlib.pyplot as plt

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder

solution_folder = getlatestfolder()

data_Al     = np.genfromtxt(solution_folder + '/aluminium.csv', delimiter=',', skip_header=True)
data_Ni     = np.genfromtxt(solution_folder + '/nickel.csv', delimiter=',', skip_header=True)
data_NiAl   = np.genfromtxt(solution_folder + '/nickel_aluminide.csv', delimiter=',', skip_header=True)

plt.style.use('science')

fig, axes_left = plt.subplots()

axes_left.set_xlabel('Temperature (K)')

axes_right = axes_left.twinx()

axes_left.set_ylabel('Specific Enthalpy (kJ / mol)')
axes_right.set_ylabel('Specific Heat Capacity (J / mol - K)')

# for data, material in zip([data_Al, data_Ni, data_NiAl], ['Al', 'Ni', 'NiAl']) : 
for data, material in zip([data_Al], ['Al']) : 

	T = data[:, 0]

	c = data[:, 1]
	h = data[:, 2]

	axes_right.plot([],[])
	axes_right.plot(T, c, ls=':', label=' Heat Capacity')
	axes_left.plot(T, h / 1000, label='Enthalpy')
	axes_left.plot([], [], ls=':', label='Heat Capacity')

# axes_right.set_ylim([0, 1500])
# axes_right.set_yscale('log')

axes_left.legend(bbox_to_anchor=(0.4,0.53))
# axes_right.legend(loc='upper left')

# axes_left.grid(which='major', color='lightgrey')
# axes_left.minorticks_on()
# axes_left.grid(which='minor', color='lightgrey', ls='--')

# axes_left.set_title('Variation of Thermodynamic Quantities with Temperature')

# fig.set_size_inches(4 * 1.6, 4)
fig.tight_layout()
# fig.set_dpi(600)

plt.savefig(solution_folder + '/Thermodynamic_Properties.png', dpi=600, transparent=True)
# plt.show()
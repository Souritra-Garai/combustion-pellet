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

fig, ax = plt.subplots()

ax.set_xlabel('Temperature (K)')
ax.set_ylabel('Thermal Conductivity (W / m - K)')

for data, material in zip([data_Al, data_Ni, data_NiAl], ['Al', 'Ni', 'NiAl']) : 

    T = data[:, 0]

    k = data[:, 3]

    ax.plot(T, k, label=material)

ax.legend()

ax.grid(which='major', color='grey')
ax.minorticks_on()
ax.grid(which='minor', color='grey', ls='--')

ax.set_title('Variation of Thermal Conductivity with Temperature')

fig.set_size_inches(10, 4.5)
fig.set_dpi(300)

plt.savefig(solution_folder + '/Thermal_Conductivity.png', dpi=600)

# plt.show()
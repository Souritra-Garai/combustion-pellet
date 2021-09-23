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

axes_left.set_xlabel('Temperature (K)')
axes_left.set_ylabel('Specific Heat Capacity (J / kg - K)')

axes_right = axes_left.twinx()

axes_right.set_ylabel('Specific Enthalpy (J / kg)')

for data, material in zip([data_Al, data_Ni, data_NiAl], ['Al', 'Ni', 'NiAl']) : 

    T = data[:, 0]

    c = data[:, 1]
    h = data[:, 2]

    axes_left.plot(T, c, '--', label=material)
    axes_right.plot(T, h, label=material)

axes_left.set_ylim([0, 2E3])

axes_left.legend()
axes_right.legend()

plt.show()
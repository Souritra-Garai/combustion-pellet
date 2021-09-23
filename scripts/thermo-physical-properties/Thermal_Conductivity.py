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

plt.show()
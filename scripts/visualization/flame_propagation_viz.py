# Script to compare flame propagation across multiple simulation
# situated in a directory
# Place the script and run from just outside the desired directory

import os
import sys

directory = 'Time Step Convergence/'

path = os.path.join(sys.path[0], directory)

subdirectory_dict = dict()

for subdirectory in os.listdir(path) :

    subdirectory_path = os.path.join(path, subdirectory)

    if os.path.isdir(subdirectory_path) :

        subdirectory_dict[subdirectory] = subdirectory_path

import numpy as np
import matplotlib.pyplot as plt

time_list = [0.05, 0.15, 0.25]
axes_list = []

fig = plt.figure(figsize=[8,6])

for i, time in enumerate(time_list) :

    axes_list.append(fig.add_subplot(len(time_list), 1, i+1))


for subdirectory in sorted(subdirectory_dict.keys(), key=int) :

    data = np.genfromtxt(os.path.join(subdirectory_dict[subdirectory], 'temperature.csv'), delimiter=',')

    t = data[1:, 0]
    x = data[0, 1:]

    T = data[1:, 1:]

    for ax, time in zip(axes_list, time_list) :

        i = np.argmin(np.abs(t - time))

        ax.plot(x, T[i], label=subdirectory)

for ax, time in zip(axes_list, time_list) :

    ax.set_title('Time, t = ' + str(time) + ' seconds')

    ax.grid(which='major', color='grey')
    ax.minorticks_on()
    ax.grid(which='minor', color='grey', ls='--')

    ax.set_ylabel('Temperatur, T in K')
    ax.legend()

ax.set_xlabel('x in m')

fig.suptitle('Flame propagation for different Time Steps (in $\mu s$)')
plt.show()


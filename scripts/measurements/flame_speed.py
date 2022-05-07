import numpy as np
import matplotlib.pyplot as plt

import os
import sys

import argparse

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder, getpath

folder = getlatestfolder()

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--folderpath", help="Path to folder containing temperature.csv")

args = parser.parse_args()

if args.folderpath :

	folder = args.folderpath

print('Processing solution at the directory ' + folder)

data = np.genfromtxt(os.path.join(folder, 'temperature.csv'), delimiter=',')

t = data[1:, 0]
x = data[0, 1:]

T = data[1:, 1:]

flame_front_temperature = 1000	# K

flame_locations = x[np.argmin(np.abs(T - flame_front_temperature), axis=1)]	# m

valid_indices = np.where(flame_locations[1:] > flame_locations[:-1])
# print(flame_locations[valid_indices])
# print(flame_locations[valid_indices].shape)

flame_locations = flame_locations[valid_indices]
t = t[valid_indices]

valid_indices = np.where(t < 0.62)

flame_locations = flame_locations[valid_indices]
t = t[valid_indices]

line, cov = np.polyfit(t, flame_locations, 1, cov=True)
# print(line)

print('Flame Speed : ', line[0]*1000, ' mm/s')

fig = plt.figure()
ax = fig.add_subplot()

ax.plot(t, flame_locations, label='Simulation Results')
ax.plot(t, np.polyval(line, t), label='Line fit', lw=1)

ax.set_title('Flame Speed : ' +  str(line[0]*100) + ' cm/s')

ax.set_xlabel('Time (seconds)')
ax.set_ylabel('Combustion Front Location')

ax.grid(which='major', color='grey')
ax.minorticks_on()
ax.grid(which='minor', color='grey', ls='--')

ax.legend()

# fig.set_size_inches(10, 4.5)
# fig.set_dpi(600)

# plt.savefig(folder + '/Flame_Position_2000.png', dpi=600)
plt.show()
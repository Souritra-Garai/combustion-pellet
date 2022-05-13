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
x = data[0, 1:] * 1E3

T = data[1:, 1:]

N = t.shape[0]

# Number of line plots
n = 10

fig = plt.figure(figsize=[4*1.6, 4], constrained_layout=True)

line_plot_axes = fig.add_subplot(1, 1, 1)

for i in range(n) :

	line_plot_axes.plot(x, T[i * N // n], label=u'$t$ = {:.0f} ms'.format(t[i * N // n] * 1E3))
	# line_plot_axes.plot([], [], label=u't = {:.0f} ms'.format(t[i * N // n] * 1E3))

line_plot_axes.plot(x, T[-1], label=u"$t$ = {:.0f} ms".format(t[-1] * 1E3))
# line_plot_axes.plot([], [], label=u"t = {:.0f} ms".format(t[-1] * 1E3))

line_plot_axes.set_xlabel('$x$ (mm)')
line_plot_axes.set_ylabel('Temperature (K)')

# line_plot_axes.set_title('Temperature Evolution in Pellet')

line_plot_axes.grid(which='major', color='lightgrey')
line_plot_axes.minorticks_on()
line_plot_axes.grid(which='minor', color='lightgrey', ls='--')

# line_plot_axes.set_axis_off()
line_plot_axes.legend()

plt.show()
# plt.savefig(folder + '/flame_propagation.png', dpi=600, transparent=True)
# plt.savefig(folder + '/flame_propagation_legend.png', dpi=600, transparent=True)
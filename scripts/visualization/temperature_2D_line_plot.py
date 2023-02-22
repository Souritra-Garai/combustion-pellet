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
parser.add_argument("-s", "--savefigure", help="Save temperature plot", action='store_true')

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

fig = plt.figure(constrained_layout=True)

line_plot_axes = fig.add_subplot(1, 1, 1)

def genLabel(time) :

	if t[-1] > 1 :

		return u't = {:.2f} s'.format(time)
	
	else :

		return u't = {:.0f} ms'.format(time * 1E3)

for i in range(n) :

	line_plot_axes.plot(x, T[i * N // n], label=genLabel(t[i * N // n]))

line_plot_axes.plot(x, T[-1], label=genLabel(t[-1]))

line_plot_axes.set_xlabel('x (mm)')
line_plot_axes.set_ylabel('T (K)')

line_plot_axes.set_title('Temperature Evolution in Pellet')

line_plot_axes.grid(which='major', color='grey')
line_plot_axes.minorticks_on()
line_plot_axes.grid(which='minor', color='grey', ls='--')

line_plot_axes.legend(bbox_to_anchor=(1.04, 0.5), loc="center left", borderaxespad=0)

plt.show()

if args.savefigure :

	fig.savefig(folder + '/temperature_plot.jpeg', dpi=600)
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

dataA = np.genfromtxt(os.path.join(folder, 'concentration_A.csv'), delimiter=',')
dataB = np.genfromtxt(os.path.join(folder, 'concentration_B.csv'), delimiter=',')

index1 = 40
index2 = 70

t = dataA[index1:index2, 0]
x = dataA[0, 1:]

C_A = dataA[index1:index2, 1:]
C_B = dataB[index1:index2, 1:]

N = t.shape[0]

# Number of line plots
n = 5

fig = plt.figure()

line_plot_axes = fig.add_subplot(2, 1, 1)

for i in range(0, n) :

	line_plot_axes.plot(x, C_A[i * N // n], label=u'Time t = {:.3f} seconds'.format(t[i * N // n]))

line_plot_axes.plot(x, C_A[-1], label=u"Time t = {:.3f} seconds".format(t[-1]))

line_plot_axes.set_xlabel('x (m)')
line_plot_axes.set_ylabel('Concentration of Aluminium ($mol. / m^3$)')

line_plot_axes.set_title('Atomic Concentration Evolution in a Particle')

line_plot_axes.grid(which='major', color='grey')
line_plot_axes.minorticks_on()
line_plot_axes.grid(which='minor', color='grey', ls='--')

line_plot_axes.legend()

line_plot_axes = fig.add_subplot(2, 1, 2)

for i in range(0, n) :

	line_plot_axes.plot(x, C_B[i * N // n], label=u'Time t = {:.3f} seconds'.format(t[i * N // n]))

# line_plot_axes.plot(x, C_A[-1], label=u"Concentration of Al Time t = {:.3f} seconds".format(t[-1]))
line_plot_axes.plot(x, C_A[-1], label=u"Time t = {:.3f} seconds".format(t[-1]))

line_plot_axes.set_xlabel('x (m)')
line_plot_axes.set_ylabel('Concentration of Nickel ($mol. / m^3$)')

# line_plot_axes.set_title('Atomic Concentration Evolution in a Particle')

line_plot_axes.grid(which='major', color='grey')
line_plot_axes.minorticks_on()
line_plot_axes.grid(which='minor', color='grey', ls='--')

line_plot_axes.legend()

plt.show()
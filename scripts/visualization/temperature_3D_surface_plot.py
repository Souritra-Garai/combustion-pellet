import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder, getpath

folder = getlatestfolder()

print('Processing solution at the directory ' + folder)

data = np.genfromtxt(os.path.join(folder, 'temperature.csv'), delimiter=',')

t = data[1:, 0]
x = data[0, 1:]

T = data[1:, 1:]

N = t.shape[0]

t_arr, x_arr = np.meshgrid(t, x, indexing='ij')

fig = plt.figure()

surface_plot_axes = fig.add_subplot(1, 1, 1, projection='3d')

temperature_surface = surface_plot_axes.plot_surface(t_arr, x_arr, T, cmap='magma')

surface_plot_axes.set_xlabel('t (s)')
surface_plot_axes.set_ylabel('x (m)')
surface_plot_axes.set_zlabel('Temperature (K)')

surface_plot_axes.set_title('Temperature Evolution in the Pellet')

plt.show()
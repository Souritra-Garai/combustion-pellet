import numpy as np
import matplotlib.pyplot as plt

from matplotlib import cm
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation, FFMpegWriter
from numpy.core.defchararray import array

from particle_viz import ParticleDiffusion3D

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.split(sys.path[0])[0]))

from scripts.utilities.solution_folder import getlatestfolder, getpath

solution_folder = getlatestfolder()

data = np.genfromtxt(solution_folder + '/temperature.csv', delimiter=',')[:-1]

x_array = data[0, 1:]
t_array = data[1:, 0]

T_matrix = data[1:, 1:]

normalizer = Normalize(np.min(T_matrix), np.max(T_matrix))

temperature_colormap = cm.get_cmap('plasma', 1000)

fig = plt.figure(figsize=[8,6])

ax = fig.add_subplot(1, 1, 1, projection='3d')

particles = []

diameter = x_array[-1] * 1.67
n = 50

for x, T in zip(x_array, T_matrix[0]) :

	r = np.sqrt(np.random.random(n)) * diameter / 4
	theta = np.random.random(n) * 2 * np.pi

	y = r * np.cos(theta)
	z = r * np.sin(theta)

	particles.append(ax.plot(x * np.ones_like(y), y, z, marker='o', linestyle='', markersize=4, c=temperature_colormap(normalizer(T)))[0])

particle1 = ParticleDiffusion3D(solution_folder + '/1', position=[x_array[1], diameter/2, diameter/2], scale=20)
particle1.setUpPlot(ax)

particle2 = ParticleDiffusion3D(solution_folder + '/333', position=[x_array[333], diameter/2, diameter/2], scale=20)
particle2.setUpPlot(ax)

particle3 = ParticleDiffusion3D(solution_folder + '/666', position=[x_array[666], diameter/2, diameter/2], scale=20)
particle3.setUpPlot(ax)

particle4 = ParticleDiffusion3D(solution_folder + '/999', position=[x_array[999], diameter/2, diameter/2], scale=20)
particle4.setUpPlot(ax)

y = np.linspace(0, diameter / 2)
z = (diameter / 4) + 2 * y
z[z > diameter/2] = diameter/2

for x in x_array[[1, 333, 666, 999]] :

	ax.plot(np.ones_like(y)*x, y, z, color='black', lw=0.5)

theta = np.linspace(0, 2 * np.pi, 1000)
ax.plot(x_array[-1] * np.ones_like(theta), diameter * np.cos(theta) / 4, diameter * np.sin(theta) / 4)

def update(i) :

	ax.set_title(u"Time t = {:.3f} seconds".format(t_array[i]), y=1.0, pad=-14)

	for particle, T in zip(particles, T_matrix[i]) :

		particle.set_color(temperature_colormap(normalizer(T)))

	arr = [particles]

	arr.extend(particle1.update(i))
	arr.extend(particle2.update(i))
	arr.extend(particle3.update(i))
	arr.extend(particle4.update(i))

	return arr

particle1.update(0)
particle2.update(0)
particle3.update(0)
particle4.update(0)

ax.set_axis_off()

e = 1.1

lim = np.array([-0.05 * diameter, 0.6 * diameter])

ax.set_xlim(lim)
ax.set_ylim(lim)
ax.set_zlim(lim)

my_animation = FuncAnimation(fig, update, frames = np.arange(0, t_array.shape[0], 1))

ax.view_init(azim=-60, elev=15)

fig.colorbar(cm.ScalarMappable(norm=normalizer, cmap=temperature_colormap), ax=ax, label='Temperature (K)')

# fig.set_size_inches(8 * 1.5, 8 * 1.5)
# fig.set_dpi(300)

plt.show()

# writervideo = FFMpegWriter(fps=30)
# my_animation.save('Pellet_Flame_Propagation.mp4', writer=writervideo)
# plt.close()
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

solution_folder = "solutions/Video Animation" # getlatestfolder()

data = np.genfromtxt(solution_folder + '/temperature.csv', delimiter=',')[:-1]

x_array = data[0, 1:]
t_array = data[1:, 0]

T_matrix = data[1:, 1:]

plt.style.use('science')

normalizer = Normalize(np.min(T_matrix), 2000)

temperature_colormap = cm.get_cmap('plasma', 1000)

fig = plt.figure(figsize=[8,6])

ax = fig.add_subplot(1, 1, 1, projection='3d')

particles = []

diameter = x_array[-1] * 1.67
n = 50

theta = np.linspace(0, 2 * np.pi, 1000)
ax.plot(x_array[0] * np.ones_like(theta), diameter * np.cos(theta) / 4, diameter * np.sin(theta) / 4, color='skyblue')

for x, T in zip(x_array, T_matrix[0]) :

	r = np.sqrt(np.random.random(n)) * diameter / 4
	theta = np.random.random(n) * 2 * np.pi

	y = r * np.cos(theta)
	z = r * np.sin(theta)

	particles.append(ax.plot(x * np.ones_like(y), y, z, marker='o', linestyle='', markersize=4, c=temperature_colormap(normalizer(T)))[0])

theta = np.linspace(0, 2 * np.pi, 1000)
ax.plot(x_array[-1] * np.ones_like(theta), diameter * np.cos(theta) / 4, diameter * np.sin(theta) / 4, color='skyblue')

def update(i) :

	# ax.set_title(u"Time t = {:.3f} seconds".format(t_array[i]), y=1.0, pad=-14)

	for particle, T in zip(particles, T_matrix[i]) :

		particle.set_color(temperature_colormap(normalizer(T)))

	return [particles]

update(80)

ax.set_axis_off()

# my_animation = FuncAnimation(fig, update, frames = np.arange(0, t_array.shape[0], 1))

ax.view_init(azim=-60, elev=15)

cb = fig.colorbar(cm.ScalarMappable(norm=normalizer, cmap=temperature_colormap), ax=ax)
cb.set_label('Temperature (K)', labelpad=30, size=20)
cb.ax.tick_params(labelsize=20)

fig.tight_layout()
fig.set_size_inches(8 * 1.5, 6 * 1.5)
fig.set_dpi(600)

# plt.show()
plt.savefig(solution_folder + '/Flame Propagation_new.png', transparent=True)

# writervideo = FFMpegWriter(fps=30)
# my_animation.save('Pellet_Flame_Propagation.mp4', writer=writervideo)
# plt.close()
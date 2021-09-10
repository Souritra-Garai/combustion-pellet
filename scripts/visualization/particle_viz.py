# To visualise the diffusion process happenning in a core shell particle

import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

class SphericalDiffusion3d :

    def __init__(self, radial_coordinates_array, num_points = 1000, position = [0, 0, 0]) :

        self.__num_points = num_points
        
        self.__center = position
        
        self.__radial_coordinates = radial_coordinates_array

        self.__theta    = np.random.random(self.__num_points) * np.math.pi
        self.__phi      = np.random.random(self.__num_points) * np.math.pi * 2

        self.__mpl_objs = None

    def __getradii(self, probability_array) :

        return np.sort(np.random.choice(
            self.__radial_coordinates,
            self.__num_points,
            True,
            probability_array
        ))

    def __transform2cartesian(self, probability_array) :

        r = self.__getradii(probability_array)

        x = self.__center[0] + r * np.sin(self.__theta) * np.cos(self.__phi)
        y = self.__center[1] + r * np.sin(self.__theta) * np.sin(self.__phi)
        z = self.__center[2] + r * np.cos(self.__theta)

        return (x, y, z)

    def set_up_plot(self, axes, colour) :

        self.__mpl_objs, = axes.plot([], [], [], "o", markersize=0.5, c = colour)

        return self.__mpl_objs,

    def update(self, probability_array) :
        
        x, y, z = self.__transform2cartesian(probability_array)

        self.__mpl_objs.set_data(x, y)
        self.__mpl_objs.set_3d_properties(z)

        return self.__mpl_objs,

    def getplotlims(self, extension_factor = 1.2) :

        return [- extension_factor * self.__radial_coordinates[-1], extension_factor * self.__radial_coordinates[-1]]

class ParticleDiffusion3D :

    def __init__(self, solution_path, num_points = 1000, position = [0, 0, 0], scale = 1) :

        conc_A_data = np.genfromtxt(solution_path + '/concentration_A.csv', delimiter=',')
        conc_B_data = np.genfromtxt(solution_path + '/concentration_B.csv', delimiter=',')

        self.__time_array = conc_A_data[1:, 0]
        self.__radial_coordinates_array = scale * conc_A_data[0, 1:]

        self.__spherical_diffusion_A = SphericalDiffusion3d(self.__radial_coordinates_array, num_points, position)
        self.__spherical_diffusion_B = SphericalDiffusion3d(self.__radial_coordinates_array, num_points, position)

        self.__concentration_matrix = {
            'A' : conc_A_data[1:, 1:],
            'B' : conc_B_data[1:, 1:]
        }

        self.__concentration_matrix['A'][conc_A_data[1:, 1:] < 0] = 0
        self.__concentration_matrix['B'][conc_B_data[1:, 1:] < 0] = 0

    def __get_shell_probability(self, material, time_index) :

        shell_prob = self.__concentration_matrix[material][time_index] * self.__radial_coordinates_array ** 2

        return shell_prob / np.sum(shell_prob)

    def set_up_plot(self, axes) :

        mpl_objs = []

        mpl_objs.extend(self.__spherical_diffusion_A.set_up_plot(axes, 'red'))
        mpl_objs.extend(self.__spherical_diffusion_B.set_up_plot(axes, 'blue'))

        return mpl_objs

    def update(self, time_index) :

        mpl_objs = []

        mpl_objs.extend(self.__spherical_diffusion_A.update(self.__get_shell_probability('A', time_index)))
        mpl_objs.extend(self.__spherical_diffusion_B.update(self.__get_shell_probability('B', time_index)))

        return mpl_objs

    def get_time(self, time_index) :

        return self.__time_array[time_index]

    def get_plot_lims(self, extension_factor = 1.2) :

        return self.__spherical_diffusion_A.getplotlims(extension_factor)

if __name__ == '__main__' :

    from matplotlib.animation import FuncAnimation, FFMpegWriter

    from utilities.solution_folder import getlatestfolder, getpath

    particle = ParticleDiffusion3D(getlatestfolder(), 10000, [0,0,0], 1E6)

    fig = plt.figure(figsize=[10,10])
    fig.suptitle('Diffusion in Ni-Coated Al Particle', fontweight='bold')

    ax = fig.add_subplot(1, 1, 1, projection='3d')

    e = 1.1

    ax.set_xlim(particle.get_plot_lims(e))
    ax.set_ylim(particle.get_plot_lims(e))
    ax.set_zlim(particle.get_plot_lims(e))

    ax.set_xlabel(r'x $\left(\mu m\right)$')
    ax.set_ylabel(r'y $\left(\mu m\right)$')
    ax.set_zlabel(r'z $\left(\mu m\right)$')

    particle.set_up_plot(ax)

    ax.plot([],[],[], 'o', markersize=5, c='blue', label='Nickel')
    ax.plot([],[],[], 'o', markersize=5, c='red', label='Aluminium')
    ax.legend()

    # title = ax.text(0.5,0.85, "", bbox={'facecolor':'w', 'alpha':0.5, 'pad':5},
    #             transform=ax.transAxes, ha="center", s=0)

    def update(i) :

        ax.set_title(u"Time t = {:.3f} seconds".format(particle.get_time(i)))

        return particle.update(i)

    my_anim = FuncAnimation(fig, update, 5000, blit=False, interval = 1)

    # plt.show()

    writervideo = FFMpegWriter(fps=60)
    my_anim.save('Core_Shell_Particle_Diffusion.mp4', writer=writervideo)
    plt.close()
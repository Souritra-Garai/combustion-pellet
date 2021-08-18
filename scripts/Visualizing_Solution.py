import os
import sys

import numpy as np

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

all_folders = [os.path.join(os.path.dirname(sys.path[0]), 'solutions/' + d) for d in os.listdir(os.path.join(os.path.dirname(sys.path[0]), 'solutions'))]

# print(all_folder)

folder = max(all_folders, key=os.path.getmtime)

print('Visualizing solution from\n' + folder)

T = np.genfromtxt(os.path.join(folder, 'temperature.csv'), delimiter=',')[:, :-1]
# X = np.genfromtxt(os.path.join(folder, 'Conversion.csv'), delimiter=',')[:, :-1]

N, M = T.shape

Dt = 1E-5
Dx = 6.35E-3 / (M-1)
L = 6.35E-3

# with open(os.path.join(folder, 'Combustion_Config.txt')) as file :

#     for line in file :

#         mylist = list(line.split('\t'))

#         if mylist[0] == 'Delta t :' :

#             Dt = float(mylist[1])

#             # print(mylist, Dt)

#         if mylist[0] == 'Delta x :' :

#             Dx = float(mylist[1])

#             # print(mylist, Dt)

#         if mylist[0] == 'Length :' :

#             L = float(mylist[1])

#             # print(mylist, Dt)        

t = np.linspace(0, (N-1)*Dt, N)
x = np.linspace(0, L, M)

t_arr, x_arr = np.meshgrid(t, x, indexing='ij')

fig = plt.figure(figsize=plt.figaspect(1/(1.6*2)))
fig.suptitle('Combustion of Packed Pellets of Ni-Coated Al Particles')

ax1 = fig.add_subplot(1, 1, 1, projection='3d')

# Plot the surface.
surf = ax1.plot_surface(t_arr, x_arr, T, cmap='magma')

ax1.set_xlabel('t (s)')
ax1.set_ylabel('x (m)')
ax1.set_zlabel('Temperature (K)')

ax1.set_title('Temperature Evolution in the Pellet')

ax1.set_zlim([250, 2500])
# ax.set_zscale('log')

# ax2 = fig.add_subplot(1, 2, 2)

# combustion_front = []

# for eta in X :

#     for i in range(M) :

#         if np.isclose(eta[i], 0.000001) :

#             break

#     combustion_front.append(i*Dx)

# ax2.plot(t, combustion_front, color='black')

# ax2.grid(which='major', color='green')
# ax2.minorticks_on()
# ax2.grid(which='minor', color='green', ls='--')

# ax2.set_title('Propagation of Combustion Front in Pellet')

# ax2.set_xlabel('t (s)')
# ax2.set_ylabel('x (m)')

plt.show()
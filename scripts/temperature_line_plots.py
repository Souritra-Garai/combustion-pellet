import os
import sys

import numpy as np

import matplotlib.pyplot as plt

all_folders = [os.path.join(os.path.dirname(sys.path[0]), 'solutions/' + d) for d in os.listdir(os.path.join(os.path.dirname(sys.path[0]), 'solutions'))]

# print(all_folder)

folder = max(all_folders, key=os.path.getmtime)

print('Visualizing solution from\n' + folder)

T = np.genfromtxt(os.path.join(folder, 'temperature.csv'), delimiter=',')[1:, 1:]
# X = np.genfromtxt(os.path.join(folder, 'Conversion.csv'), delimiter=',')[:, :-1]

N, M = T.shape

Dt = 1E-4
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

ax1 = fig.add_subplot(1, 1, 1)

# Plot the surface.
n = 5

for i in range(n):

    ax1.plot(x, T[i*N//n], label=u"Time t = {:.3f} seconds".format(t[i*N//n]))

# ax1.plot(x, T[0], label=u"Time t = {:.3f} seconds".format(0))
# ax1.plot(x, T[N//4], label=u"Time t = {:.3f} seconds".format(t[N//4]))
# ax1.plot(x, T[2*N//4], label=u"Time t = {:.3f} seconds".format(t[2*N//4]))
# ax1.plot(x, T[3*N//4], label=u"Time t = {:.3f} seconds".format(t[3*N//4]))
ax1.plot(x, T[N-1], label=u"Time t = {:.3f} seconds".format(t[N-1]))

ax1.set_xlabel('x (m)')
ax1.set_ylabel('T (K)')
# ax1.set_zlabel('Temperature (K)')

ax1.set_title('Temperature Evolution in the Pellet')

ax1.grid(which='major', color='green')
ax1.minorticks_on()
ax1.grid(which='minor', color='green', ls='--')

ax1.legend()
# ax1.set_zlim([250, 2500])
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
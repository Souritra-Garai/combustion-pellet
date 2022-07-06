from cmath import nan
from time import time
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'high-vis'])

num_grid_points = np.array([50, 100, 250, 500, 1000])

flame_speed = np.array([10.9972, 10.9795, 10.9740, 10.9696, 10.9692])

plt.plot(num_grid_points, flame_speed, '-o')

plt.xscale('log')

plt.xlabel('Number of Grid Points, $M$')
plt.ylabel('Combustion Velocity (mm /s)')

# plt.legend()
# plt.show()

plt.savefig('Grid Size Convergence.png', dpi=600, transparent=True)

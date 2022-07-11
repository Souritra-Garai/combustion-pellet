from cmath import nan
from time import time
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'high-vis'])

num_grid_points = np.array([26, 31, 41, 51, 101, 251, 501, 1001, 1501, 2001])

flame_speed = np.array([8.4777, 8.5523, 8.6571, 8.6802, 8.6825, 8.6878, 8.6850, 8.6856, 8.6857, 8.6857])

scatter, = plt.plot(num_grid_points, flame_speed)

# x = np.linspace(50, 2000)
# poly = np.poly1d(np.polyfit(np.log(num_grid_points), flame_speed, 2))

# plt.plot(x, poly(np.log(x)), c=scatter.get_color())

plt.xlim(left=10)
plt.ylim([8.45, 8.8])
plt.xscale('log')

plt.xlabel('Number of Grid Points, $M$')
plt.ylabel('Combustion Velocity (mm /s)')

# plt.legend()
# plt.show()

plt.savefig('Grid Size Convergence.png', dpi=600, transparent=True)

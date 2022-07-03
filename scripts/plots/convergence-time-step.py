from cmath import nan
from time import time
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'high-vis'])

time_step = np.array([10, 5, 1, 0.5, 0.1])

flame_speed_explicit	= np.array([np.nan,	np.nan,	8.6023,	8.6023,	8.6024])
flame_speed_implicit	= np.array([8.5973,	8.6000,	8.6027,	8.6030,	8.6037])
flame_speed_crank		= np.array([8.5991,	8.6009,	8.6026,	8.6029,	8.6032])

plt.plot(time_step, flame_speed_implicit, label='Implicit Method')
plt.plot(time_step, flame_speed_crank, label='Crank Nicolson Method')
plt.plot(time_step[2:], flame_speed_explicit[2:], label='Explicit Method')

plt.xscale('log')

plt.xlabel('Time Step, $\Delta t$ ($\mu s$)')
plt.ylabel('Combustion Velocity (mm /s)')

plt.legend()
# plt.show()

plt.savefig('Time Step Convergence.png', dpi=600, transparent=True)

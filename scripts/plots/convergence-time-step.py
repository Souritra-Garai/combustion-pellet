from cmath import nan
from time import time
import numpy as np
import matplotlib.pyplot as plt

plt.style.use(['science', 'high-vis'])

time_step   = np.array([1000, 500, 100, 50, 10, 5, 1, 0.5, 0.1])
time_step_e = np.array([3, 2, 1, 0.5, 0.1])

flame_speed_explicit	= np.array([8.6868,	8.6884,	8.6883, 8.6883, 8.6889])
flame_speed_implicit	= np.array([8.4809, 8.5586, 8.6513, 8.6723 ,8.6793, 8.6820, 8.6823, 8.6836, 8.6839])
flame_speed_crank		= np.array([8.4997, 8.5869, 8.6680, 8.6825, 8.6831, 8.6846, 8.6856, 8.6861, 8.6854])

scatter_implicit, = plt.plot(time_step,   flame_speed_implicit, label='Implicit scheme')
scatter_crank,    = plt.plot(time_step,   flame_speed_crank,    label='Crank Nicolson scheme')
scatter_explicit, = plt.plot(time_step_e, flame_speed_explicit, label='Explicit scheme')

# x = np.linspace(0.1, 10)
# poly_implicit = np.poly1d(np.polyfit(time_step, flame_speed_implicit, 2))
# poly_crank    = np.poly1d(np.polyfit(time_step, flame_speed_crank, 2))

# x_e = np.linspace(0.1, 3)
# indices = [1, 3, 4]
# poly_explicit = np.poly1d(np.polyfit(time_step_e[indices], flame_speed_explicit[indices], 1))

# plt.plot(x, poly_implicit(x), c=scatter_implicit.get_color())
# plt.plot(x, poly_crank(x),    c=scatter_crank.get_color())
# plt.plot(x_e, poly_explicit(x), c=scatter_explicit.get_color())

# plt.axvspan(0.01, 3, facecolor=scatter_explicit.get_color(), alpha=0.2)

plt.annotate("Explicit\nscheme\nunstable\nfor $\Delta t > 3\mu s$",
            xy=(3, 8.6), xycoords='data',
            xytext=(10, 8.6), textcoords='data',
            arrowprops=dict(arrowstyle="<-",
                            connectionstyle="arc3"),
			verticalalignment='center_baseline'
            )

plt.plot([3, 3], [8.69, 8.57], c='black', alpha=0.8)

plt.xscale('log')

plt.xlim(left=0.09)

plt.xlabel('Time Step, $\Delta t$ ($\mu s$)')
plt.ylabel('Combustion Velocity (mm /s)')

plt.legend()
# plt.show()

plt.savefig('Time Step Convergence.png', dpi=600, transparent=True)

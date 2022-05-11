import numpy as np
import matplotlib.pyplot as plt

from particle_viz import ParticleDiffusion2D


Al_color = '#ADEFD1FF'
Ni_color = '#00203FFF'

particle = ParticleDiffusion2D('solutions/Video Animation/333', 10000, [0,0,0], 1E6)

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)

particle.setUpPlot(ax, Al_color, Ni_color)

# ax.plot([],[], 'o', markersize=5, color=Ni_color, label='Nickel')
# ax.plot([],[], 'o', markersize=5, color=Al_color, label='Aluminium')
# ax.legend()

ax.set_axis_off()

ax.set_xlim(particle.getPlotLims())
ax.set_ylim(particle.getPlotLims())

particle.update(75)

fig.set_size_inches(4, 4)
fig.set_dpi(600)

# plt.show()

plt.savefig('solutions/Video Animation/333/particle_3.png', transparent=True)

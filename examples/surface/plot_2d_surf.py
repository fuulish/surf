#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ase.units import Bohr

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

dat = np.loadtxt('out_0_2dsurfdown.dat')*Bohr

x = dat[0,1:].copy()
y = dat[1:,0].copy()

X, Y = np.meshgrid(x,y)

#ax.plot_surface(X, Y, dat[1:,1:])
ax.plot_wireframe(X, Y, dat[1:,1:])

dat = np.loadtxt('out_0_2dsurfup.dat')*Bohr

x = dat[0,1:].copy()
y = dat[1:,0].copy()

X, Y = np.meshgrid(x,y)

#ax.plot_surface(X, Y, dat[1:,1:])
ax.plot_wireframe(X, Y, dat[1:,1:])

plt.show()

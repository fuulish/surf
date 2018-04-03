#!/usr/bin/env python

import numpy as np
from cube import Cube
import matplotlib.pyplot as plt
from ase.units import Bohr

nstrt = 0
nstop = 500
strde = 10

nsnap = (nstop - nstrt ) / strde

zprof = None
for i in range(nstrt, nstop, strde):
    c = Cube('out_%i_instant-surface.cube' %i)

    z = np.average(np.average(c.data, axis=0), axis=0)

    if zprof is None:
        zprof = z.copy()
    else:
        zprof += z

zprof *= 1. / (Bohr**3 * nsnap)

np.savetxt('zprof_avdens.dat', zprof)

plt.plot(zprof)
plt.show()

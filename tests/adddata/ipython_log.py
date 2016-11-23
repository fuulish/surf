# IPython log file

import os
import numpy as np

os.system("awk 'NF == 5' dip.lammps > dip.out")

dat = np.loadtxt('dip.out')
np.savetxt('dip_clean.dat', dat, fmt="%14.8f")

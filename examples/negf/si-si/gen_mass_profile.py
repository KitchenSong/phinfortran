import numpy as np
import matplotlib.pyplot as plt
import random

m1 = 28.0855
m2 = 72.64 
nat = 4

nblochx = 2 
nblochy = 2

nlay = 8

f = open('mass_profile.dat','w')

f.write(str(nat*nblochy*nblochx*(nlay//4+6))+'\n')

nuc = nat*nblochx*nblochy
rndlist = np.zeros((nuc))
for i in range(nuc): 
    # Si0.2Ge0.8
    if np.random.rand()>=0.2:
        rndlist[i] = 1 # Ge
    else:
        rndlist[i] = 2 # Si



for i in range(3+nlay//4+3):
    count = 0
    for kk in range(nblochy):
        for ll in range(nblochx):
            for j in range(nat):
                f.write(str(int(2))+'\n')



f.close()

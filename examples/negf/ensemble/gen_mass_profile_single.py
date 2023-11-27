import numpy as np
import sys 
import matplotlib.pyplot as plt

# masses of binary compounds
m1 = 28.0855
m2 = 72.64 
az = 5.46616
# number of atoms per uc
nat = 4
# number of transverse unit cells
nblochx = 2 
nblochy = 2

isd = sys.argv[1] # seed number
nld = int(sys.argv[2]) # layer number
nmax = int(sys.argv[3]) # maxlay number, this determines the max cell length
nmin = 1 # maxlay number, minimum cell size

nsl = nld  #int(sys.argv[2])
numsl = np.zeros((nsl,),dtype='int')
num = np.zeros((5000,),dtype='int')

import random

random.seed(isd)
np.random.seed(seed=int(isd))
for i in range(5000):
    num[i] = nat*(random.randint(nmin,nmax))
np.save(str(isd)+'-lengthlist.npy',num)
numsl = num[0:nsl]


print('mean uc rep number:')
print(np.mean(numsl)/nat)
print(numsl[0:nld])
print('nz:',np.sum(numsl[0:nld]//nat)+6)

#(np.sum(ll[0:sl[i]])//4)+6

numz = numsl
typez = 0*numsl
for i in range(len(typez)):
    if np.mod(i,2) == 0:
        typez[i] = 1
    else:
        typez[i] = 2

nlay = np.sum(numsl[0:nld])
print(nlay*az/nat)
numsltil = np.zeros((nlay,),dtype=int)
itemp = 0
for i in range(len(typez)):
    for j in range(numz[i]):
        numsltil[itemp] = i
        itemp += 1


f = open('mass_profile.dat','w')
f.write(str(nat*nblochy*nblochx*(nlay//nat+6))+'\n')

# left contact
for i in range(3):
    for kk in range(nblochy):
        for ll in range(nblochx):
            for j in range(nat):
                f.write(str(2)+'\n')

# device region
for i in range(nlay//nat):
    for kk in range(nblochy):
        for ll in range(nblochx):
            for jj in range(nat):
#                if jj == 3 or jj==0:
#                    if np.random.rand()>0.5:
#                        f.write(str(1)+'\n')
#                    else:
#                        f.write(str(2)+'\n')
#                else:
                if typez[numsltil[i*nat+jj]]==1:
                    if np.random.rand()>0.2:
                        f.write(str(1)+'\n') # ge
                    #f.write(str(1)+'\n') # ge
                    else:
                        f.write(str(2)+'\n') # si
                else:
                    f.write(str(2)+'\n')

# right contact
for i in range(3):
    for kk in range(nblochy):
        for ll in range(nblochx):
            for j in range(nat):
                f.write(str(2)+'\n')

#import matplotlib.pyplot as plt
#plt.plot(numsl,'-o')         
#plt.show()

f.close()

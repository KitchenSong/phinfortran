import numpy as np
import matplotlib.pyplot as plt
from gen_k2d_points import Kdist,labels
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 14})

def readband(fn): 
    f = open(fn,'r')
    lines = f.readlines()
    f.close()
    n = len(lines)
    m = len(lines[0].split())
    data = np.zeros((n,m))
    for i in range(n):
        lin = lines[i].split()
        for j in range(m):
            data[i,j] = float(lin[j])
    return data

output = 'output'
intpm = 'bilinear'
shk = 0.4
ymax = 16
cmp = 'inferno'
pad = 0.02

fig,axes = plt.subplots(nrows=1,ncols=3,constrained_layout=True,figsize=(20,10))
data = readband('./'+output+'/sdosl.dat')
data = np.log(abs(data)+1e-25)
vmax = 5
vmin = -10
im = axes[0].imshow(data,interpolation=intpm,extent=[0,1,0,ymax],origin='lower',cmap=cmp,vmin=vmin,vmax=vmax)  
plt.colorbar(im,shrink=shk,ax=axes[0],format='%2.0f',pad=pad)
axes[0].set_aspect(1.0/ymax)
plt.sca(axes[0])
axes[0].tick_params(direction='out')
plt.xticks(Kdist/Kdist[-1],labels)  # Set label locations.
plt.ylabel('Frequency (Thz)')

data = readband('./'+output+'/sdosr.dat')
data = np.log(abs(data)+1e-25)
im = axes[1].imshow(data,interpolation=intpm,extent=[0,1,0,ymax],origin='lower',cmap=cmp ,vmin=vmin,vmax=vmax)
axes[1].set_aspect(1.0/ymax)
plt.colorbar(im,shrink=shk,ax=axes[1],format='%2.0f',pad=pad)
plt.sca(axes[1])
axes[1].tick_params(direction='out')
plt.xticks(Kdist/Kdist[-1],labels)  # Set label locations.
plt.ylabel('Frequency (Thz)')


data1 = readband('./'+output+'/transr_s.dat')
data = abs(data1)
vmax = 3.0
vmin = 0
im = axes[2].imshow(data,interpolation=intpm,extent=[0,1,0,ymax],origin='lower',cmap=cmp,vmin=vmin,vmax=vmax)  
axes[2].set_aspect(1.0/ymax)
plt.sca(axes[2])
axes[2].tick_params(direction='out')
plt.xticks(Kdist/Kdist[-1],labels)  # Set label locations.
plt.ylabel('Frequency (Thz)')
plt.colorbar(im,shrink=shk,ax=axes[2],pad=pad)

plt.savefig('sdos_t.png',dpi=400,transparent=True)
plt.show()


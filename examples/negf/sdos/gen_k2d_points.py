import numpy as np
import sys

labels = ['M','Y','$\Gamma$','M','X','$\Gamma$']
kpoints = "0.5 0.5 0.0 0.0 0.5 0.0 0.0 0.0 0.0 0.5 0.5 0.0 0.5 0.0 0.0 0.0 0.0 0.0"
nps = 150
kpoints = np.array(kpoints.split(),dtype=float)
nkt = len(kpoints)//3
if nkt * 3 != len(kpoints):
    print("Please provide correct k points!")
    sys.exit(1)
if nkt == 1:
    print("Please provide at least two k points!")
    sys.exit(1)

K = []
totdist = 0
Kdist = [totdist]
kpath = np.sqrt((kpoints[3]-kpoints[0])**2+(kpoints[4]-kpoints[1])**2+(kpoints[5]-kpoints[2])**2)
for i in range(0,nkt-1):
    kstart = [kpoints[i*3],kpoints[i*3+1],kpoints[i*3+2]]
    kend = [kpoints[(i+1)*3],kpoints[(i+1)*3+1],kpoints[(i+1)*3+2]]
    dist = np.sqrt((kstart[0]-kend[0])**2+(kstart[1]-kend[1])**2+(kstart[2]-kend[2])**2)
    totdist += dist
    Kdist.append(totdist)
    nps1 = int(nps*dist/kpath)
    for j in range(nps1):
        x = kstart[0]+(kend[0]-kstart[0])*j/nps1
        y = kstart[1]+(kend[1]-kstart[1])*j/nps1
        z = kstart[2]+(kend[2]-kstart[2])*j/nps1
        K.append([x,y,z])
K = np.array(K)

if __name__ == "__main__":
    # writing the k points
    f = open('kpoint.dat_2d','w')
    f.write(str(K.shape[0]+1)+'\n')
    for i in range(K.shape[0]):
        x, y, z = K[i,0], K[i,1], K[i,2]
        f.write('{0:16.13f} {1:16.13f} {2:16.13f}\n'.format(x,y,z))
    f.write('{0:16.13f} {1:16.13f} {2:16.13f}\n'.format(kpoints[-3],kpoints[-2],kpoints[-1]))
    f.close()

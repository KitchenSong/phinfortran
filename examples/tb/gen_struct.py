import numpy as np
from numpy.linalg import inv
import sys


a = float(sys.argv[1])

latvec = np.zeros((3,3))
latvec = np.array([[a, -a,  0.00000000],
[ a , a,  0.00000000],
[ 0.00000000,  0.00000000,  a*8]])
pos = np.zeros((32,4))
pos_uc = np.array([[ 0.0000000000000 , 0.0000000000000,  0.0000000000000],
                   [a/4,a/4,a/4],
                   [a/2 , 0.0000000000000, a/2],
                   [a/4, -a/4, a/4*3]])

f1 = open('input.fdf','w')
nx = 2
ny = 2
nz = 8
newcell =latvec

pp = np.zeros((3,))
natm = nz*nx*ny*4
vuc = np.zeros(newcell.shape)
vuc[:,:] = newcell[:,:]
vuc[0,:] = vuc[0,:]/float(2)
vuc[1,:] = vuc[1,:]/float(2)
vuc[2,:] = vuc[2,:]/float(nz)
natoms = 4*nx*ny*nz
ntp = 2
unitcell = np.dot(np.diag([nx,ny,nz]),vuc)

f1.write("NumberOfAtoms "+ str(natoms)+'\n'+"NumberOfSpecies "+ str(ntp)+'\n'+"LatticeConstant 1.0 Ang"+'\n'+"%block LatticeVectors\n")
for i in range(3):
    f1.write('{0:16.13f} {1:16.13f} {2:16.13f}'.format(unitcell[i,0],unitcell[i,1],unitcell[i,2])+'\n')
f1.write("%endblock LatticeVectors"+'\n'+"AtomicCoordinatesFormat  Ang"+'\n'+"%block AtomicCoordinatesAndAtomicSpecies"+'\n')



count = 0
for i in range(nz):
    for j in range(ny):
        for k in range(nx):
            for l in range(4):
                if count < natm/2:
                    ns = 2
                else:
                    ns = 1
                pp = np.dot(np.array([float(k),float(j),float(i)]),vuc)+pos_uc[l,:]
                f1.write('{0:16.13f} {1:16.13f} {2:16.13f} {3:d}'.format(pp[0],pp[1],pp[2],int(ns))+'\n')
                count += 1
f1.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

f1.write("NumberOfSpecies "+ str(ntp) +'\n')
f1.write("%block ChemicalSpeciesLabel\n")
f1.write("1 32 Ge\n")
f1.write("2 14 Si\n")
f1.write("%endblock ChemicalSpeciesLabel\n")
f1.close()

f1 = open('trans_input_uc.fdf','w')
nx = 1
ny = 1
nz = 8

pp = np.zeros((3,))
natm = nz*nx*ny*4
natoms = 4*nx*ny*nz
ntp = 2
unitcell = newcell

f1.write("NumberOfAtoms "+ str(natoms)+'\n'+"NumberOfSpecies "+ str(ntp)+'\n'+"LatticeConstant 1.0 Ang"+'\n'+"%block LatticeVectors\n")
for i in range(3):
    f1.write('{0:16.13f} {1:16.13f} {2:16.13f}'.format(vuc[i,0],vuc[i,1],unitcell[i,2])+'\n')
f1.write("%endblock LatticeVectors"+'\n'+"AtomicCoordinatesFormat  Ang"+'\n'+"%block AtomicCoordinatesAndAtomicSpecies"+'\n')



count = 0
for i in range(nz):
    for j in range(ny):
        for k in range(nx):
            for l in range(4):
                if count < natm/2:
                    ns = 2
                else:
                    ns = 1
                pp = np.dot(np.array([float(k),float(j),float(i)]),vuc)+pos_uc[l,:]
                f1.write('{0:16.13f} {1:16.13f} {2:16.13f} {3:d}'.format(pp[0],pp[1],pp[2],int(ns))+'\n')
                count += 1
f1.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

f1.write("NumberOfSpecies "+ str(ntp) +'\n')
f1.write("%block ChemicalSpeciesLabel\n")
f1.write("1 32 Ge\n")
f1.write("2 14 Si\n")
f1.write("%endblock ChemicalSpeciesLabel\n")
f1.close()

f1 = open('input_uc_left.fdf','w')
nx = 1
ny = 1
nz = 1

pp = np.zeros((3,))
natm = nz*nx*ny*4
natoms = 4*nx*ny*nz
ntp = 2
unitcell = newcell

f1.write("NumberOfAtoms "+ str(natoms)+'\n'+"NumberOfSpecies "+ str(ntp)+'\n'+"LatticeConstant 1.0 Ang"+'\n'+"%block LatticeVectors\n")
for i in range(3):
    f1.write('{0:16.13f} {1:16.13f} {2:16.13f}'.format(vuc[i,0],vuc[i,1],vuc[i,2])+'\n')
f1.write("%endblock LatticeVectors"+'\n'+"AtomicCoordinatesFormat  Ang"+'\n'+"%block AtomicCoordinatesAndAtomicSpecies"+'\n')

count = 0
for i in range(nz):
    for j in range(ny):
        for k in range(nx):
            for l in range(4):
                ns = 2
                pp = np.dot(np.array([float(k),float(j),float(i)]),vuc)+pos_uc[l,:]
                f1.write('{0:16.13f} {1:16.13f} {2:16.13f} {3:d}'.format(pp[0],pp[1],pp[2],int(ns))+'\n')
                count += 1
f1.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

f1.write("NumberOfSpecies "+ str(ntp) +'\n')
f1.write("%block ChemicalSpeciesLabel\n")
f1.write("1 32 Ge\n")
f1.write("2 14 Si\n")
f1.write("%endblock ChemicalSpeciesLabel\n")
f1.close()
f1 = open('input_uc_right.fdf','w')
nx = 1
ny = 1
nz = 1

pp = np.zeros((3,))
natm = nz*nx*ny*4
natoms = 4*nx*ny*nz
ntp = 2
unitcell = newcell

f1.write("NumberOfAtoms "+ str(natoms)+'\n'+"NumberOfSpecies "+ str(ntp)+'\n'+"LatticeConstant 1.0 Ang"+'\n'+"%block LatticeVectors\n")
for i in range(3):
    f1.write('{0:16.13f} {1:16.13f} {2:16.13f}'.format(vuc[i,0],vuc[i,1],vuc[i,2])+'\n')
f1.write("%endblock LatticeVectors"+'\n'+"AtomicCoordinatesFormat  Ang"+'\n'+"%block AtomicCoordinatesAndAtomicSpecies"+'\n')

count = 0
for i in range(nz):
    for j in range(ny):
        for k in range(nx):
            for l in range(4):
                ns = 1
                pp = np.dot(np.array([float(k),float(j),float(i)]),vuc)+pos_uc[l,:]
                f1.write('{0:16.13f} {1:16.13f} {2:16.13f} {3:d}'.format(pp[0],pp[1],pp[2],int(ns))+'\n')
                count += 1
f1.write("%endblock AtomicCoordinatesAndAtomicSpecies\n")

f1.write("NumberOfSpecies "+ str(ntp) +'\n')
f1.write("%block ChemicalSpeciesLabel\n")
f1.write("1 32 Ge\n")
f1.write("2 14 Si\n")
f1.write("%endblock ChemicalSpeciesLabel\n")
f1.close()

vpc = np.array([[ 0,a/2,a/2],
                [a/2,0, a/2],
                [a/2 , a/2,0]])
pos_uc = np.array([[0.0 ,0.0, 0.0],
                   [0.25, 0.25, 0.25],
                   [0.0, 1.0, 0.0],
                   [0.25, 1.25, -0.75]])

f1 = open('pc1.fdf','w')
nx = 1
ny = 1
nz = 1

pp = np.zeros((3,))
natm = nz*nx*ny*4
natoms = 4*nx*ny*nz
ntp = 1
unitcell = newcell

f1.write("NumberOfAtoms "+ str(natoms)+'\n'+"NumberOfSpecies "+ str(ntp)+'\n'+"LatticeConstant 1.0 Ang"+'\n'+"%block LatticeVectors\n")
for i in range(3):
    f1.write('{0:16.13f} {1:16.13f} {2:16.13f}'.format(vpc[i,0],vpc[i,1],vpc[i,2])+'\n')
f1.write("%endblock LatticeVectors"+'\n'+"%block AtomicScaledCoordinatesAndAtomicSpecies"+'\n')
count = 0
for i in range(nz):
    for j in range(ny):
        for k in range(nx):
            for l in range(4):
                ns = 2
                pp = np.dot(np.array([float(k),float(j),float(i)]),vuc)+pos_uc[l,:]
                f1.write('{0:16.13f} {1:16.13f} {2:16.13f} {3:d}'.format(pp[0],pp[1],pp[2],int(ns))+'\n')
                count += 1
f1.write("%endblock AtomicScaledCoordinatesAndAtomicSpecies\n")

f1.close()
f1 = open('pc2.fdf','w')
nx = 1
ny = 1
nz = 1

pp = np.zeros((3,))
natm = nz*nx*ny*4
natoms = 4*nx*ny*nz
ntp = 1
unitcell = newcell

f1.write("NumberOfAtoms "+ str(natoms)+'\n'+"NumberOfSpecies "+ str(ntp)+'\n'+"LatticeConstant 1.0 Ang"+'\n'+"%block LatticeVectors\n")
for i in range(3):
    f1.write('{0:16.13f} {1:16.13f} {2:16.13f}'.format(vpc[i,0],vpc[i,1],vpc[i,2])+'\n')
f1.write("%endblock LatticeVectors"+'\n'+"%block AtomicScaledCoordinatesAndAtomicSpecies"+'\n')

count = 0
for i in range(nz):
    for j in range(ny):
        for k in range(nx):
            for l in range(4):
                ns = 1
                pp = np.dot(np.array([float(k),float(j),float(i)]),vuc)+pos_uc[l,:]
                f1.write('{0:16.13f} {1:16.13f} {2:16.13f} {3:d}'.format(pp[0],pp[1],pp[2],int(ns))+'\n')
                count += 1
f1.write("%endblock AtomicScaledCoordinatesAndAtomicSpecies\n")

f1.close()

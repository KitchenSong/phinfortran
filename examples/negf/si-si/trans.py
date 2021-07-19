import numpy as np
import glob
import os
import matplotlib.pyplot as plt

fnn = 'output'
direc = (os.getcwd())
names =  glob.glob(direc+'/'+fnn+'/trl'+'*.dat')
nn = 400 #float(len(names))

keywd = 'transr_t'
names =  glob.glob(direc+'/'+fnn+'/'+keywd+'*.dat')
print(names)

for fname in names:
    f = open(fname,'r')
    lines = f.readlines()
    data = np.zeros((len(lines),))
    e = np.zeros((len(lines),))

    for i in range(len(lines)):
        e[i] = float(lines[i].split()[0])
        for j in range(len(lines[i].split())-1):
            data[i] += float(lines[i].split()[j+1])
    plt.plot(np.sqrt(e),data/nn,label=fname[:-4],ls='--',lw=2)


keywd = 'transl_s'
names =  glob.glob(direc+'/'+fnn+'/'+keywd+'*.dat')

for fname in names:
    f = open(fname,'r')
    lines = f.readlines()
    data = np.zeros((len(lines),))
    e = np.zeros((len(lines),))

    for i in range(len(lines)):
        e[i] = float(lines[i].split()[0])
        for j in range(len(lines[i].split())-1):
            data[i] += float(lines[i].split()[j+1])
    plt.plot(np.sqrt(e),(data)/nn,label=fname[:-4],alpha=0.4,ls='--')
keywd = 'transr_ns'
names =  glob.glob(direc+'/'+fnn+'/'+keywd+'*.dat')

for fname in names:
    f = open(fname,'r')
    lines = f.readlines()
    data = np.zeros((len(lines),))
    e = np.zeros((len(lines),))

    for i in range(len(lines)):
        e[i] = float(lines[i].split()[0])
        for j in range(len(lines[i].split())-1):
            data[i] += float(lines[i].split()[j+1])
    plt.plot(np.sqrt(e),data/nn,label=fname[:-4])


#np.save('d0.npy',data/nn)
#d0 = np.load('d0.npy')
#plt.plot(np.sqrt(e),d0)




plt.legend()
#plt.savefig('phtr.pdf',dpi=200)
plt.show()



import numpy as np
import glob
import os
import matplotlib.pyplot as plt

nk = 49

A = (5.46616*np.sqrt(2))**2*1.0e-20

kb = 1.380649e-23
hbar = 1.0545718e-34
thz2rads = 1e12*np.pi*2

def condt(w,T,trans):
    omega = w
    itd = w*0.0
    expf = w*0.0
    dfdw = w*0.0
    f = w*0.0

    expf = np.exp(hbar*omega/kb/T)

    dfdw = expf/(expf-1)**2*hbar*omega/kb/T**2
    f = 1/(expf-1)


    itd = 1/A/2/np.pi*(hbar*w*dfdw*trans)

    cd = (np.trapz(itd,x=w))
    return cd


def pltran(temp):
    fnn = './output'
    direc = (os.getcwd())
    keywd = 'transr_t'
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

    return condt(np.sqrt(e)*thz2rads,temp,data/nk)


plt.figure(figsize=(6,4))
templist = np.array([300, 400, 500,600])
templist = np.linspace(10,700,100)
c = np.zeros(templist.shape)
count = 0
for count,tt in enumerate(templist): 
    c[count] = pltran(tt)
    count+=1
plt.plot(templist,c)
plt.xlabel('Temperature (K)')
plt.ylabel(r'Two-probe conductance ($W/K-m^2$)')
#l = (llist-6)*4*1.366540979010257
#plt.plot(l,clist*l*1.0e-10,marker='o',alpha=0.8,markersize=10,label='T='+str(temp)+'K')
#plt.xlabel('Length [A]')
#plt.ylim([0,0.6])
#plt.xlim([0,10])
#plt.ylim([0,0.5])
#plt.figure(figsize=(6,4))

#pltran(50)

#plt.savefig('phtr.pdf',dpi=200)
plt.show()



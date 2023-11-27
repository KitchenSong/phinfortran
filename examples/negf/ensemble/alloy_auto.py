import os
import numpy as np
import time
import psutil
import re

def return_n_proc():
    # return number of actual running processes
    process_names = [proc.name() for proc in psutil.process_iter()]
    count = 0
    for i in range(len(process_names)):
        if re.findall(r"phinfortran",process_names[i]):
            #print(process_names[i])
            count += 1
    return count


def return_job_status(direct):
    if os.path.isfile(direct+'/transl_tot.dat'):
        return True
    else:
        return False

def return_all_job_status(sd,nz):
    finished = []
    running = []
    unfinished = []
    count = 0
    for i in range(len(nz)):
        if return_job_status(str(sd)+'/'+str(nz[i])+'/'+'output'):
            finished.append(i)
            count = count + 1
        elif os.path.isdir(str(sd)+'/'+str(nz[i])+'/output'):
            running.append(i)
        else:
            unfinished.append(i)
    return count,finished,running,unfinished

def doing_calculation(ncore,lmax,sd,ll,sl):
    nat = 4 # uc atom number
    # number of calculations
    ncal = len(sl)
    nz = np.zeros(sl.shape,dtype=int)

    nd = sl*1+6

    for i in range(len(sl)):
        nz[i] = (np.sum(ll[0:sl[i]])//nat)+6

    count,f,rn,uf = return_all_job_status(sd,nz)
    #print(nz)
    #print(nz[17])
    #print(count,uf,rn)
    start_time = time.monotonic()
    tag = 0

    while return_all_job_status(sd,nz)[0] <= ncal:
        count,f,rn,uf = return_all_job_status(sd,nz)
        while len(rn)<ncore and len(uf)>0:
            iuf = uf[0]
            os.system('sed -i "/nz=/cnz='+str(nz[iuf])+'" ./running.sh') 
            os.system('sed -i "/nsl=/cnsl='+str(sl[iuf])+'" ./running.sh') 
            os.system('sed -i "/sd=/csd='+str(sd)+'" ./running.sh') 
            os.system('sed -i "/lmax=/clmax='+str(lmax)+'" ./running.sh') 
            os.system('sh running.sh')
            time.sleep(3)
            print('starting 1 job ...')
            count,f,rn,uf = return_all_job_status(sd,nz)
        time.sleep(5)
        count,f,rn,uf = return_all_job_status(sd,nz)
        if len(rn) == 0 and len(uf) ==0:
            print(round((time.monotonic() - start_time)/3600,3),'h')
            tag = 1
            break
        if len(rn)>0 and len(rn)<ncore and return_n_proc()==0:
            print(round((time.monotonic() - start_time)/3600,3),'h')
            print('there are still unfinished tasks')
            tag = 2
            break
    return tag

ncore = 1 # num of available cpus
lmax = 8
sd = 0

# dummy run
os.system('python gen_mass_profile_single.py '+ str(sd) +' 20 '+str(lmax)+';')

ll = np.load(str(sd)+'-lengthlist.npy')


N = 22
sstart = 3
sl = np.arange(sstart,2*N+3,2)

doing_calculation(ncore,lmax,sd,ll,sl)



#isl = [12,14]
#sl = sl[isl]
    
#ncal = len(sl)
#nz = np.zeros(sl.shape,dtype=int)
#nd = sl*1+6
#for i in range(len(sl)):
#        nz[i] = (np.sum(ll[0:sl[i]])//4)+6

#count,f,rn,uf = return_all_job_status(sd,nz)
#nz = np.zeros(sl.shape,dtype=int)
#print(f,rn,uf)
# in this step, we need to resubmit those tasks (identified as running but not actually running) that are not properly finished.



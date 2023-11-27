import psutil
import re
from psutil._common import bytes2human


def pprint_ntuple(nt):
    for name in nt._fields:
        value = getattr(nt, name)
        if name != 'percent':
            value = bytes2human(value)
        print('%-10s : %7s' % (name.capitalize(), value))

def return_n_proc():
    process_names = [proc.name() for proc in psutil.process_iter()]
    print('MEMORY\n------')
    pprint_ntuple(psutil.virtual_memory())
    print('\nSWAP\n----')
    pprint_ntuple(psutil.swap_memory())

    #psutil.virtual_memory
    count = 0
    for i in range(len(process_names)):
        if re.findall(r"phinfortran",process_names[i]):
            #print(process_names[i])
            count += 1
    return count
if return_n_proc()==0:
    print(return_n_proc())


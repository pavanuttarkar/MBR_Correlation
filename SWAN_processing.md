### Download data from drive

```
https://drive.google.com/drive/u/1/folders/1LXoIMgWsf9PmoJQxtfDnxJ9_x-qqyGVo
```


```
from numpy import *
from matplotlib.pyplot import *
```


```
dt                  =    np.dtype([('header', 'S8'),
                                  ('Source', 'S10'),
                                  ('Attenuator_1', '>u1'),
                                  ('Attenuator_2', '>u1'),
                                  ('Attenuator_3', '>u1'),
                                  ('Attenuator_4', '>u1'), 
                                  ('LO', '>u2'),
                                  ('FPGA', '>u2'),
                                  ('GPS', '>u2'),
                                  ('Packet', '>u4'),
                                  ('data', '>i1', 1024)])
```


### Read data using memcopy


`data = np.memmap('ch01_test.mbr',  dtype = dt, mode = 'c')`

### Example output

```
f2py3 -c  internal_standalone_WORKING_with_RFI.f90 -L/usr/lib/x86_64-linux-gnu/ -lfftw3 -lm -m internal_loop_RFI -DF2PY_REPORT_ON_ARRAY_COPY=1 --f90flags='-fcheck=bounds -frepack-arrays -C' --fcompiler='gfo
rtran'
f2py3 -c  internal_standalone_WORKING_with_RFI_for_single_file_dev.f90 -L/usr/lib/x86_64-linux-gnu/ -lfftw3 -lm -m single_file_combined_RFInoRFI -DF2PY_REPORT_ON_ARRAY_COPY=1 --f90flags='-fcheck=bounds -fre
pack-arrays -C'
mkdir OBJECT_FILE
cp -aur *.so ../OBJECT_FILE/.
```

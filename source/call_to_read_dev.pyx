cimport numpy as np
import numpy as np
import binascii
from os.path import *
import time
from _header_Fring_cy import *
import datetime
from datetime import datetime
import internal_loop2 
import _header_Fring_cy
import _header_geometric_cy
import multiprocessing
import mmap
import glob
from scipy.stats import linregress
from sympy import S, symbols
import sympy
from os.path import getsize
#cython: profile=True




###############################################
###      ## ######  #####   ####    ##  ###### # 
# ## ## ##  ######  ##  #   #####  ##   ###### #
#  # ## #   ##  ##  #####   ##  ####      ##   #
#   #  #    ##  ##  ##  #   ##  ####    ###### #
#   #  #    ##  ##  ##  #   ##  ####    ###### #
###############################################

#
#fromfile is highly unpredictable interms of reading time!!
#
from cpython cimport Py_buffer
from cpython.buffer cimport PyBUF_SIMPLE, PyBUF_WRITEABLE
from libcpp.vector cimport vector


# Adding debug flag
debug = False

# Exception
class LOMismatch(Exception):
    """
    Exception Raised when LO is different.
    """
    def __init__(self, message="LO Mismatch"):
        self.message = message
        super().__init__(self.message)
    pass

#Address location
ADDR_LOC = "./SAMPLING_INFO/"

cdef class SimplestBuffer:
    cdef:
        vector[char] buf   # We're using vector from C++ to manage memory allocation
     
    # in cython, methods defined with 'def' are slower but accessible from python
    # extend will add bytes to our internal memory
    def extend(self, input_bytes):
        self.add_bytes(input_bytes, len(input_bytes))   
    
    # methods defined with 'cdef' may be faster but are accessible only from cython
    cdef add_bytes(self, char *b, int num_bytes):  
        self.buf.insert(self.buf.end(), b, b + num_bytes)
    
    def __getbuffer__(self, Py_buffer *buffer, int flags):
        # if the requested buffer type is not PyBUF_SIMPLE then error out
        # we will allow either readonly or writeable buffers however
        if flags != PyBUF_SIMPLE and flags != PyBUF_SIMPLE | PyBUF_WRITEABLE:
            raise BufferError
            
        buffer.buf = &self.buf[0]            # points to our buffer memory
        buffer.format = NULL                 # NULL format means bytes 
        buffer.internal = NULL               # this is for our own use if needed
        buffer.itemsize = 1                  # size in bytes of a single element   
        buffer.len = self.buf.size()
        buffer.ndim = 1
        buffer.obj = self
        buffer.readonly = not (flags & PyBUF_WRITEABLE)
        buffer.shape = NULL                  # none of shapes, strides or suboffsets
        buffer.strides = NULL                # are used with PyBUF_SIMPLE
        buffer.suboffsets = NULL    

    # the buffer protocol requires this method
    def __releasebuffer__(self, Py_buffer *buffer):
        pass       

cdef linear_fit(a,b):
        z                       =       np.polyfit(a,b,1)
        print("ATTENTION!!")
        print(z[0], z[1])
        r                       =       np.poly1d(z)
        a_new                   =       np.linspace(a[0],a[-1],len(a))
        b_new                   =       r(a_new)
        x = symbols("*x")#Change of * made
        poly = sum(S("{:6.2f}".format(v))*x**i for i, v in enumerate(z[::-1]))
        eq_latex = sympy.printing.latex(poly)
        return  a_new,b_new, eq_latex    

cpdef fil_list(fil1):
    f1_1  =   open(fil1, 'r').readlines()
    f1  =   ''.join(f1_1).split('>>>')[0]
    f1  =   f1.split('\n')[:-1]

    f2  =   ''.join(f1_1).split('>>>')[1]
    f2  =   f2.split('\n')[1:-1]
    

    cdef int i  =   0

    cdef list line =   []
    cdef list slope1  =  []  
    cdef list slope2  =  []
    cdef list inter1  =  []
    cdef list inter2  =  []


    for i in range(len(f1)):
        line.append(decrypy_get_gpssync(f1[i], f2[i], 60))  # 60 is insignificant ason 23 DEC 2020.
        slope1.append(line[i][0].slope)
        slope2.append(line[i][1].slope)
        inter1.append(line[i][0].intercept)
        inter2.append(line[i][1].intercept)
    slope   =   (np.mean(slope1)+np.mean(slope2))/2
    inter   =   np.mean(inter1)/2+np.mean(inter2)/2
    for i in range(len(f1)):
        tf1    =   open(str(ADDR_LOC + "Info_on_straight_line"+str(f1[i].split('/')[-1])), 'rw')
        tf2    =   open(str(ADDR_LOC + "Info_on_straight_line"+str(f2[i].split('/')[-1])), 'rw')
        temp1  =   tf1.readlines()
        temp2  =   tf2.readlines()
        if(inter1 >0):
            temp1[2] =   str(slope)+'*x'+str(inter)+'\n'
            temp2[2] =   str(slope)+'*x'+str(inter)+'\n'
        else:
            temp1[2] =   str(slope)+'*x+'+str(slope)
            temp2[2] =   str(slope)+'*x+'+str(slope)
        tf1.close()
        tf2.close()
        tf1    =   open(str(ADDR_LOC + "Info_on_straight_line"+str(f1[i].split('/')[-1])), 'w+')
        tf2    =   open(str(ADDR_LOC + "Info_on_straight_line"+str(f2[i].split('/')[-1])), 'w+')

        tf1.writelines(temp1)
        tf2.writelines(temp2)
        tf1.close()
        tf2.close()
    return 0

cpdef tuple decrypy_get_gpssync(file_name, file_name1, avg, s_mbr='000.mbr'):
    '''
    Getting GPS sync data..

    Args:
        file_name (str): File Name of the 1st File.
        file_name1 (str): File Name of the 2nd File.
        avg (int): # Packets to Average
        s_mbr (str) [000.mbr]: Series
    '''
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    comf     = np.memmap(file_name,  dtype = dt, mode = 'c')
    comf1    = np.memmap(file_name1, dtype = dt, mode = 'c')
    cdef list accblipg1   =   []
    cdef list accblipp1   =   []
    cdef list accblipg2   =   []
    cdef list accblipp2   =   []
    cdef list blipp1       =   []
    cdef list blipp2       =   []
    cdef list blipg1       =   []
    cdef list blipg2       =   []

    cdef int fact        =   1
    cdef int i           =   0
    cdef str series      = file_name[-7:-4]+'.mbr'
    cdef long long int len1        = 0#comf['Packer'][-1] - comf['Packet'][0] 
    cdef long long int len2        = 0#comf1['Packet'][-1] - comf1['Packet'][0]
    
    for i in range(len(comf)):
        temp = comf['FPGA'][i] & fact
        if(temp == 0):
            accblipp1.append(comf['Packet'][i])
            accblipg1.append(comf['GPS'][i])
    for i in range(len(comf1)):
        temp=comf1['FPGA'][i] & fact
        if(temp == 0):
            accblipp2.append(comf1['Packet'][i])
            accblipg2.append(comf1['GPS'][i])
    slp1    =    linregress(accblipg1, accblipp1)
    slp2    =    linregress(accblipg2, accblipp2)
    slp_1   =    linear_fit(accblipg1, accblipp1)
    slp_2   =    linear_fit(accblipg2, accblipp2)
    
    x   =   max(comf['GPS'][0], comf1['GPS'][0])+1

    fillist         =               np.loadtxt('DEPEND/files_naming.txt', dtype=str)
    series_code     =               int(np.where(fillist==series)[-1][0])#fillist.index(series)+1
    print('\n\n\n\n\n') 
    print(series_code)

    fileeq1         =               open(str(ADDR_LOC + "Info_on_straight_line"+str(file_name.split('/')[-1])), 'w+')
    fileeq2         =               open(str(ADDR_LOC + "Info_on_straight_line"+str(file_name1.split('/')[-1])), 'w+')

    fileeq1.write('First GPS Value == '+str(comf['GPS'][10])+'\n'+str(comf['Packet'][10])+'\n')
    fileeq2.write('First GPS Value == '+str(comf1['GPS'][10])+'\n'+str(comf1['Packet'][10])+'\n')

    ctr     =   1
    temp    =   0

    if(int(float(eval(slp_1[-1]))) < 0):
        temp    =   1
    len1    =   int(float(eval(slp_1[-1])))-comf['Packet'][0]#len(comf)   -   int(float(eval(slp_1[-1])))
    len2    =   int(float(eval(slp_2[-1])))-comf1['Packet'][0]#len(comf1)  -   int(float(eval(slp_2[-1])))
    print(x, len(comf), str(int(float(eval(slp_1[-1])))), len(comf[int(float(eval(slp_1[-1]))):]))
    print(x, len(comf1), str(int(float(eval(slp_2[-1])))), len(comf1[int(float(eval(slp_1[-1]))):]))
    rem_len =   abs(len(comf[len1:]) - len(comf1[len2:]))#abs(len(comf[int(float(eval(slp_1[-1]))):])-len(comf1[int(float(eval(slp_1[-1]))):]))

    print('######################################')
    print(rem_len)
    print('######################################')
    if(len1 > len2):#if(str(eval(slp_1[-1])) > str(eval(slp_2[-1])) and str(eval(slp_2[-1])) > 0):
        fileeq1.write(str(slp_2[-1])+'\n')
        fileeq2.write(str(slp_2[-1])+'\n')
        
        file_to_write       =   file_name.split('_')
        file_to_write_fil   =   file_name.split('_')
        file_to_write[-1]   =   fillist[int(series_code)+1]
        file_to_write_fil[-1]=  fillist[series_code]
        file_to_write       =   '_'.join(file_to_write)
        print(file_to_write)
        file_to_write_fil   =   '_'.join(file_to_write_fil)

        filrem          =               open(ADDR_LOC + str(file_to_write_fil.split('/')[-1])+'_rem.data', 'w+')
        filrem.write(str(rem_len)+'\n')#str(eval(slp_2[-1])-comf1['Packet'][0]-eval(slp_1[-1])+comf['Packet'][0])+'\n')
        filrem.write(str(file_to_write)+'\n')

        print(fillist[int(series_code)+1])
        fileeq2.write('N,0\n')
        fileeq1.write('Y,'+str(file_to_write)+','+str(rem_len)+'\n')

    else:
        fileeq1.write(str(slp_1[-1])+'\n')
        fileeq2.write(str(slp_1[-1])+'\n')


        file_to_write       =   file_name1.split('_')
        file_to_write_fil   =   file_name1.split('_')
        file_to_write[-1]   =   fillist[int(series_code)+1]
        file_to_write_fil[-1]=  fillist[series_code]
        file_to_write       =   '_'.join(file_to_write)
        print('=============')
        print(file_to_write)
        print('=============')
        file_to_write_fil   =   '_'.join(file_to_write_fil)


        filrem          =               open(ADDR_LOC + str(file_to_write_fil.split('/')[-1])+'_rem.data', 'w+')
        filrem.write(str(rem_len)+'\n')#str(eval(slp_1[-1])-comf['Packet'][0]-eval(slp_2[-1])+comf1['Packet'][0])+'\n')
        filrem.write(str(file_to_write)+'\n')

        fileeq1.write('N,0\n')
        fileeq2.write('Y,'+str(file_to_write)+','+str(rem_len)+'\n')

    fileeq1.write(str(comf['GPS'][0])+','+str(comf['GPS'][-1]))
    fileeq2.write(str(comf1['GPS'][-1])+',' + str(comf1['GPS'][-1]))
        
    print(str(slp2.slope)+'*x'+str(slp2.intercept))
    print(str(slp1.slope)+'*x'+str(slp1.intercept))

    #Generating RFI mask#
    if(series == s_mbr):
        comf_X, comf_Y, comf1_X, comf1_Y, LO = call_to_read(file_name, file_name1, 60, '1', np.zeros((30)))
        le      = len(comf_X)/(int(int(256)*2)*int(avg))
        le1     = len(comf1_X)/(int(int(256)*2)*int(avg))
        le      =   min(le, le1)
        creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 = internal_loop2.external_loop(comf_X, comf_Y, comf1_X, comf1_Y, avg, le, 255)
        RFI1    =   RFI_Reject(creal1, 3, 256, 60)
        RFI2    =   RFI_Reject(creal2, 3, 256, 60)
        RFI3    =   RFI_Reject(creal3, 3, 256, 60)
        RFI4    =   RFI_Reject(creal4, 3, 256, 60)
        
        #Writing RFI masks#
        np.savetxt(str(ADDR_LOC + "RFI_mask_X_"+str(file_name.split('/')[-1])), RFI1)
        np.savetxt(str(ADDR_LOC + "RFI_mask_Y_"+str(file_name.split('/')[-1])), RFI2)
        np.savetxt(str(ADDR_LOC + "RFI_mask_X_"+str(file_name1.split('/')[-1])), RFI3)
        np.savetxt(str(ADDR_LOC + "RFI_mask_Y_"+str(file_name1.split('/')[-1])), RFI4)

    return slp1, slp2, slp_1, slp_2#, RFI_mask1, RFI_mask2, RFI_mask3, RFI_mask4

cpdef read_RFI(file_name, file_name1):
    """
    Function to read saved RFI Matrix.

    Args:
        file_name (str): 1st File Name
        file_name1 (str): 2nd File Name 
    
    Returns: 
        - **RFI_mask1** array(int) RFI Mask X1 
        - **RFI_mask2** array(int) RFI Mask Y1 
        - **RFI_mask3** array(int) RFI Mask X2
        - **RFI_mask4** array(int) RFI Mask Y2
    """

    fil_000         =               file_name[:-7]+'000.mbr'
    fil1_000        =               file_name1[:-7]+'000.mbr'

    RFI_mask1       =               np.loadtxt(str(ADDR_LOC + "RFI_mask_X_"+str(fil_000.split('/')[-1])))#+fillist[0])
    RFI_mask2       =               np.loadtxt(str(ADDR_LOC + "RFI_mask_Y_"+str(fil_000.split('/')[-1])))
    RFI_mask3       =               np.loadtxt(str(ADDR_LOC + "RFI_mask_X_"+str(fil1_000.split('/')[-1])))#+fillist[0])
    RFI_mask4       =               np.loadtxt(str(ADDR_LOC + "RFI_mask_Y_"+str(fil1_000.split('/')[-1])))

    return RFI_mask1, RFI_mask2, RFI_mask3, RFI_mask4

cpdef tuple decrypy_file_new_SWAN_onhold(file_name, file_name1, ch):
    '''
    ch should be the channel number of first file..
    Takes the read file and sorts the X and Y polarizartion in the file into
    comf_X, comf1_X, comf_Y, comf1_Y.

    Args:
        file_name (str): 1st File Name
        file_name1 (str): 2nd File Name
        ch (int): Channel No.
    
    Returns:
        - **tempcomf_X** array(float)
        - **tempcomf_Y** array(float)
        - **tempcomf1_X** array(float)
        - **tempcomf1_Y** array(float)
        - **LO1** (int)
    '''

    #If available get the earlier data set#
    cdef str series      = file_name[-7:-4]
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    if(series   !=  '000'):
        rem     =   open(ADDR_LOC + 'Info_on_straight_line'+str(file_name.split('/')[-1]), 'r').readlines()   #open('/data/swan/pavan/SAMPLING_INFO/'+str(file_name.split('/')[-1])+'_rem.data', 'r').readlines()
        print(rem)
        rem_l   =   rem[-2].split(',')
        
        if(rem_l[0] =='Y'):
            comf_rem             = np.memmap(rem_l[1], dtype = dt, mode='c')
            comf_rem             = comf_rem[:int(float(rem_l[2][:-1]))]
        
        else:
            rem     =   open(ADDR_LOC + 'Info_on_straight_line'+str(file_name1.split('/')[-1]), 'r').readlines()
            print(rem)
            rem_l   =   rem[-2].split(',')
            comf_rem             = np.memmap(rem_l[1], dtype = dt, mode='c')
            comf_rem             = comf_rem[:int(float(rem_l[2][:-1]))]
     
    cdef int LO1, LO2

    cdef double timea    = time.time()
    cdef str fil_tag     = file_name.split('_')[-4]
    cdef int    i        = 0
    #op1                  = open(file_name, 'rb')
    #op2                  = open(file_name1, 'rb')
    print('Caching memory, for faster read..')    
    #os.system('./comp '+str(file_name)+' '+str(file_name1))
    comf     = np.memmap(file_name,  dtype = dt, mode = 'c')
    comf1    = np.memmap(file_name1, dtype = dt, mode = 'c')
    print('Time required to read..'+str(time.time()-timea))
    print(file_name)

    if(series != '000' and file_name[-39:-35] == rem_l[1][-39:-35]):
        print('\n\n\n\n\n\n\n\nIn 1\n\n\n\n')
        comf    =   np.hstack((comf, comf_rem))
    elif(series != '000' and file_name1[-39:-35] == rem_l[1][-39:-35]):
        print('\n\n\n\n\n\n\n\nIn 2\n\n\n\n')
        comf1   =   np.hstack((comf1, comf_rem))
    print(len(comf), len(comf1))
    LO1                  =   comf['LO'][100]
    LO2                  =   comf1['LO'][100]

    if(LO1!=LO2):
        print('LO1 = ' +str(LO1)+'\n')
        print('LO2 = ' +str(LO2)+'\n')
        #raise RuntimeError ('Both LOs are different!')
        raise LOMismatch

    cdef long long int templen    =   0
    cdef long long int templen1   =   0
    
    ft      =   np.dtype('>i1')
    st      =   max(comf['GPS'][0], comf1['GPS'][0])
    op      =   min(comf['GPS'][-1], comf1['GPS'][-1])

    len1_fil1    =   comf['Packet'][10]
    len1_fil2    =   comf1['Packet'][10]
    len2_fil1    =   comf['Packet'][-1]
    len2_fil2    =   comf1['Packet'][-1]

    start_point1 =   len(comf['data']) -10
    start_point2 =   len(comf1['data'])-10
    templen =   len2_fil1  -len1_fil1 + 10
    templen1=   len2_fil2  -len1_fil2 + 10
    
    print('len2_fil1-len1_fil1 '+str(len2_fil1)+'-'+str(len1_fil1)+'='+str(templen))
    print('len2_fil2-len1_fil2 '+str(len2_fil2)+'-'+str(len1_fil2)+'='+str(templen1))


    print(templen, templen1)
    cdef np.ndarray  lin1       =   comf['Packet'][10:]  -len1_fil1#comf['Packet'][100]
    cdef np.ndarray  lin2       =   comf1['Packet'][10:] -len1_fil2# comf1['Packet'][100]
    cdef np.ndarray tempcomf    =   np.zeros((templen, 1024), dtype = np.int8)
    cdef np.ndarray tempcomf1   =   np.zeros((templen1, 1024), dtype = np.int8)
    tempcomf[lin1]    =   np.memmap.copy(comf['data'][10:])#data_temp#np.frombuffer(data_temp, dtype=ft, count=1024)
    
    print('Time required to decrypt one file..'+str(time.time()-timea))
    tempcomf1[lin2]   =   np.memmap.copy(comf1['data'][10:])#data_temp#np.frombuffer(data_temp, dtype=ft, count=1024)
    tempcomf    = tempcomf.ravel()
    tempcomf1   = tempcomf1.ravel()

    tempcomf_X  			=   np.array(tempcomf[1::2], order = 'F')
    tempcomf_Y  			=   np.array(tempcomf[0::2],  order = 'F')
    tempcomf1_X 			=   np.array(tempcomf1[1::2], order = 'F')
    tempcomf1_Y 			=   np.array(tempcomf1[0::2], order = 'F')
   
    print('Time required to decrypt both files..'+str(time.time()-timea))
    #GPS_len =   [op, st]
    return tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, LO1

cpdef tuple decrypy_file_new_SWAN_onhold_hold(file_name, file_name1, avg):
    """
    Takes the read file and sorts the X and Y polarizartion in the file into
    comf_X, comf1_X, comf_Y, comf1_Y.

    Args:
        file_name (str): 1st File Name
        file_name1 (str): 2nd File Name
        avg (int): # Packets to Average
    
    Returns:
        - **tempcomf_X** array(int)
        - **tempcomf_Y** array(int)
        - **tempcomf1_X** array(int)
        - **tempcomf1_Y** array(int)
        - **LO1** (int) LO
    """

    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    cdef int LO1, LO2
    cdef str series      = file_name[-7:-4]
    cdef double timea    = time.time()
    cdef str fil_tag     = file_name.split('_')[-4]
    cdef int    i        = 0

    comf     = np.memmap(file_name,  dtype = dt, mode = 'c')
    comf1    = np.memmap(file_name1, dtype = dt, mode = 'c')

    print('Time required to read..'+str(time.time()-timea))

    LO1                  =   comf['LO'][10]
    LO2                  =   comf1['LO'][10]

    if(LO1!=LO2):
        print('Both LOs are different!')
        #raise RuntimeError ('Both LOs are different!')
        raise LOMismatch

    cdef long long int templen    =   0
    cdef long long int templen1   =   0

    ft      =   np.dtype('>i1')

    len1_fil1    =   comf['Packet'][10]
    len1_fil2    =   comf1['Packet'][10]
    len2_fil1    =   comf['Packet'][-1]
    len2_fil2    =   comf1['Packet'][-1]

    start_point1 =   len(comf['data']) -10
    start_point2 =   len(comf1['data'])-10
    templen =   len2_fil1  -len1_fil1 + 10
    templen1=   len2_fil2  -len1_fil2 + 10

    print('\n')
    print('len1_fil1, len1_fil2, len2_fil1, len2_fil2')
    print(len1_fil1, len1_fil2, len2_fil1, len2_fil2)

    cdef np.ndarray  lin1       =   comf['Packet'][10:] - len1_fil1#comf['Packet'][100]
    cdef np.ndarray  lin2       =   comf1['Packet'][10:] -len1_fil2# comf1['Packet'][100]
    cdef np.ndarray tempcomf    =   np.zeros((templen, 1024), dtype = np.int8)
    cdef np.ndarray tempcomf1   =   np.zeros((templen1, 1024), dtype = np.int8)

    tempcomf[lin1]    =   np.memmap.copy(comf['data'][10:])#data_temp#np.frombuffer(data_temp, dtype=ft, count=1024)

    print('Time required to decrypt one file..'+str(time.time()-timea))

    tempcomf1[lin2]   =   np.memmap.copy(comf1['data'][10:])#data_temp#np.frombuffer(data_temp, dtype=ft, count=1024)
    tempcomf    = tempcomf.ravel()
    tempcomf1   = tempcomf1.ravel()

    tempcomf_X  			=   np.array(tempcomf[1::2], order = 'F')
    tempcomf_Y  			=   np.array(tempcomf[0::2],  order = 'F')
    tempcomf1_X 			=   np.array(tempcomf1[1::2], order = 'F')
    tempcomf1_Y 			=   np.array(tempcomf1[0::2], order = 'F')

    print('Time required to decrypt both files..'+str(time.time()-timea))
    #Tryring to increase speed of this module#
    return tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, LO1

cpdef fix_time(sec, mi, hr, day):
    """
    Function to Time

    Args:
        sec (int): Seconds
        mi (int): Minute
        hr (int): Hour
        day (int): Day

    Returns:
        - **sece** (int) Seconds
        - **minu** (int) Minute
        - **hour** (int) Hour
        - **day** (int) Day
        
    """
    cdef int trk =   0
    cdef int fact=   1
    
    hour    =   int(sec/3600.0)+hr
    minu    =   (sec/3600.0%1)*60+mi
    sece    =   (minu%1)*60
    
    if(sece > 59):
        while(fact):
            print('In sec') 
            sece =   sece-60
            if(sece<60 and sece > 0):
                fact=0
            trk =   trk+1   
            mi  =   mi+1
        fact=1
        trk=0

    if(minu>59):   
        while(fact):
            
            print('In min')
            minu =   minu-60
            if(minu<60 and minu > 0):
                fact=0
            trk =   trk+1
            hr  =   hr+1
        fact=1
        trk=0
    
    if(hour>23):
        while(fact):
            hour =   hour-24
            if(hour<24 and hour > 0):
                fact=0
            trk =   trk+1
            day =  day+1
    
        fact=1
        trk=0
    return sece, minu, hour, day

cpdef shift_geo(np.ndarray comf_X, np.ndarray comf_Y, np.ndarray comf1_X, np.ndarray comf1_Y, np.ndarray delay):
    """
    Compensate for geometric delay.

    Args:
        comf_X
        comf_Y
        comf1_X
        comf1_Y
        delay
    
    Returns:
        - **comf_X**
        - **comf_Y**
        - **comf1_X**
        - **comf1_Y**

    """
    cdef i                  =   0
    cdef np.ndarray loc     =   np.linspace(0, len(comf_X), len(delay))

    cdef double delay_r1           =    0
    cdef double delay_trc          =    0
    for i in range(len(loc)-1):
        if(delay_r1 > 1.0):
            print('In delay_resv')
            print(int(delay_r1))
            #delay_trc                 =   loc[i]  +   sum(delay[:i])
            comf_X                    =   np.insert(comf_X, int(loc[i]) , np.zeros((int(delay_r1)), dtype=np.int8))    
            comf_Y                    =   np.insert(comf_Y, int(loc[i]) , np.zeros(int((delay_r1)), dtype=np.int8))
            print(len(comf_X), len(comf_Y))

            delay_r1   =   delay_r1-int(delay_r1)
        else:
            delay_r1         =  delay_r1    +    delay[i] 
    return comf_X, comf_Y, comf1_X, comf1_Y


"""
cdef post_gps_test(part_X, part_Y, part1_X, part1_Y):
    '''
        Considering only 1/4th of the file, we check for the correlation between the files, thereby deriving the
        delay parameters, to optimize delay calculation..
    '''
    
    return 
"""


cpdef call_to_read(str file_name, str file_name1, avg, sysargv, delay):#, number, number1):
    """
    Function to read and apply geo shift.

    Args:
        file_name (str): File 1 Name
        file_name1 (str): File 2 Name 
        avg (int): # Packets to average
        sysargv: sysargv
        delay: Delay
    
    Returns:
        - **comf_X** array(float)
        - **comf_Y** array(float)
        - **comf1_X** array(float)
        - **comf1_Y** array(float)
        - **LO** (int)
    """

    #Binary file structure#
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    timea    = time.time()

    if(sysargv==str(1)):
        comf_X, comf_Y, comf1_X, comf1_Y,LO        =       decrypy_file_new_SWAN_onhold(file_name, file_name1, avg)
        print('__________________________________________________')
        print(len(comf_X), len(comf_Y), len(comf1_X), len(comf1_Y))
        print('__________________________________________________')
        comf_X, comf_Y, comf1_X, comf1_Y        =       shift_geo(comf_X, comf_Y, comf1_X, comf1_Y, delay)
    if(sysargv==str(2)):
        print('MBRDSP without packet compensation..')
    return comf_X, comf_Y, comf1_X, comf1_Y, LO

def gen_RFI_matrix(str file_name, str file_name1):
    '''
        .. warning::
            Should be run only one, not for every file..time constraint :/
        
        Function to generate RFI Matrix

        Args:
            file_name (str): File 1 Name
            file_name1 (str): File 2 Name

        Returns:
            - **RFI_X1** array(int) X1 RFI Reject Array 
            - **RFI_Y1** array(int) Y1 RFI Reject Array
            - **RFI_X2** array(int) X2 RFI Reject Array
            - **RFI_Y2** array(int) Y2 RFI Reject Array 
    '''
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])

    comf_X, comf_Y, comf1_X, comf1_Y,LO        =       decrypy_file_new_SWAN_onhold(file_name, file_name1, 60)

    creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 = internal_loop_RFI.external_loop(comf_X_1, comf_Y_1, comf1_X_1, comf1_Y_1, avg, le, 255)
    RFI_X1  =   RFI_Reject(creal1, 3, 256, 60)
    RFI_Y1  =   RFI_Reject(creal2, 3, 256, 60)
    RFI_X2  =   RFI_Reject(creal3, 3, 256, 60)
    RFI_Y2  =   RFI_Reject(creal4, 3, 256, 60)
    
    return RFI_X1, RFI_Y1, RFI_X2, RFI_Y2

cdef np.ndarray rms(a,meana, fftsize):
        """
        Function to Calcualte RMS for Efficiency Calculation.

        Args:
            a array(float): Input Array
            meana (float): Mean of the Input Array
            fftsize (int): FFT Size
        
        Returns:
            - **rmsa** array(float) RMS of Input
        """

        cdef np.ndarray rmsa       =       np.zeros(int(fftsize), dtype='double')
        cdef int i              =       0
        for i in range(int(fftsize)):
                #for j in range(len(a[0])):
                 rmsa[i] = np.sqrt(np.mean((a[i] - meana[i])**2))#    rmsa[i] += (a[i][j] - meana[i])**2
                #rmsa[i]=sqrt(rmsa[i]/len(a[0]))
        return rmsa

cpdef np.ndarray RFI_Reject(spec, sig, fftsize, avg):
    '''
    RFI Rejection module,
    
    Args:      
        spec, sigma deviation, fftsize
    
    Returns:
        FLAGS array(int)
    '''
    
    cdef int i                  = 0
    cdef np.ndarray SNR         =       np.zeros(256, dtype=float)
    cdef np.ndarray xmean       =       np.zeros(256, dtype=float)
    cdef np.ndarray arms        =       np.zeros(256, dtype=float)
    cdef float SNR2             =       avg
    cdef np.ndarray efficiency_x=       np.zeros(256, dtype=float)
    cdef list RFI_list          =       []
    cdef np.ndarray FLAGS       =       np.ones((256), dtype=int)   
    xmean                       =       np.mean(spec, axis=1)
    xrms                        =       rms(spec, xmean, fftsize)
    SNR                         =       xmean/xrms
    efficiency_x                =       SNR/SNR2
    for i in range(256):
        if(efficiency_x[i] > (np.mean(efficiency_x)+int(sig)*np.std(efficiency_x)) or  efficiency_x[i] < (np.mean(efficiency_x)-int(sig)*np.std(efficiency_x))):
            RFI_list.append(i)
            FLAGS[i]    =   0   
    
    
    return FLAGS

cdef Correlation_FORT(comf_X_1, comf_Y_1, comf1_X_1, comf1_Y_1, avg, le, file_name, file_name1, fftsize, external_num, c9, fil, chX1, chX2, chY1, chY2):#, delay):
    """
    Function to compute Cross-Correlations
    
    Args:
        comf_X_1
        comf_Y_1
        comf1_X_1
        comf1_Y_1
        avg
        le
        file_name
        file_name1
        fftsize
        external_num
        c9
        fil
        chX1
        chX2
        chY1
        chY2
    
    Returns:
        - **band**
        - **time1**
        - **c9**
        - **creal8**
        - **creal9**
    """

    time1          = np.zeros((2, le), dtype=complex)
    band           = np.zeros((2, fftsize), dtype=complex)
    timin      =       time.time()
    creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 = internal_loop2.external_loop(comf_X_1, comf_Y_1, comf1_X_1, comf1_Y_1, avg, le-1, 256-1)#, delay, 60000)
    print('Time for internal loop : ' +str(time.time()-timin))

    #####Generating New Folders, if not available#######################
    savenow     =   time.time()
    os.system('mkdir CORRELATION')
    st          =   'DELAY_Sample_len'+str(external_num)
    fg  =   str(file_name).split('/')
    fg1 =   fg[len(fg)-1]

    fh  =   str(file_name1).split('/')
    fh1 =   fh[len(fh)-1]
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX1)+str(chX2)+'_'+str(fg1)+'_'+str(fh1)+'.txt', creal8)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chY1)+str(chY2)+'_'+str(fg1)+'_'+str(fh1)+'.txt', creal9)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X1X1_'+str(fg1)+'_'+str(fh1)+'.txt', creal1)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_Y2Y2_'+str(fg1)+'_'+str(fh1)+'.txt', creal2)

    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX1)+str(chY1)+str(fg1)+'_'+str(fh1)+'.txt', pollkg1)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX2)+str(chY2)+str(fg1)+'_'+str(fh1)+'.txt', pollkg2)

    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX1)+str(chY2)+str(fg1)+'_'+str(fh1)+'.txt', creal5)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX2)+str(chY1)+str(fg1)+'_'+str(fh1)+'.txt', creal6)

    #####Checking Correlation and saving all stuff######################

    time1[0]    =   np.nanmean(creal9[20:230], axis =0)#/deno240
    band[0]     =   np.nanmean(creal9, axis =1)#/deno241

    time1[1]    =   np.nanmean(creal8[20:230], axis =0)#/deno80
    band[1]     =   np.nanmean(creal8, axis =1)#/deno81
    c9[int(external_num)]      =       np.nanmean(abs(band[0]))
    print('Mean  Y1Y2...'+str(np.nanmean(abs(np.nanmean(creal9[30:230], axis =1)))))
    print('Mean  X1X2...'+str(np.nanmean(abs(np.nanmean(creal8[30:230], axis =1)))))
    #np.save('CORRELATION/'+str(fil)+'/Correlation_X1X2_Y1Y2_Band'+str(fg1)+'_'+str(fh1)+'.txt',band)
    #np.save('CORRELATION/'+str(fil)+'/Correlation_X1X2_Y1Y2_Time'+str(fg1)+'_'+str(fh1)+'.txt',time1)

    print('Time for python function = '+ str(time.time()-timin))

    return band, time1, c9, creal8, creal9


class bcolors:
    """
    Class for colored cli Output
    """
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    
    pass

cdef PART_CORR(file_name, file_name1, sysgps,avg, RA, Dec, T1, T2, fftsiz):
    """
    Function for Correlation.

    Args:
        file_name
        file_name1
        sysgps,avg
        RA
        Dec
        T1
        T2
        fftsiz
    
    Returns:
        - **YYb**
        - **YYt**
        - **c9**
        - **creal8**
        - **creal9**
    """
    #Extracting date and time tags from file name and generating time jump matrix#
    series, sec, minu, hour, dat, month, year	        =	_header_Fring_cy.extract_time(file_name)
    series1, sec1, minu1, hour1, dat1, month1, year1	=	_header_Fring_cy.extract_time(file_name1)
    sec                                                 =       max(sec, sec1)

    #GPS_Synchronization starting#
    Memfactor = _header_Fring_cy.gps_sync(file_name, file_name1, sysgps)
    slope     = [['64453.125']]#get_slope(file_name, file_name1)

    #Generating Delay#

    print(sec, minu, hour, dat, month, year)
    print(sec1, minu1, hour1, dat1, month1, year1)
    lst1    =   _header_geometric_cy.selflst(int(sec), int(minu), int(hour), int(dat), int(month), int(year))
    sece    =   int(sec)    +   30
    mint    =   minu
    hourt   =   hour
    if(sece>59):
        mint    =   minu+1
        sece    =   sece-sec
        if(mint > 59):
            hourt   =   hour+1
            mint    =   0

    lst2    =   _header_geometric_cy.selflst(int(sece), int(mint), int(hourt), int(dat), int(month), int(year))
    series_int  =   float(series)*31+int(sec)
    delay                       =       _header_geometric_cy.Cal_time_onhold(float(sec), float(minu), float(hour), float(dat), float(month), float(year), float(RA),  float(Dec), float(avg), 1.0, 31, T1, T2)
    if(np.mean(delay) < 0):
        delay                       =       _header_geometric_cy.Cal_time_onhold(float(sec), float(minu), float(hour), float(dat), float(month), float(year), float(RA),  float(Dec), float(avg), 1.0, 40, T2, T1)

    timea    = time.time()


    comf_X, comf_Y, comf1_X, comf1_Y, LO    =  call_to_read(file_name, file_name1, avg, sys.argv[1], delay)
    print('Time taken to read files and decrypt...' +str(time.time()-timea))


    #Warning about swapping#
    #print(bcolors.WARNING+'Remember!! keep avg below 60000, or swapping will happen!!'+bcolors.ENDC)

    #GPS Compensation to find the correlation#

    if(sysgps == str(1)):
        comf_X    =   comf_X[int(round(Memfactor[0])):]#-int(sys.argv[5])/2:]#+int(round(float(sys.argv[5])))*1/2:] #This is to ensure the correlation peak is seen at the mid-point.
        comf1_X   =   comf1_X[int(round(Memfactor[1])):]

        comf_Y    =   comf_Y[int(round(Memfactor[0])):]#-int(sys.argv[5])/2:]#+int(round(float(sys.argv[5])))*1/2:]
        comf1_Y   =   comf1_Y[int(round(Memfactor[1])):]

    #Generatinf end point#

    le      = len(comf_X)/(int(int(sys.argv[-1])*2)*int(avg))
    le1     = len(comf1_X)/(int(int(sys.argv[-1])*2)*int(avg))
    le      =   min(le, le1)

    #Building file to save#
    tim         =   str('{0:02d}'.format((datetime.now().hour)))+str('{0:02d}'.format((datetime.now().minute)))+str('{0:02d}'.format((datetime.now().second)))
    dat         =   str('{0:04d}'.format((datetime.now().year)))+str('{0:02d}'.format((datetime.now().month)))+str('{0:02d}'.format((datetime.now().day)))
    corrfile    =   dat+'_'+tim+'.corr'

    #Range of samples to shift#
    #If input number if 10, then the samples are shifted from -10 to 10.

    c9          =       np.zeros((int(round(float(sys.argv[3])))*2+1))
    samp_range  =       np.linspace(0, int(round(float(sys.argv[3]))), int(round(float(sys.argv[3]))), dtype =int)#linspace(-int(round(float(sys.argv[5]))), int(round(float(sys.argv[5]))), int(round(float(sys.argv[5])))*2+1)
    count_c9    =   0
    #Always only comf will be shifted and not comf1#



    fill1    =   file_name.split('_')
    fill2    =   file_name1.split('_')


    chX1      =   'X'+str((fill1[-5])[-1])
    chX2      =   'X'+str((fill2[-5])[-1])
    chY1      =   'Y'+str((fill1[-5])[-1])
    chY2      =   'Y'+str((fill2[-5])[-1])







    #Writing Metadata#
    filll                         =   fill1[-4]+'_'+fill1[-3]+'_'+fill1[-3]#dat+'_'+tim
    os.system('mkdir CORRELATION/'+filll)
    point   =   open('CORRELATION/'+str(filll)+'/Metadata_'+str(file_name.split('/')[-1])+'_'+str(file_name1.split('/')[-1])+'.txt', 'w+')
    point.write('#SOURCE Obs_Frequency\tRA\tDec\tN-FFT\tIntegration_Time\tStart_Time_LST\tStop_Time_LST\n')
    point.write(str(filll)+'\t'+str(LO)+'\t'+str(RA)+'\n'+str(Dec)+'\t'+str(sys.argv[-1])+'\t'+str(avg)+'\t'+str(lst1)+'\t'+str(lst2))
    point.close()







    print('Time till loop..' +str(time.time()-timea))
    #comf_X1     =       comf_X.copy()
    #comf_Y1     =       comf_Y.copy()
    #comf_X  =   comf_X1[int(jjj)*1:]
    #comf_Y  =   comf_Y1[int(jjj)*1:]
    st                          =       str(filll)+"/DELAY_Sample_len0"
    os.system("mkdir CORRELATION/"+str(filll)+"/DELAY_Sample_len0")
    YYb, YYt, c9, creal8, creal9        =       Correlation_FORT(comf_X, comf_Y, comf1_X, comf1_Y, avg, le,  file_name, file_name1, int(sys.argv[-1]), 0, c9, filll, chX1, chX2, chY1, chY2)#, delay)
    count_c9    =   count_c9+1

    return YYb, YYt, c9, creal8, creal9

cdef plot_all(file_name, file_name1):
    """
    Function for plotting spec XX and YY.

    Args:
        file_name (str): File 1 Name
        file_name1 (str): File 2 Name
    
    Returns:
        - **spXX** array spec XX
        - **spYY** array spec YY
    """

    fill1    =   file_name.split('_')
    fill2    =   file_name1.split('_')

    cdef str filll                         =   fill1[-4]+'_'+fill1[-3]+'_'+fill1[-3]#dat+'_'+tim
    

    cdef str chX1      =   'X'+str((fill1[-5])[-1])
    cdef str chX2      =   'X'+str((fill2[-5])[-1])
    cdef str chY1      =   'Y'+str((fill1[-5])[-1])
    cdef str chY2      =   'Y'+str((fill2[-5])[-1])

    ltXX      =   glob.glob('CORRELATION/'+filll+'/Correlation_'+str(chX1)+str(chX2)+'*')
    ltYY      =   glob.glob('CORRELATION/'+filll+'/Correlation_'+str(chY1)+str(chY2)+'*')
    spXX    =   'spec =  np.hstack(('
    spYY    =   'specY = np.hstack((' 
    cdef list sp1     =   []
    cdef list sp2     =   []
    cdef int i        =   0


    for i in range(len(ltXX)):
        sp1.append(load(ltXX[i]))
        sp2.append(load(ltYY[i]))
        spXX  =   spXX+'sp1['+str(i)+'], '
        spYY  =   spYY+'sp2['+str(i)+'], '
    spXX  =   spXX+'))'
    spYY  =   spYY+'))'
    exec(spXX)
    exec(spYY)

    
    savetxt('CORRELATION/Correlation_filter_'+str(chX1)+str(chX2), spec)
    savetxt('CORRELATION/Correlation_filter_'+str(chY1)+str(chY2), spec1)

    return spXX, spYY

cpdef CORR(sysgps, fil1, fil2, avg, RA, Dec, T1, T2, fftsiz):
    """
    Args:
        sysgps
        fil1
        fil2
        avg
        RA
        Dec
        T1
        T2
        fftsiz

    Returns:
        0

    Raises:
        RuntimeError
    """

    cdef unsigned int j = 0

    if(len(fil1)!=len(fil2)):
        raise RuntimeError ('No. of files are not equal!!')


    for j in range(len(fil1)):
        file_name                       =   fil1[j]
        file_name1                      =   fil2[j]
        YYb, YYt, c9, creal8, creal9    =   PART_CORR(file_name, file_name1, sysgps, avg, RA, Dec, T1, T2, fftsiz)
    plot_all(file_name, file_name1)
    return 0


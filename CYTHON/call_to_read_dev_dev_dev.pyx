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
import sys
import _header_gps_cy
from cython.view cimport array as cvarray
import cython
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
    #cdef list line2 =   []
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
    slope_1  =   np.mean(slope1)
    slope_2  =   np.mean(slope2)
    inter_1  =   np.mean(inter1)
    inter_2  =   np.mean(inter2)
    for i in range(len(f1)):
        tf1    =   open(str("SAMPLING_INFO/Info_on_straight_line"+str(f1[i].split('/')[-1])), 'rw')
        tf2    =   open(str("SAMPLING_INFO/Info_on_straight_line"+str(f2[i].split('/')[-1])), 'rw')
        temp1  =   tf1.readlines()
        temp2  =   tf2.readlines()
        if(inter >0):
            temp1[2] =   str(slope)+'*x+'+str(inter)+'\n' #Replaced slope_1 and inter_1 with slope and inter..Remember!!
            temp2[2] =   str(slope)+'*x+'+str(inter)+'\n'
        else:
            temp1[2] =   str(slope)+'*x-'+str(inter)
            temp2[2] =   str(slope)+'*x-'+str(inter)
        tf1.close()
        tf2.close()
        tf1    =   open(str("SAMPLING_INFO/Info_on_straight_line"+str(f1[i].split('/')[-1])), 'w+')
        tf2    =   open(str("SAMPLING_INFO/Info_on_straight_line"+str(f2[i].split('/')[-1])), 'w+')

        tf1.writelines(temp1)
        tf2.writelines(temp2)
        #Not Applicable..#Replaced here from temp2 to temp1 for both files to have same synchronization equation#
        #tf2.writelines(temp1)
        tf1.close()
        tf2.close()
    return 0

cpdef tuple decrypy_get_gpssync(file_name, file_name1, avg, s_mbr='000.mbr'):
    '''
        Getting GPS sync data..

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
    

    for i in range(10, len(comf)-1):
        if(comf['GPS'][i+1]- comf['GPS'][i]  ==  1 and comf['Packet'][i+1] - comf['Packet'][i] == 1):
            accblipp1.append(float(comf['Packet'][i+1]))
            accblipg1.append(float(comf['GPS'][i+1]))
    #Important message to anyone opening this file.. 
    #The condition of the diff of GPS and the Packet was added as any packet loss during the transition of GPS would
    #skew the resulting straight line equation towards more error..hence the two condition of both Packet and GPS..

    for i in range(10, len(comf1)-1):
        if(comf1['GPS'][i+1] - comf1['GPS'][i] ==1 and comf1['Packet'][i+1] - comf1['Packet'][i] == 1):
            accblipp2.append(float(comf1['Packet'][i+1]))
            accblipg2.append(float(comf1['GPS'][i+1]))
    np.savetxt('Packet_GPS', [accblipp1, accblipg1])
    slp1    =    linregress(accblipg1, accblipp1)
    slp2    =    linregress(accblipg2, accblipp2)
    print(slp1)
    print(slp2)
    slp_1   =    linear_fit(accblipg1, accblipp1)
    slp_2   =    linear_fit(accblipg2, accblipp2)
    
    slp1_eqn=       str(slp1.slope)+'*x+'+str(slp1.intercept) 
    slp2_eqn=       str(slp2.slope)+'*x+'+str(slp2.intercept)
    x   =   max(comf['GPS'][0], comf1['GPS'][0])+1
    

    fillist         =               np.loadtxt('DEPEND/files_naming.txt', dtype=str)
    series_code     =               int(np.where(fillist==series)[-1][0])#fillist.index(series)+1
    print('\n\n\n\n\n') 
    print(series_code)

    
    
    #fil_000         =               file_name[:-7]+'000.mbr'
    #fil1_000        =               file_name1[:-7]+'000.mbr'
    


    fileeq1         =               open(str("SAMPLING_INFO/Info_on_straight_line"+str(file_name.split('/')[-1])), 'w+')
    fileeq2         =               open(str("SAMPLING_INFO/Info_on_straight_line"+str(file_name1.split('/')[-1])), 'w+')
    #RFI_mask1       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil_000.split('/')[-1])))#+fillist[0])
    #RFI_mask2       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil_000.split('/')[-1])))
    #RFI_mask3       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil1_000.split('/')[-1])))#+fillist[0])
    #RFI_mask4       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil1_000.split('/')[-1])))

    #fileeq1.write('First GPS Value == '+str(accblipg1[0])+'\n'+str(accblipp1[0])+'\n')
    #fileeq2.write('First GPS Value == '+str(accblipg2[0])+'\n'+str(accblipp2[0])+'\n')

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
        fileeq1.write(str(slp2_eqn)+'\n')
        fileeq2.write(str(slp2_eqn)+'\n')
        
        file_to_write       =   file_name.split('_')
        file_to_write_fil   =   file_name.split('_')
        file_to_write[-1]   =   fillist[int(series_code)+1]
        file_to_write_fil[-1]=  fillist[series_code]
        file_to_write       =   '_'.join(file_to_write)
        print(file_to_write)
        file_to_write_fil   =   '_'.join(file_to_write_fil)

        filrem          =               open('SAMPLING_INFO/'+str(file_to_write_fil.split('/')[-1])+'_rem.data', 'w+')
        filrem.write(str(rem_len)+'\n')#str(eval(slp_2[-1])-comf1['Packet'][0]-eval(slp_1[-1])+comf['Packet'][0])+'\n')
        filrem.write(str(file_to_write)+'\n')

        print(fillist[int(series_code)+1])
        fileeq2.write('N,0\n')
        fileeq1.write('Y,'+str(file_to_write)+','+str(rem_len)+'\n')


        #try:
        #    fil_to_rem          =   (file_to_write.split('/')[-1])[0:4]
        #    file_to_write       =   (file_to_write.split('/')[-1])[-7:-4]
        #except:
        #    fil_to_rem          =   file_to_write[0:4]
        #    file_to_write       =   file_to_write[-7:-4]
        #filrem.write(str(fil_to_rem))
        #print(file_to_write)
        #print(fil_to_rem)
    else:
        fileeq1.write(str(slp1_eqn)+'\n')
        fileeq2.write(str(slp1_eqn)+'\n')


        file_to_write       =   file_name1.split('_')
        file_to_write_fil   =   file_name1.split('_')
        file_to_write[-1]   =   fillist[int(series_code)+1]
        file_to_write_fil[-1]=  fillist[series_code]
        file_to_write       =   '_'.join(file_to_write)
        print('=============')
        print(file_to_write)
        print('=============')
        file_to_write_fil   =   '_'.join(file_to_write_fil)


        filrem          =               open('SAMPLING_INFO/'+str(file_to_write_fil.split('/')[-1])+'_rem.data', 'w+')
        filrem.write(str(rem_len)+'\n')#str(eval(slp_1[-1])-comf['Packet'][0]-eval(slp_2[-1])+comf1['Packet'][0])+'\n')
        #filrem.write(str(min(eval(slp_1[-1]), eval(slp_2[-1]))-comf['Packet'][0])+'\n')
        #filrem(str())
        filrem.write(str(file_to_write)+'\n')

        #try:
        #    fil_to_rem          =   (file_to_write.split('/')[-1])[0:4]
        #    file_to_write       =   (file_to_write.split('/')[-1])[-7:-4]
        #except:
        #    fil_to_rem          =   file_to_write[0:4]
        #    file_to_write       =   file_to_write[-7:-4]
        #print(file_to_write)
        #print(fil_to_rem)
        #filrem.write(str(fil_to_rem))

        fileeq1.write('N,0\n')
        fileeq2.write('Y,'+str(file_to_write)+','+str(rem_len)+'\n')



    fileeq1.write(str(comf['GPS'][0])+','+str(comf['GPS'][-1]))
    fileeq2.write(str(comf1['GPS'][0])+',' + str(comf1['GPS'][-1]))

    num =   comf['Packet'][int(np.where(comf['GPS']==comf['GPS'][-1])[0][0])]
    num1=   comf1['Packet'][int(np.where(comf1['GPS']==comf1['GPS'][-1])[0][0])]

    
    fileeq1.write('\n'+str(comf['Packet'][-1])+'-'+str(num)+' = '+str(comf['Packet'][-1]-num)+','+ str(comf['Packet'][0]))
    fileeq2.write('\n'+str(comf1['Packet'][-1])+'-'+str(num)+' = '+str(comf1['Packet'][-1]-num1)+','+ str(comf1['Packet'][0]))

    print(str(slp1.slope)+'*x'+str(slp1.intercept))
    
    #Generating RFI mask#
    if(series == s_mbr):
        comf_X, comf_Y, comf1_X, comf1_Y, LO, Memfactor_not_used_here = call_to_read(file_name, file_name1, 60, '1', np.zeros((30)))
        le      = len(comf_X)/(int(int(256)*2)*int(avg))
        le1     = len(comf1_X)/(int(int(256)*2)*int(avg))
        le      =   min(le, le1)
        creal1, creal2, creal3, creal4, creal8, creal9= internal_loop2.external_loop(comf_X, comf_Y, comf1_X, comf1_Y, avg, le, 255)
        RFI1    =   RFI_Reject(creal1, 1, 256, 60)
        RFI2    =   RFI_Reject(creal2, 1, 256, 60)
        RFI3    =   RFI_Reject(creal3, 1, 256, 60)
        RFI4    =   RFI_Reject(creal4, 1, 256, 60)
      
       #Writing RFI masks#
        np.savetxt(str("SAMPLING_INFO/RFI_mask_X_"+str(file_name.split('/')[-1])), RFI1)
        np.savetxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(file_name.split('/')[-1])), RFI2)
        np.savetxt(str("SAMPLING_INFO/RFI_mask_X_"+str(file_name1.split('/')[-1])), RFI3)
        np.savetxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(file_name1.split('/')[-1])), RFI4)

    
    return slp1, slp2, slp_1, slp_2#, RFI_mask1, RFI_mask2, RFI_mask3, RFI_mask4

cpdef get_memfactor_est(str file_name, str file_name1):
    '''
         Memfactor estimation using packet numbering...
    '''
    dt  =       np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    comf        =       np.memmap(file_name, dtype=dt, mode='c')
    comf1       =       np.memmap(file_name1, dtype=dt, mode='c')
    series      =       file_name[-7:]
    if(series != '000.mbr'):
        gp1_0   =       comf['GPS'][0]
        gp2_0   =       comf1['GPS'][0]
        gp_mx   =       max(gp1_0, gp2_0)+1
    else:
        gp1_0    =       comf['GPS'][100]
        gp2_0    =       comf1['GPS'][100]
        gp_mx    =       max(gp1_0, gp2_0)+1

    print(gp_mx)
    tag_gp1      =       np.where(comf['GPS'] == gp_mx)[0]
    tag_gp2      =       np.where(comf1['GPS'] == gp_mx)[0]
    print(tag_gp1, tag_gp2)
    return tag_gp1, tag_gp2


cpdef read_RFI(file_name, file_name1):


    fil_000         =               file_name[:-7]+'000.mbr'
    fil1_000        =               file_name1[:-7]+'000.mbr'


    RFI_mask1       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil_000.split('/')[-1])))#+fillist[0])
    RFI_mask2       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil_000.split('/')[-1])))
    RFI_mask3       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil1_000.split('/')[-1])))#+fillist[0])
    RFI_mask4       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil1_000.split('/')[-1])))
    
    aff1            =               np.asarray((np.unique(RFI_mask1, return_counts=True)))
    aff2            =               np.asarray((np.unique(RFI_mask2, return_counts=True)))
    aff3            =               np.asarray((np.unique(RFI_mask3, return_counts=True)))
    aff4            =               np.asarray((np.unique(RFI_mask4, return_counts=True)))
    ch01            =               1
    ch02            =               1
    ch03            =               1
    ch04            =               1

    if(aff1[:,0][-1] > 100):
        ch01        =               -1

    if(aff2[:,0][-1] > 100):
        ch02        =               -1

    if(aff3[:,0][-1] > 100):
        ch03        =               -1

    if(aff4[:,0][-1] > 100):
        ch04        =               -1

    return RFI_mask1, RFI_mask2, RFI_mask3, RFI_mask4, [ch01, ch02, ch02, ch04]

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef compensate_pack_loss(comf, comf1, tempcomf, tempcomf1, lin1, lin2, Mem1, Mem2):
     
    search_time       =   time.time()
    cdef np.ndarray pack_loss         =   comf['Packet'][1:] - comf['Packet'][0:len(comf['Packet'])-1]  
    cdef np.ndarray pack_loss1        =   comf1['Packet'][1:] - comf1['Packet'][0:len(comf1['Packet'])-1]
    cdef np.ndarray pack_loss_t1      =   np.where(pack_loss != 1)[0]
    cdef np.ndarray pack_loss_t2      =   np.where(pack_loss1 != 1)[0]
    cdef np.ndarray pack_loss_val1    =   pack_loss[pack_loss_t1]
    cdef np.ndarray pack_loss_val2    =   pack_loss1[pack_loss_t2]



    time_marker       =   time.time()
    marker_str1       =     'np.hstack(('
    marker_str2       =     'np.hstack(('
    for i in range(len(pack_loss_t1)):
        marker_str1      =           marker_str1+'np.linspace(lin1[pack_loss_t1['+str(i)+']]+1, (pack_loss_val1['+str(i)+']+lin1[pack_loss_t1['+str(i)+']]-1), pack_loss_val1['+str(i)+']-1, dtype = int) -Mem1 , '
    for j in range(len(pack_loss_t2)):
        marker_str2      =           marker_str2+'np.linspace(lin2[pack_loss_t2['+str(j)+']]+1, (pack_loss_val2['+str(j)+']+lin2[pack_loss_t2['+str(j)+']]-1), pack_loss_val2['+str(j)+']-1, dtype = int) -Mem2 , '
    marker_str1          =           marker_str1+'))'
    marker_str2          =           marker_str2+'))' 
    cdef np.ndarray marker1=np.array([])
    cdef np.ndarray marker2=np.array([])
    if(len(pack_loss_t1)>0):
        print('Packet loss in the file 1 was seen..')
        marker1 =           eval(marker_str1)
    if(len(pack_loss_t2)>0):
        print('Packet loss in the file 2  was seen..')
        marker2 =           eval(marker_str2)

    print('Marker values are \n\n\n')
    #print(marker1)
    #print(marker1)
    #Removing any negative index numbers, arising due to the subtraction of Mem2, i.e packet loss before the GPS+pack_loss jump
    marker1	=marker1[len(marker1[np.where(marker1<0)]):]	
    marker2	=marker2[len(marker2[np.where(marker2<0)]):]
    #print(marker_str1)
    #print(marker_str2)
    cdef np.ndarray marker              =           np.hstack((marker1, marker2))
    #np.save('marker1', marker1)
    #np.save('marker2', marker2)


    print('Marker Time..'+str(time.time() - time_marker))
    print(Mem1, Mem2)


    time_delete       =   time.time()
    #np.save('tempcomf', tempcomf)
    #np.save('tempcomf1', tempcomf1)
    #tempcomf		 =   compensate_packet_loss_delete_function(tempcomf[Mem1:], marker)
    #tempcomf1		 =   compensate_packet_loss_delete_function(tempcomf1[Mem2:], marker)



    #tempcomf_1           =   np.delete(tempcomf[Mem1:],  marker , axis = 0)
    #tempcomf_2           =   np.delete(tempcomf1[Mem2:], marker , axis = 0)
    tempcomf           	=   np.delete(tempcomf,  marker , axis = 0)
    tempcomf1          =   np.delete(tempcomf1,  marker , axis = 0)

    #tempcomf            =   np.vstack((tempcomf[:Mem1], tempcomf_1))
    ##tempcomf1           =   np.vstack((tempcomf1[:Mem2], tempcomf_2))

    #print('Search time..'+str(time.time() - search_time))
    #print('Delete Time..'+str(time.time() - time_delete))





    tempcomf    = tempcomf.ravel()
    tempcomf1   = tempcomf1.ravel()
    #print('Compensated for the packet loss, the final length of tempcomf and tempcomf1 are..'+str(len(tempcomf))+' and '+str(len(tempcomf1)))
    return tempcomf, tempcomf1

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef compensate_packet_loss_delete_function(temp_comf, marker):
    
    cdef np.ndarray tempcomf_1           =   np.delete(temp_comf,  marker , axis = 0)
    #tempcomf_2           =   np.delete(tempcomf1[Mem2:], marker , axis = 0)
    cdef np.ndarray tempcomf             =   np.vstack((temp_comf, tempcomf_1)).ravel()
    #tempcomf1           =   np.vstack((tempcomf1[:Mem2], tempcomf_2))
    return tempcomf

cpdef compensate_pack_loss_rem(comf_rem, tempcomf_rem, lin1, val):
    
    search_time       				  =   time.time()
    cdef np.ndarray pack_loss         =   comf_rem['Packet'][1:val] - comf_rem['Packet'][0:val-1]  
    cdef np.ndarray pack_loss_t1      =   np.where(pack_loss != 1)[0]
    cdef np.ndarray pack_loss_val1    =   pack_loss[pack_loss_t1]
    print(pack_loss_t1)
    print(pack_loss_val1)


    time_marker       =   time.time()
    marker_str1       =     'np.hstack(('
    if(len(pack_loss_t1) == 0):
        print('No data packets loss seen in the next file..hence skipping the packet loss compensation for the non intergral file..')
        return tempcomf_rem.ravel()
    for i in range(len(pack_loss_t1)):
        marker_str1      =           marker_str1+'np.linspace(lin1[pack_loss_t1['+str(i)+']]+1, (pack_loss_val1['+str(i)+']+lin1[pack_loss_t1['+str(i)+']]-1), pack_loss_val1['+str(i)+']-1, dtype = int) , '
    marker_str1          =           marker_str1+'))'
    print(marker_str1)
    cdef np.ndarray marker             =           eval(marker_str1)



    tempcomf_rem          =   np.delete(tempcomf_rem,  marker , axis = 0)
    tempcomf_rem    		= 	tempcomf_rem.ravel()
    
    return tempcomf_rem

cpdef tuple decrypy_file_new_SWAN_onhold_without_wierd_comp_for_testing_ONLY(file_name, file_name1, ch, Memfactor=np.array([0.0, 0.0])):


    '''
                ch should be the channel number of first file..

		Takes the read file and sorts the X and Y polarizartion in the file into

		comf_X, comf1_X, comf_Y, comf1_Y.

    '''
    #If available get the earlier data set#
    cdef str series      = file_name[-7:-4]
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
     
    cdef int LO1, LO2

    cdef double timea    = time.time()
    cdef str fil_tag     = file_name.split('_')[-4]
    cdef int    i        = 0
    print('Caching memory, for faster read..')    
    #os.system('./comp '+str(file_name)+' '+str(file_name1))
    comf     = np.memmap(file_name,  dtype = dt, mode = 'c')
    comf1    = np.memmap(file_name1, dtype = dt, mode = 'c')
    
    
    
    #Finding the number of packets to jump to get to the synchronization from the linear equation result#
    tim_find = time.time()
    cdef int Mem1_index =   0
    cdef int Mem2_index =   0
    cdef int Mem1_fact      =   0
    cdef int Mem2_fact      =   0
    cdef int Mem1           =   0
    cdef int Mem2           =   0
    for i in range(len(comf)):
        Mem1_fact           =   int(float(round(Memfactor[0]/512.0)))
        Mem2_fact           =   int(float(round(Memfactor[1]/512.0)))

        if(comf['Packet'][i] == Mem1_fact and Mem1_index==0):
            Mem1             =  i
            Mem1_index       =  1
        if(comf1['Packet'][i] == Mem2_fact and Mem2_index==0):
            Mem2             =  i
            Mem2_index       =  1
        if(Mem1_index & Mem2_index  ==  1):
            break;
    print('Memfactor 1.....'+str(Mem1))
    Memfactor[0]    =   Mem1*512 + (Memfactor[0]%1)*512
    Memfactor[1]    =   Mem2*512 + (Memfactor[1]%1)*512
    
    
    cdef double ploss1          =       comf['Packet'][-1]  -  comf['Packet'][0] -len(comf['Packet'])
    cdef double ploss2          =       comf1['Packet'][-1] -  comf1['Packet'][0] - len(comf1['Packet'])
    cdef np.ndarray Memfactor1  =       np.zeros((2), dtype=float)
    Memfactor1                  =       gps_alingn_files(ploss1, ploss2, Memfactor[0], Memfactor[1], file_name, file_name1, series)
    #Changing original Memfactor with the new Memfactor
    print('Before Memfactor change..Memfactor[0], Memfactor[1]..'+str(Memfactor[0])+','+str(Memfactor[1]))
    Memfactor                    =  Memfactor1
    print('After Memfactor change..Memfactor[0], Memfactor[1]..'+str(Memfactor[0])+','+str(Memfactor[1]))
    
    print('Time required to read..'+str(time.time()-timea))
    print('Time required to find the Memfactor..'+str(time.time() - tim_find))

    
    
    print(len(comf), len(comf1))
    LO1                  =   comf['LO'][100]
    LO2                  =   comf1['LO'][100]

    if(LO1!=LO2):
        print('LO1 = ' +str(LO1)+'\n')
        print('LO2 = ' +str(LO2)+'\n')
        #raise RuntimeError ('Both LOs are different!')


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
    cdef np.ndarray  lin1       =   comf['Packet'][10:]  -len1_fil1
    cdef np.ndarray  lin2       =   comf1['Packet'][10:] -len1_fil2
    tim_allot                   =   time.time()
    cdef np.ndarray tempcomf    #=   np.zeros((templen, 1024), dtype = np.int8) 
    cdef np.ndarray tempcomf1   #=   np.zeros((templen1, 1024), dtype = np.int8)
    print('Time for array array allotment...'+str(time.time()-tim_allot))
    
    tempcomf          =   np.memmap.copy(comf['data'][10:])

    print('Time required to decrypt one file..'+str(time.time()-timea))
    tempcomf1         =   np.memmap.copy(comf1['data'][10:])
    

    tempcomf        =       tempcomf.ravel()
    tempcomf1       =       tempcomf1.ravel()

    tempcomf_X                          =   np.array(tempcomf[1::2],  order = 'F')
    tempcomf_Y                          =   np.array(tempcomf[0::2],  order = 'F')
    tempcomf1_X                         =   np.array(tempcomf1[1::2], order = 'F')
    tempcomf1_Y                         =   np.array(tempcomf1[0::2], order = 'F')
    
    print('Time required to decrypt both files..'+str(time.time()-timea))
    print(len(tempcomf_X), len(tempcomf1_X))
    return tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, LO1, Memfactor1

cpdef get_last_gps_timing(comf, comf1, file_name, file_name1, series, Memfactor):
    '''
         Module to calculate last GPS value in the files
         This will help constrain the simulated delay error.
         INPUT:
                 Memmaped array file1 (comf),  Memmaped array file1 (comf1), file path (including name) of first file,
                     file path (including name) of second file, series code, Memfactor array 
         OUTPUT:
                 Total time of file1, Total time of file2 (in sec)  
    '''

    cdef str file_name_1, file_name1_1 
    
    try:
        file_name_1 =       file_name.split('/')[-1]
        file_name1_1=       file_name1.split('/')[-1]
    except:
        file_name_1 =       file_name
        file_name1_1=       file_name1
    file_name_1		=		file_name_1[:-7]+'000.mbr'
    file_name1_1	=		file_name1_1[:-7]+'000.mbr'
    

    #Last GPS value recorded in the packet
    cdef long long int last_gps1	=	0
    cdef long long int last_gps2        =       0
    cdef long long int penul_gps1	=	0
    cdef long long int penul_gps2	=	0
    cdef long long int last_pack1       =       0
    cdef long long int last_pack2       =       0

 
    cdef double	       comp1		=	0
    cdef double        comp2            =       0
   
    cdef long long int first_gps1	=	0
    cdef long long int first_gps2	=	0



    #Getting last GPS Value associated packet number 
    last_gps1				=	comf['GPS'][len(comf)-1]
    last_gps2                           =       comf1['GPS'][len(comf1)-1]

    last_pack1				=	comf['Packet'][len(comf)-1]
    last_pack2                          =       comf1['Packet'][len(comf1)-1]

    #Getting penultimate GPS Value associated packetnumper
    penul_gps1				=	comf['Packet'][np.where(comf['GPS'] == last_gps1-1)[0][0]]
    penul_gps2                          =       comf1['Packet'][np.where(comf1['GPS'] == last_gps2-1)[0][0]]
    
    #Getting difference 
    
    cdef int gps_check			=	0

    delta_gps1				=	last_pack1	-	penul_gps1
    delta_gps2				=	last_pack2	-	penul_gps2

    #Getting total packet number
    cdef long long int pack1		
    cdef long long int pack2            

    if(series=='000'):
        first_gps1			=	comf['GPS'][100]
        first_gps2			=	comf1['GPS'][100]
        pack1            		=       comf['Packet'][np.where(comf['Packet'] == penul_gps1)[0][0]] - comf['Packet'][100]    +delta_gps1
        pack2                           =       comf1['Packet'][np.where(comf1['Packet'] == penul_gps2)[0][0]] -comf1['Packet'][100]   +delta_gps2
    else:
        first_gps1                      =       comf['GPS'][0]
        first_gps2                      =       comf1['GPS'][0]
        pack1                           =       comf['Packet'][np.where(comf['Packet'] == penul_gps1)[0][0]] - comf['Packet'][0] +delta_gps1
        pack2                           =       comf1['Packet'][np.where(comf1['Packet'] == penul_gps2)[0][0]]- comf1['Packet'][0]+delta_gps2


    #IMP: Calculated slope values are only taken from the 000 series..for consistancy and ease of operation. 
    cdef float slope1			=	gps_slope(file_name_1)
    cdef float slope2			=	gps_slope(file_name1_1)
    cdef float delta_time1		=	(pack1-Memfactor[0]/512.0)/slope1	
    cdef float delta_time2             	=       (pack2-Memfactor[1]/512.0)/slope2
    
    #Saving these values if not already saved..
        

 
    return np.array([[first_gps1+(Memfactor[0]/512)/slope1, delta_time1], [first_gps2+(Memfactor[1]/512)/slope2, delta_time2]])

cpdef get_gps_info_000_file(comf, comf1, Memfactor):
    cdef long long int Mem1_index =   0
    cdef long long int Mem2_index =   0
    cdef long long int Mem1_fact      =   0
    cdef long long int Mem2_fact      =   0
    cdef long long int Mem1           =   0
    cdef long long int Mem2           =   0
    
    print('Memfactor in get_gps_info_000_file..'+str(Memfactor[0])+','+str(Memfactor[1]))
    Mem1_fact           =   int(float(Memfactor[0]/512.0))
    Mem2_fact           =   int(float(Memfactor[1]/512.0))
    print('Mem1_fact, Mem2_fact'+str(Mem1_fact)+','+str(Mem2_fact))

    for i in range(len(comf)):
        if(comf['Packet'][i] == Mem1_fact and Mem1_index==0):
            Mem1             =  i
            Mem1_index       =  1
        if(comf1['Packet'][i] == Mem2_fact and Mem2_index==0):
            Mem2             =  i
            Mem2_index       =  1
        if(Mem1_index & Mem2_index  ==  1):
            break;
    cdef long long int pre_ploss	=	comf['Packet'][Mem1] - comf['Packet'][0] - Mem1#len(comf['Packet'][:Mem1])+1
    cdef long long int pre_ploss1	=	comf1['Packet'][Mem2] - comf1['Packet'][0] - Mem2#len(comf1['Packet'][:Mem2])+1
    print('Memfactor 1.....'+str(Mem1))
    print('Memfactor 2.....'+str(Mem2))
   
    print('pre_ploss, pre_ploss1 in get_gps_info_000_file is..'+str(pre_ploss)+','+str(pre_ploss1)) 
    Mem1	    =	Mem1 + pre_ploss #+1
    Mem2	    =   Mem2 + pre_ploss1#+1   

    Memfactor[0]    =   Mem1*512 + (Memfactor[0]%512)#*512
    Memfactor[1]    =   Mem2*512 + (Memfactor[1]%512)#*512

    print('Memfactor in get_gps_info_000_file..'+str(Memfactor[0])+','+str(Memfactor[1]))

    return Mem1, Mem2, np.array(Memfactor)

cpdef tuple decrypy_file_new_SWAN_onhold(file_name, file_name1, ch, Memfactor=[0, 0]):


    '''
         Module to generate the decrypted version of the files, with separate set of X and Y Pol array.
         Input:
              (file_name, file_name1, ch, Memfactor=[0, 0])
              file_name	=	File name of the first file.
              file_name1=	File name of the second file.
              ch 	=	Avegrage
              Memfactor	=	Memory factor intented for the correction, only required for the 000 file and not required for the
              			non-000 files.
        Output:
              1. File One X-Pol 1D np array, 2. File Two Y-Pol 1D np array, 3. File Two X-Pol 1D np array, 4. File Two Y-Pol 1D np array, 5. Local Oscillator value derived
              from the header (int), 6. Memfactor 1D np array (Jump factor to use).
              
    '''
    #If available get the earlier data set#
    cdef str series      = file_name[-7:-4]
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
     
    cdef int LO1, LO2


    cdef np.ndarray Memfactor1  =       np.zeros((2), dtype=float)
    cdef double timea    = time.time()
    cdef str fil_tag     = file_name.split('_')[-4]
    cdef int    i        = 0
    print('Caching memory, for faster read..')    
    #os.system('./comp '+str(file_name)+' '+str(file_name1))
    #comf     = np.memmap(file_name,  dtype = dt, mode = 'c')
    #comf1    = np.memmap(file_name1, dtype = dt, mode = 'c')
    
    #Temporary additions...

    #comf11	=	np.memmap(file_next, dtype = dt, mode = 'c')
    #comf1	=	np.concatenate((comf1, comf11[0:600000]))
    time_read		=		time.time() 
    cdef np.ndarray file_time
    #Calling reading_function_spinoff
    tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, Memfactor, Mem1, Mem2, LO1, LO2, file_time	=	read_spinoff(series, Memfactor, file_name, file_name1, dt)  
    print('Time for readspinoff...'+str(time.time()-time_read))

    Memfactor1		 =   Memfactor

    #Now reading next file, to make up for the packet loss and GPS shift#
    cdef int le      			= 	(len(tempcomf_X)-int(round(Memfactor[0])))/(int(512)*int(ch)) #- int(round(Memfactor[0]))/512.0
    cdef int le1     			= 	(len(tempcomf1_X)-int(round(Memfactor[1])))/(int(512)*int(ch))#- int(round(Memfactor[1]))/512.0
    file_name_1, file_name1_1           =       break_file_name(file_name, file_name1)   


    
    time_jumpX, time_jumpY, time_flag	=	 get_fparam(file_name_1, file_name1_1) 
    le      				=   	 min(le, le1)
    RFI 				=	 _header_gps_cy.read_RFI(file_name, file_name1) 
    '''
        dt_Memfactor composition
            Pol-X   Pol-Y
        0
        1
        
        
        dt_Memfactor_flag composition
            Pol-X   Pol-Y
        0
        1



    '''

    dt_Memfactor			=	np.dtype([('MemfactorX', np.float64), ('MemfactorY', np.float64)])#, [('Memfacto1X', '>u4'), ('Memfactor1Y', '>u4')]])
    new_Memfactor			=	np.zeros((2, 2), dtype=dt_Memfactor)
    dt_Memfactor_flag                   =       np.dtype([('MemfactorX_flag', '>i1'), ('MemfactorY_flag', '>i')])
    new_Memfactor_flag			=	np.zeros((2), dtype=dt_Memfactor_flag)

    print('Memfactor last..')
    print(Memfactor)
    if(Memfactor[0]!=0): 
        new_Memfactor['MemfactorX'][0][0]	     =	    	Memfactor[0]+time_jumpX
        new_Memfactor['MemfactorX'][0][1]	     =		0
        new_Memfactor['MemfactorY'][0][0]	     =		Memfactor[0]+time_jumpY
        new_Memfactor['MemfactorY'][0][1]	     =		0
        new_Memfactor['MemfactorX'][1][0]            =          (Memfactor[0]%512+time_jumpX)
        if(new_Memfactor['MemfactorX'][1][0]<0):
            print("In new_Memfactor['MemfactorX'][1][0]<0 condition")
            new_Memfactor['MemfactorX'][1][0]        =          (Memfactor[0]%512+time_jumpX)%512
        
        new_Memfactor['MemfactorX'][1][1]            =          0
        
        new_Memfactor['MemfactorY'][1][0]            =          (Memfactor[0]%512+time_jumpY)
        if(new_Memfactor['MemfactorY'][1][0]<0):
            print("In new_Memfactor['MemfactorY'][1][0]<0 condition")
            new_Memfactor['MemfactorY'][1][0]        =          (Memfactor[0]%512+time_jumpY)%512
        
        new_Memfactor['MemfactorY'][1][1]            =          0


        if( int((Memfactor[0]+time_jumpY)/512) < int(Memfactor[0]/512)):
            tempcomf1_Y                       		 =          tempcomf1_Y[int(512-(Memfactor[0]%512+time_jumpY))%512:]
            new_Memfactor['MemfactorY'][1][0]            =          0#(Memfactor[0]%512+time_jumpY)%512
            print('\n\n\n\n')
        if( int((Memfactor[0]+time_jumpX)/512) < int(Memfactor[0]/512)):
            tempcomf1_X                       		 =          tempcomf1_X[int(512-(Memfactor[0]%512+time_jumpX))%512:]
            new_Memfactor['MemfactorX'][1][0]            =          0#(Memfactor[0]%512+time_jumpX)%512
            print('\n\n\n\n')
        print(time_jumpY, time_jumpX)
        print('new_Memfactor at the end is..')
        print(new_Memfactor)
    else:



        new_Memfactor['MemfactorX'][0][0]  =       0#Memfactor[0] - time_jumpX
        new_Memfactor['MemfactorY'][0][0]  =       0#Memfactor[0] - time_jumpY

        new_Memfactor['MemfactorX'][0][1]  =       Memfactor[1] - time_jumpX#X
        new_Memfactor['MemfactorY'][0][1]  =       Memfactor[1] - time_jumpY#Y



        new_Memfactor['MemfactorX'][1][0]  =       0#Memfactor[0] - time_jumpX
        new_Memfactor['MemfactorY'][1][0]  =       0#Memfactor[0] - time_jumpY

        new_Memfactor['MemfactorX'][1][1]  =       (Memfactor[1]%512 - time_jumpX)#X
        if(new_Memfactor['MemfactorX'][1][1] < 0):
            new_Memfactor['MemfactorX'][1][1]  =       (Memfactor[1]%512 - time_jumpX)%512#X
        new_Memfactor['MemfactorY'][1][1]  =       (Memfactor[1]%512 - time_jumpY)#Y
        print(Memfactor[1]%512-time_jumpX, Memfactor[1]%512-time_jumpY)
        if(new_Memfactor['MemfactorY'][1][1]<0):
            new_Memfactor['MemfactorY'][1][1]  =       (Memfactor[1]%512 - time_jumpY)%512


        if( int((Memfactor[1]-time_jumpY)/512) < int(Memfactor[1]/512)):
            tempcomf_Y                                  =          tempcomf_Y[512:]#tempcomf_Y[int(512-(Memfactor[1]%512+time_jumpY))%512:]
            #new_Memfactor['MemfactorY'][1][1]           =          (Memfactor[0]%512-time_jumpY + 512
            print('In int((Memfactor[1]-time_jumpY)/512) < int(Memfactor[1]/512) condition..')
        if( int((Memfactor[1]-time_jumpX)/512) < int(Memfactor[1]/512)):
            tempcomf_X                                  =          tempcomf_X[512:]#tempcomf_X[int(512-(Memfactor[1]%512+time_jumpX))%512:]
            #new_Memfactor['MemfactorX'][1][1]           =          (Memfactor[0]%512-time_jumpX + 512
            print('In int((Memfactor[1]-time_jumpX)/512) < int(Memfactor[1]/512) condition..')
        print(time_jumpY, time_jumpX)
        print('new_Memfactor at the end is..')
        print(new_Memfactor)
    if(time_flag==0):
        print('Time flag is low..generating new fparam file..')
        print(file_name, file_name1, tempcomf_X[int((Memfactor[0]%512)):], tempcomf_Y[int((Memfactor[0]%512)):], tempcomf1_X[int((Memfactor[1]%512)):], tempcomf1_Y[int((Memfactor[1]%512)):], ch, le, RFI)
        fparam_delayX, fparam_delayY	=	_header_gps_cy.gen_finer_shift_parameter(file_name, file_name1, tempcomf_X[int((Memfactor[0]%512)):], tempcomf_Y[int((Memfactor[0]%512)):], tempcomf1_X[int((Memfactor[1]%512)):], tempcomf1_Y[int((Memfactor[1]%512)):], ch, le, RFI) 
        #Fitting a parabolic curve, we have 
        print('fparam are X, Y..'+str(fparam_delayX)+','+str(fparam_delayY))
        
        if(Memfactor[0]!=0):
            new_Memfactor['MemfactorX'][0][0]            =          Memfactor[0]+time_jumpX
            new_Memfactor['MemfactorX'][0][1]            =          0
            new_Memfactor['MemfactorY'][0][0]            =          Memfactor[0]+time_jumpY
            new_Memfactor['MemfactorY'][0][1]            =          0

            new_Memfactor['MemfactorX'][1][0]            =          (Memfactor[0]%512+time_jumpX)
            if(new_Memfactor['MemfactorX'][1][0]<0):
                print("In new_Memfactor['MemfactorX'][1][0]<0 condition")
                new_Memfactor['MemfactorX'][1][0]        =          (Memfactor[0]%512+time_jumpX)%512

            new_Memfactor['MemfactorX'][1][1]            =          0

            new_Memfactor['MemfactorY'][1][0]            =          (Memfactor[0]%512+time_jumpY)
            if(new_Memfactor['MemfactorY'][1][0]<0):
                print("In new_Memfactor['MemfactorY'][1][0]<0 condition")
                new_Memfactor['MemfactorY'][1][0]        =          (Memfactor[0]%512+time_jumpY)%512

            new_Memfactor['MemfactorY'][1][1]            =          0


            if( int((Memfactor[0]+time_jumpY)/512) < int(Memfactor[0]/512)):
                tempcomf1_Y                                  =          tempcomf1_Y[int(512-(Memfactor[0]%512+time_jumpY))%512:]
                new_Memfactor['MemfactorY'][1][0]            =          0#(Memfactor[0]%512+time_jumpY)%512
                print('\n\n\n\n')
            if( int((Memfactor[0]+time_jumpX)/512) < int(Memfactor[0]/512)):
                tempcomf1_X                                  =          tempcomf1_X[int(512-(Memfactor[0]%512+time_jumpX))%512:]
                new_Memfactor['MemfactorX'][1][0]            =          0#(Memfactor[0]%512+time_jumpX)%512
                print('\n\n\n\n')
            print(time_jumpY, time_jumpX)
            print('new_Memfactor at the end is..')
            print(new_Memfactor)

        else:



            new_Memfactor['MemfactorX'][0][0]  =       0#Memfactor[0] - time_jumpX
            new_Memfactor['MemfactorY'][0][0]  =       0#Memfactor[0] - time_jumpY

            new_Memfactor['MemfactorX'][0][1]  =       Memfactor[1] - time_jumpX#X
            new_Memfactor['MemfactorY'][0][1]  =       Memfactor[1] - time_jumpY#Y



            new_Memfactor['MemfactorX'][1][0]  =       0#Memfactor[0] - time_jumpX
            new_Memfactor['MemfactorY'][1][0]  =       0#Memfactor[0] - time_jumpY

            new_Memfactor['MemfactorX'][1][1]  =       (Memfactor[1]%512 - time_jumpX)#X
            if(new_Memfactor['MemfactorX'][1][1] < 0):
                new_Memfactor['MemfactorX'][1][1]  =       (Memfactor[1]%512 - time_jumpX)%512#X
            new_Memfactor['MemfactorY'][1][1]  =       (Memfactor[1]%512 - time_jumpY)#Y
            print(Memfactor[1]%512-time_jumpX, Memfactor[1]%512-time_jumpY)
            if(new_Memfactor['MemfactorY'][1][1]<0):
                new_Memfactor['MemfactorY'][1][1]  =       (Memfactor[1]%512 - time_jumpY)%512


            if( int((Memfactor[1]-time_jumpY)/512) < int(Memfactor[1]/512)):
                tempcomf_Y                                  =          tempcomf_Y[int(512-(Memfactor[1]%512+time_jumpY))%512:]
                new_Memfactor['MemfactorY'][1][1]           =          0#(Memfactor[0]%512+time_jumpY)%512

            if( int((Memfactor[1]-time_jumpX)/512) < int(Memfactor[1]/512)):
                tempcomf_X                                  =          tempcomf_X[int(512-(Memfactor[1]%512+time_jumpX))%512:]
                new_Memfactor['MemfactorX'][1][1]           =          0#(Memfactor[0]%512+time_jumpX)%512
            print(time_jumpY, time_jumpX)
            print('new_Memfactor at the end is..')
            print(new_Memfactor)




        Memfactor	=	Memfactor1
        print('Memfactor...'+str(Memfactor[0])+' , '+str(Memfactor[1]))
    
    print('Time required to decrypt both files..'+str(time.time()-timea))
    print(len(tempcomf_X), len(tempcomf1_X))
    print(new_Memfactor_flag)
    return tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, LO1, file_time, Memfactor1, new_Memfactor 

cdef sub_brute_force(comf_search, val):
    #cdef unsigned int i	=	0
    #for i in range(len(comf_search)):
    #    if(comf_search['GPS'][i] == val):
    #        return i
    return int(np.argmax(comf_search['GPS']> val-1))
cpdef trans_flag_brute_force(comf_end, comf, comf1, Memfact):
    '''
        This module is for brute force compensation..

    '''
    print('Memfact in trans_flag_brute_force is..'+str(Memfact))
    cdef long long int swt	=	int(comf['GPS'][0] > comf1['GPS'][0])
    cdef long long int max_val	=	max(comf['GPS'][0], comf1['GPS'][0])+1
    cdef long long int min_val  =       0
    cdef int Mem_swt		=	0	
    if(swt):
        min_val			=	sub_brute_force(comf, max_val)#np.where(comf['GPS'] == max_val)[0][0]
    else:
        min_val			=	sub_brute_force(comf1, max_val)#np.where(comf1['GPS'] == max_val)[0][0]
    print('Max_val in trans_flag_brute_force..'+str(max_val))
    print('Min_val in trans_flag_brute_force..'+str(min_val))
    cdef long long int Mem1		=	0
    try:
        print('Trying to find the GPS value in modified Memfact..'+str(Memfact))
        Mem1			=	sub_brute_force(comf_end[int(Memfact)], max_val) - 1#np.where(comf_end[int(Memfact)]['GPS']==max_val)[0][0]
    except:
        print('Did not find the  GPS value in modified Memfact..'+str(Memfact)+' hence looking in next file!')
        Mem1                    =       sub_brute_force(comf_end[int(Memfact)+1], max_val)#np.where(comf_end[int(Memfact)+1]['GPS']==max_val)[0][0]
        Mem_swt			=	1
    print('Mem1, Mem_swt..'+str(Mem1)+','+str(Mem_swt))
    return Mem1+(len(comf)-1)*(int(Memfact)+Mem_swt) - min_val  + Mem_swt

cpdef get_baseline(file_name, file_name1):

    try:
        file_name_1 =       file_name.split('/')[-1]
        file_name1_1=       file_name1.split('/')[-1]
    except:
        file_name_1 =       file_name
        file_name1_1=       file_name1
    #Reading Baseline_connect_file.txt file
    
    cdef np.ndarray Baseline_DAS=	np.loadtxt('Baseline_connect_file.txt', dtype=str, delimiter='--')
    cdef str das1		=	(file_name_1.split('_')[0])[0:4]+'X'
    cdef str das2		=	(file_name1_1.split('_')[0])[0:4]+'Y'
   
    print(das1, das2) 
    cdef str cdas1		=	(Baseline_DAS[np.argwhere(Baseline_DAS==das1)[0][0]][1])[2:4]
    cdef str cdas2		=	(Baseline_DAS[np.argwhere(Baseline_DAS==das2)[0][0]][1])[2:4]

    return int(cdas1)-1, int(cdas2)-1

cpdef break_file_name(file_name, file_name1):

    try:
        file_name_1 =       file_name.split('/')[-1]
        file_name1_1=       file_name1.split('/')[-1]
    except:
        file_name_1 =       file_name
        file_name1_1=       file_name1

    return file_name_1, file_name1_1

cdef get_next_file(file_name, series):

    return file_name[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'


@cython.boundscheck(False)
@cython.wraparound(False)
cpdef comf_comf1_compensate(comf, comf1, file_name, file_name1, Memfactor, series):


    '''
        Important flags to remeber!!
        missing_file_comf1_1
        missing_file_comf_1
    '''
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    missing_file_comf1_1		=	0
    missing_file_comf_1			=	0

    cdef long long int ploss1                                      =       comf['Packet'][len(comf)-1]  -  comf['Packet'][0] -len(comf['Packet']) + 1
    cdef long long int ploss2                                      =       comf1['Packet'][len(comf1)-1] -  comf1['Packet'][0] - len(comf1['Packet']) + 1
    cdef long long int internal_jump
    file_name_1, file_name1_1          =       break_file_name(file_name, file_name1)

    Mem1, Mem2, Memfactor              =            get_gps_info_000_file(comf, comf1, Memfactor)
    Memfactor1                         =            gps_alingn_files(ploss1, ploss2, Memfactor[0]/512.0, Memfactor[1]/512.0, file_name_1, file_name1_1, series)

    print(Memfactor1)
    Mem1=int(Memfactor1[0]/512.0)
    Mem2=int(Memfactor1[1]/512.0)
    #Memfactor1                        =    np.array([0,0])
    if(Memfactor1[0]!=0):

        internal_jump                  =    int(Memfactor1[0]/512)# - (comf['Packet'][(int((Memfactor1[0]/512.0)))]   -   comf['Packet'][0]- len(comf[0:(int((Memfactor1[0]/512.0)))])  )#int(int(Memfactor1[0])/512))
        internal_read                  =    (comf['Packet'][(int((Memfactor1[0]/512.0)))]   -   comf['Packet'][0]- len(comf[0:(int((Memfactor1[0]/512.0)))]) )
        internal_read1                 =    (comf1['Packet'][len(comf1)-1]   -   comf1['Packet'][0]- len(comf1) + 1)

        #internal_read		       =    0

        file_next                      =    get_next_file(file_name, series)#file_name[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'
        file_next1                     =    get_next_file(file_name1, series)#file_name1[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'

        try:
            #Trying to read next file, if it even exsists!
            comf1_1                        =    np.memmap(file_next1,  dtype = dt, mode = 'c')[0:internal_read+internal_read1]
        except:
            #Did not find the next file!
            missing_file_comf1_1       =    1
        try:
            #Trying to read next file, if it even exsists!
            comf_1                         =    np.memmap(file_next,  dtype = dt, mode = 'c')[0:int(Memfactor1[0]/512.0)+internal_read]
        except:
            #Did not find the next file!
            missing_file_comf_1       =    1
       


        if(missing_file_comf1_1==0 and missing_file_comf_1 ==0):
            comf                           =    np.concatenate((comf[internal_jump:], comf_1))
            comf1                          =    np.concatenate((comf1, comf1_1))
            print('Extra read packets is..'+str(len(comf_1)))
        elif(missing_file_comf1_1==1 and missing_file_comf_1 ==0):
            comf                           =    np.concatenate((comf[internal_jump:], comf_1))
        elif(missing_file_comf1_1==0 and missing_file_comf_1 ==1):
            comf1                          =    np.concatenate((comf1, comf1_1))
        else:
            comf                           =    comf[internal_jump:]
        print('Internal loss read is..'+str(internal_read))
        print('Memfactor in 000 series..'+str(Memfactor1[0])+','+str(Memfactor1[1]))
        print('Internal GPS+Pack1 loss jump is..'+str(internal_jump))
        print('Total length is comf, comf1..'+str(len(comf))+','+str(len(comf1)))

    else:
        internal_jump                  =            int(Memfactor1[1]/512)# - (comf1['Packet'][(int((Memfactor1[1]/512.0)))] -comf1['Packet'][0] - len(comf1[0:(int((Memfactor1[1]/512.0)))]) )#-int(int(Memfactor1[1])/512))
        internal_read                  =            (comf1['Packet'][(int((Memfactor1[1]/512.0)))] -comf1['Packet'][0] - len(comf1[0:(int((Memfactor1[1]/512.0)))]) )
        internal_read1                 =            (comf['Packet'][len(comf)-1]   -   comf['Packet'][0]- len(comf) + 1)


        #internal_read		       =	    0

        file_next                      =            get_next_file(file_name1, series)#file_name1[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'
        file_next1                     =            get_next_file(file_name, series)#file_name[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'
        
        try: 
            #Trying to read next file, if it even exsists!
            comf1_1                        =            np.memmap(file_next,  dtype = dt, mode = 'c')[0:int(Memfactor1[1]/512.0)+internal_read]
        except:
            #Did not find the next file!
            missing_file_comf1_1       =    1
        try:
            #Trying to read next file, if it even exsists!
            comf_1                         =            np.memmap(file_next1,  dtype = dt, mode = 'c')[0:internal_read+internal_read1]
        except:
            #Did not find the next file!
            missing_file_comf_1       =    1

        print('internal_read1..'+str(internal_read1))

        if(missing_file_comf1_1==0 and missing_file_comf_1 ==0):
            comf                           =            np.concatenate((comf, comf_1))
            comf1                          =            np.concatenate((comf1[internal_jump:], comf1_1))
            print('Extra read packets is..'+str(len(comf1_1)))
        elif(missing_file_comf1_1==1 and missing_file_comf_1 ==0):
            comf                           =            np.concatenate((comf, comf_1))
        elif(missing_file_comf1_1==0 and missing_file_comf_1 ==1):
            comf1                          =            np.concatenate((comf1[internal_jump:], comf1_1))
        else:
             comf1                          =            comf1[internal_jump:]
        print('Internal loss read is..'+str(internal_read))
        print('Memfactor in 000 series..'+str(Memfactor1[0])+','+str(Memfactor1[1]))
        print('Internal GPS+Pack2 loss jump is..'+str(internal_jump))
        print('Total length is comf, comf1..'+str(len(comf))+','+str(len(comf1)))


    return comf, comf1, Memfactor1


cpdef read_spinoff(series, Memfactor, file_name, file_name1, dt):

    time_init				=	time.time()
    file_name_1, file_name1_1		=	break_file_name(file_name, file_name1)

    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])

    cdef np.ndarray comf
    cdef np.ndarray comf1

    comf     = np.memmap(file_name,  dtype = dt, mode = 'c')#[0:]
    comf1    = np.memmap(file_name1, dtype = dt, mode = 'c')#[1968466:]
    LO1                  =   comf['LO'][100]
    LO2                  =   comf1['LO'][100]

    if(LO1!=LO2):
        print('LO1 = ' +str(LO1)+'\n')
        print('LO2 = ' +str(LO2)+'\n')
        #raise RuntimeError ('Both LOs are different!')


    cdef long long int ploss1                                      =       comf['Packet'][-1]  -  comf['Packet'][0] -len(comf['Packet']) + 1
    cdef long long int ploss2                                      =       comf1['Packet'][-1] -  comf1['Packet'][0] - len(comf1['Packet']) + 1
    cdef long long int internal_jump
    cdef list comf_end                       			   =            []
    cdef int trans_flag						   =		0
    cdef int missing_file_flag					   =		0
    if(series=='000'):
        comf, comf1, Memfactor1			   =		comf_comf1_compensate(comf, comf1, file_name, file_name1, Memfactor, series)     
    else:

        #Pre GPS align
        Memfactor			   =		pre_gps_align(file_name_1, file_name1_1)   
        print('Early Memfactor is..'+str(Memfactor[0])+','+str(Memfactor[1]))         
        Memfactor1			   =		gps_alingn_files(ploss1, ploss2, Memfactor[0], Memfactor[1], file_name_1, file_name1_1, series)
        Mem1=int(Memfactor1[0]/512.0)
        Mem2=int(Memfactor1[1]/512.0)
        if(Memfactor1[0]!=0):
            #trying to read next files..#
            Mem_fact			   =		(Memfactor1[0]/512)/2027520
            comf_end.append(np.memmap(file_name,  dtype = dt, mode = 'c'))
            print(series, Mem_fact)
            print(int(series)+1, int(series)+int(Mem_fact)+2)
            print('Mem_fact..'+str(Mem_fact))
            for i in range(int(series)+1, int(series)+int(Mem_fact)+2):
                print('itration is..'+str(i))
                file_next_mem_fact		   =		file_name[:-7]+str(int(i)).rjust(3, '0')+'.mbr'
                try:
                    comf_end.append(np.memmap(file_next_mem_fact,  dtype = dt, mode = 'c'))
                except: 
                    missing_file_flag              =            missing_file_flag+1 
                    print('The file series '+str(i)+'is not found in the given path..')
                
                print(file_next_mem_fact)
            
            #print( comf_end[int(Mem_fact)]['Packet'][int(Memfactor1[0]/512)%2027520], comf_end[0]['Packet'][0], int(Memfactor1[0]/512))
            if(Mem_fact>1 and Mem_fact<2):
                if(missing_file_flag==0):
                    internal_jump_mem_last     =            comf_end[int(Mem_fact)-1]['Packet'][len(comf)-1] - comf_end[0]['Packet'][0] - (len(comf)-1)*(int(Mem_fact)) #- 1*(int(Mem_fact)-1)
                else:
                    internal_jump_mem_last     =            comf_end[int(Mem_fact)-1-missing_file_flag]['Packet'][len(comf)-1] - comf_end[0]['Packet'][0] - (len(comf)-1-missing_file_flag)*(int(Mem_fact))
                print('internal_jump_mem_last is..'+str(internal_jump_mem_last))
                print('In Mem_fact>1 and Mem_fact<2 condition')
                trans_flag			=		int((int(Memfactor1[0]/512) - internal_jump_mem_last)/len(comf_end[0])<int(Mem_fact))
            elif(Mem_fact>2):
                #Check if Mem_fact changes by an interger, for transition files only!
                #Calculating packet loss of first file, to see the translatory effect..
                if(missing_file_flag==0):
                    internal_jump_mem_last	=            comf_end[int(Mem_fact)-1]['Packet'][len(comf_end[0])-1] - comf_end[0]['Packet'][0] - (len(comf_end[0])-1)*int(Mem_fact) -  1*(int(Mem_fact)-1)
                    print(comf_end[int(Mem_fact)-1]['Packet'][len(comf_end[0])-1] , comf_end[0]['Packet'][0] , (len(comf_end[0])-1)*int(Mem_fact) ,  1*(int(Mem_fact)-1))
                else:
                    internal_jump_mem_last      =            comf_end[int(Mem_fact)-1-missing_file_flag]['Packet'][len(comf_end[0])-1] - comf_end[0]['Packet'][0] - (len(comf_end[0])-1)*int(Mem_fact) -  1*(int(Mem_fact)-1-missing_file_flag)
                #Translatory effect found!
                trans_flag			=		int((int(Memfactor1[0]/512) - internal_jump_mem_last)/len(comf)<int(Mem_fact))
                if(trans_flag!=1):
                    internal_jump_mem_last	=		comf_end[int(Mem_fact)-1]['Packet'][len(comf)-1] - comf_end[0]['Packet'][0] - (len(comf)-1)*(int(Mem_fact)) - 1*(int(Mem_fact)-1)
                print('internal_jump_mem_last is..'+str(internal_jump_mem_last))
                print('In Mem_fact>2 condition')
            else:
                internal_jump_mem_last	   	=            0#int(Memfactor1[0]/512)
            if(trans_flag==1):
                #internal_jump_mem              =            comf_end[int(Mem_fact)-1-missing_file_flag-1]['Packet'][len(comf_end[0])-1] - comf_end[0]['Packet'][0] - (len(comf_end[0])-1)*int(Mem_fact) -  1*(int(Mem_fact)-1-missing_file_flag)
                #internal_jump_mem		=		int(Memfactor1[0]/512)  -       internal_jump_mem_last
                #print(comf_end[int(Mem_fact)]['Packet'][internal_jump_mem_last%2027520], comf_end[int(Mem_fact)]['Packet'][0], internal_jump_mem_last%2027520)
                #print('Intermediate internal_jump_mem..'+str(internal_jump_mem))
                
               
                #We are in not so comfortable situation of transitory effect! hence using the brute force, inefficient method of synchronization!
                internal_jump_mem	       =	    trans_flag_brute_force(comf_end, comf, comf1, (int(Memfactor1[0]/512) - internal_jump_mem_last)/len(comf))
            else:
                Memfact_internal_jump          =            int(Memfactor1[0]/512) - internal_jump_mem_last
                print('Memfact_internal_jump is..'+str(Memfact_internal_jump))


                #internal_jump_mem		=	0
                #if(float(Memfact_internal_jump)/len(comf)>1 and float(Memfact_internal_jump)/len(comf)<2):
                internal_jump_mem		   =		comf_end[int(Mem_fact)]['Packet'][Memfact_internal_jump%2027520]-comf_end[int(Mem_fact)]['Packet'][0] - Memfact_internal_jump%2027520 


                print(comf_end[int(Mem_fact)]['Packet'][Memfact_internal_jump%2027520], comf_end[int(Mem_fact)]['Packet'][0], Memfact_internal_jump%2027520)
                print('Intermediate internal_jump_mem..'+str(internal_jump_mem))
                internal_jump_mem		   =		int(Memfactor1[0]/512)	-	internal_jump_mem	-	internal_jump_mem_last	
            print('Length of comf_end is..'+str(len(comf_end)))
            print('internal_jump_mem is..'+str(internal_jump_mem))
            if(missing_file_flag>0):
                comf			   =		comf_end[int(internal_jump_mem/2027520)][internal_jump_mem%2027520:]
            else:
                comf			   =		np.concatenate((comf_end[int(internal_jump_mem/2027520)][internal_jump_mem%2027520:], comf_end[int(internal_jump_mem/2027520)+1][:internal_jump_mem%2027520]))
                print('internal_jump_mem is..'+str(internal_jump_mem))
        else:
            #trying to read next files..#
            Mem_fact                       =            (Memfactor1[1]/512)/2027520
            comf_end.append(np.memmap(file_name1,  dtype = dt, mode = 'c'))
            print(series, Mem_fact)
            print(int(series)+1, int(series)+int(Mem_fact)+2)
            print('Mem_fact..'+str(Mem_fact))
            for i in range(int(series)+1, int(series)+int(Mem_fact)+2):
                print('itration is..'+str(i))
                file_next_mem_fact                 =            file_name1[:-7]+str(int(i)).rjust(3, '0')+'.mbr'
                try:
                    comf_end.append(np.memmap(file_next_mem_fact,  dtype = dt, mode = 'c'))
                except:
                    missing_file_flag              =            missing_file_flag+1
                    print('The file series '+str(i)+' is not found in the given path..')
                
                print('Missing file flag is..'+str(missing_file_flag))
                print(file_next_mem_fact)
            #print( comf_end[int(Mem_fact)]['Packet'][int(Memfactor1[1]/512)%2027520], comf_end[0]['Packet'][0], int(Memfactor1[1]/512))
            if(Mem_fact>1 and Mem_fact<2 and missing_file_flag==0):
                internal_jump_mem_last     =            comf_end[int(Mem_fact)-1]['Packet'][len(comf)-1] - comf_end[0]['Packet'][0] - (len(comf)-1)*(int(Mem_fact)) #- 1*(int(Mem_fact)-1)
                print('internal_jump_mem_last is..'+str(internal_jump_mem_last))
                print('In Mem_fact>1 and Mem_fact<2 condition')
                trans_flag	=	int((int(Memfactor1[1]/512) - internal_jump_mem_last)/len(comf_end[0])<int(Mem_fact))
            elif(Mem_fact>2 and missing_file_flag ==0):
                #Check if Mem_fact changes by an interger, for transition files only!
                #Calculating packet loss of first file, to see the translatory effect..
                internal_jump_mem_last          =            comf_end[int(Mem_fact)-1]['Packet'][len(comf_end[0])-1] - comf_end[0]['Packet'][0] - (len(comf_end[0])-1)*int(Mem_fact) -  1*(int(Mem_fact)-1)
                print(comf_end[int(Mem_fact)-1]['Packet'][len(comf_end[0])-1] , comf_end[0]['Packet'][0] , (len(comf_end[0])-1)*int(Mem_fact) ,  1*(int(Mem_fact)-1))
                trans_flag                  	=	     int((int(Memfactor1[1]/512) - internal_jump_mem_last)/len(comf_end[0])<int(Mem_fact))               
                if(trans_flag!=1):
                    internal_jump_mem_last      =               comf_end[int(Mem_fact)-1]['Packet'][len(comf)-1] - comf_end[0]['Packet'][0] - (len(comf)-1)*(int(Mem_fact)) - 1*(int(Mem_fact)-1)
                print('internal_jump_mem_last is..'+str(internal_jump_mem_last))
                print('In Mem_fact>2 condition')
            else:
                internal_jump_mem_last     =            0#int(Memfactor1[0]/512)
            




            if(trans_flag==1):
                internal_jump_mem                  =            0
                internal_jump_mem                  =            int(Memfactor1[1]/512)  -       internal_jump_mem_last
                print(comf_end[int(Mem_fact)]['Packet'][internal_jump_mem_last%2027520], comf_end[int(Mem_fact)]['Packet'][0], internal_jump_mem_last%2027520)
                print('Intermediate internal_jump_mem..'+str(internal_jump_mem))

            else:
                Memfact_internal_jump          =            int(Memfactor1[1]/512) - internal_jump_mem_last
                print('Memfact_internal_jump is..'+str(Memfact_internal_jump))


                #internal_jump_mem               =       0
                #if(float(Memfact_internal_jump)/len(comf)>1 and float(Memfact_internal_jump)/len(comf)<2):
                internal_jump_mem                  =            comf_end[int(Mem_fact)]['Packet'][Memfact_internal_jump%2027520]-comf_end[int(Mem_fact)]['Packet'][0] - Memfact_internal_jump%2027520 


                print(comf_end[int(Mem_fact)]['Packet'][internal_jump_mem_last%2027520], comf_end[int(Mem_fact)]['Packet'][0], internal_jump_mem_last%2027520)
                print('Intermediate internal_jump_mem..'+str(internal_jump_mem))
                internal_jump_mem                  =            int(Memfactor1[1]/512)  -       internal_jump_mem       -       internal_jump_mem_last
            print('Length of comf_end is..'+str(len(comf_end)))
            print('internal_jump_mem is..'+str(internal_jump_mem))
           
             
            if(missing_file_flag>0):
                comf1                           =            comf_end[int(internal_jump_mem/2027520)][internal_jump_mem%2027520:]#np.concatenate((comf_end[int(internal_jump_mem/2027520)-missing_file_flag][internal_jump_mem%2027520:], comf_end[int(internal_jump_mem/2027520)+1-missing_file_flag][:internal_jump_mem%2027520]))
            else:
                comf1                           =            np.concatenate((comf_end[int(internal_jump_mem/2027520)][internal_jump_mem%2027520:], comf_end[int(internal_jump_mem/2027520)+1][:internal_jump_mem%2027520]))
            print('internal_jump_mem is..'+str(internal_jump_mem))
 
        #print('Concatenating the second file with..'+str(len(comf1_1)))
        #print('Internal GPS+Pack loss jump is..'+str(internal_jump))
        print('Memfactor..'+str(Memfactor[0])+','+str(Memfactor[1]))    
    print('Time till getting Memfactor..'+str(time.time()-time_init))


    cdef float dploss1		#			=		tot_len_end-tot_len_start
    cdef float dploss2		#			=		tot_len_end1-tot_len_start1
    Mem1=int(Memfactor1[0]/512.0)
    Mem2=int(Memfactor1[1]/512.0)

    Memfactor                    =  Memfactor1
    print('After Memfactor change..Memfactor[0], Memfactor[1]..'+str(Memfactor[0])+','+str(Memfactor[1]))
   
    
    time_a					=		time.time()
    #dploss					=		(dploss1-Memfactor[0]/512.0)-(dploss2-Memfactor[1]/512.0)
    
    #Adjusting for the packet loss..        
    #print('dploss1 and dploss2..'+str(dploss1)+','+str(dploss2)+','+str(dploss))
    #cdef np.ndarray lin1       		=   	comf['Packet'][0:]  -comf['Packet'][0]
    #cdef np.ndarray lin2       		=	comf1['Packet'][0:] -comf1['Packet'][0]

    tim_allot                      	=   	time.time()
        
    #cdef long long int templen	   	=	comf['Packet'][-1]  -comf['Packet'][0]+10
    #cdef long long int templen1    	=   	comf1['Packet'][-1] -comf1['Packet'][0]+10

    ##Allocoating memory
    ##tempcomf, tempcomf1 		=	mem_alloc(templen, templen1)
    
    #cdef np.ndarray 	tempcomf	=   	np.zeros((templen, 1024), dtype = np.int8)
    #cdef np.ndarray	tempcomf1	=	np.zeros((templen1, 1024), dtype = np.int8)
     
    print('Time from Memfactor to memmap_copy...'+str(time.time()-time_a))
    cdef np.ndarray file_time		
    file_time				=	get_last_gps_timing(comf, comf1, file_name, file_name1, series, Memfactor)
    tempcomf, tempcomf1			=	memmap_copy(comf, comf1, Mem1, Mem2)
    #tempcomf, tempcomf1     		=	compensate_pack_loss(comf, comf1, tempcomf, tempcomf1, lin1, lin2, Mem1, Mem2)
    #tempcomf, tempcomf1			=	tempcomf.ravel(), tempcomf1.ravel()
    dploss1				=	len(tempcomf)
    dploss2				=	len(tempcomf1)
    dploss				=	(dploss1 - Memfactor[0])-(dploss2 - Memfactor[1])
    dploss				=	dploss/512.0#+sum(Memfactor)/1024.0

    print('Current dploss is...'+str(dploss))
    #with cython.boundscheck(False):
    #    #Now copying the memory mapped view to volatile memory
    #    tempcomf[lin1]          =   np.memmap.copy(comf['data'][0:])
    #    tempcomf1[lin2]         =   np.memmap.copy(comf1['data'][0:])
    #print('Time for array array allotment...'+str(time.time()-tim_allot))
    
    #tempcomf, tempcomf1              = compensate_pack_loss(comf, comf1, tempcomf, tempcomf1, lin1, lin2, Mem1, Mem2)#.ravel()
    #tempcomf, tempcomf1             = tempcomf.ravel(), tempcomf1.ravel()

    #Now getting the number of packets to be read from the next file..to compensate for the extra length of the data arising
    #from the correlation GPS compensation..
    print('tempcomf and tempcomf1 length..'+str(len(tempcomf))+','+str(len(tempcomf1))+' Memfactor..'+str(Memfactor1[0])+','+str(Memfactor1[1]))


    print('The length of tempcomf...'+str(len(tempcomf)))
    print('The length of tempcomf1...'+str(len(tempcomf1)))
    print('The length of dploss..'+str(dploss))

    #tempcomf, tempcomf1			=	memmap_copy(comf, comf1, Mem1, Mem2)
    #tempcomf, tempcomf1                 =       tempcomf.ravel(), tempcomf1.ravel()       
    #if(dploss>0):
    #    print('Now reading file2 for extra packets..')
    #    file_next   			=       file_name1[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'
    #    tempcomf_rem1			=	read_comf_rem(file_next, dt, int(round(abs(dploss))))
    #    tempcomf1			=	np.concatenate((tempcomf1, tempcomf_rem1[0:]))#np.vstack((comf1, tempcomf_rem[0:abs(dploss)]))
    #    print(len(tempcomf_rem1))
    #else:
    #    print('Now reading file1 for extra packets..')
    #    file_next           		=       file_name[:-7]+str(int(series)+1).rjust(3, '0')+'.mbr'
    #    tempcomf_rem1       		=       read_comf_rem(file_next, dt, int(round(abs(dploss))))
    #    tempcomf	   		=       np.concatenate((tempcomf, tempcomf_rem1[0:]))#np.vstack((comf1, tempcomf_rem[0:abs(dploss)]))
    #    print(len(tempcomf_rem1))



    time_allot				=   time.time() 
    tempcomf_X                          =   np.array(tempcomf[1::2],  order = 'F')
    tempcomf_Y                          =   np.array(tempcomf[0::2],  order = 'F')
    tempcomf1_X                         =   np.array(tempcomf1[1::2], order = 'F')
    tempcomf1_Y                         =   np.array(tempcomf1[0::2], order = 'F')
    print('Time for Pol array allotment...'+str(time.time()-tim_allot))
    print(len(tempcomf_X), len(tempcomf1_X))
    return tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, Memfactor, Mem1, Mem2, LO1, LO2, file_time

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef memmap_copy(comf, comf1, Mem1, Mem2):
    time_allot		    =   time.time()
    
    cdef np.ndarray lin1                =       np.zeros((len(comf)), dtype = np.int32)
    cdef np.ndarray lin2                =       np.zeros((len(comf1)), dtype = np.int32)
    print('Time to allot lin1 and lin2...'+str(time.time()-time_allot))
    cdef int i
    #for i in range(len(comf)): 
    lin1				=	np.memmap.copy(comf['Packet'][0:len(comf)])# comf['Packet'][i]#np.memmap.copy(comf['Packet']) 
    lin2                		=	np.memmap.copy(comf1['Packet'][0:len(comf1)]) #comf1['Packet'][i]#np.memmap.copy(comf1['Packet'])
    print('Time to allot lin1 and lin2...'+str(time.time()-time_allot))
    lin2				=	lin2 - comf1['Packet'][0]
    lin1				=	lin1 - comf['Packet'][0]
    

    cdef long long int templen          =       comf ['Packet'][len(comf ['Packet'])-1] -comf ['Packet'][0]+10
    cdef long long int templen1         =       comf1['Packet'][len(comf1['Packet'])-1] -comf1['Packet'][0]+10
    
    print('templen, templen1')
    print(templen, templen1) 
    #Allocoating memory
    #tempcomf, tempcomf1                =       mem_alloc(templen, templen1)

    cdef np.ndarray     tempcomf        =       np.zeros((templen, 1024),  dtype =np.int8)
    cdef np.ndarray     tempcomf1       =       np.zeros((templen1, 1024), dtype =np.int8)

    time_allot              =   time.time()

    #Now copying the memory mapped view to volatile memory
    tempcomf[lin1]          =   np.memmap.copy(comf['data'][0:])
    tempcomf1[lin2]         =   np.memmap.copy(comf1['data'][0:])
    

    tempcomf, tempcomf1     = compensate_pack_loss(comf, comf1, tempcomf, tempcomf1, lin1, lin2, 0, 0)#Mem1, Mem2)
    #tempcomf, tempcomf1	    =   tempcomf.ravel(), tempcomf1.ravel() 
    print('Time for copying memmap array to numpy array...'+str(time.time()-time_allot))
    return tempcomf, tempcomf1

'''
cpdef mem_alloc(templen, templen1):

    tempcomf =   np.zeros((templen, 1024), dtype = np.int8)
    tempcomf1=   np.zeros((templen1, 1024), dtype = np.int8)
    cdef int [:,:] tempcomf_view	=	tempcomf
    cdef int [:,:] tempcomf_view1       =       tempcomf1
    cdef int i
    with cython.boundscheck(False):
        with nogil:
            print('Time for array array allotment in 000 series...'+str(time.time()-tim_allot))

            #Now copying the memory mapped view to volatile memory
            for i in lin1:
                tempcomf_view[i] = comf['data'][i]

        tempcomf[lin1]          =   np.memmap.copy(comf['data'][0:])
        tempcomf1[lin2]         =   np.memmap.copy(comf1['data'][0:])
    print('Time for array array allotment...'+str(time.time()-tim_allot))

    #tempcomf, tempcomf1              = compensate_pack_loss(comf, comf1, tempcomf, tempcomf1, lin1, lin2, Mem1, Mem2)#.ravel()
    tempcomf, tempcomf1             = tempcomf.ravel(), tempcomf1.ravel()
'''


cpdef read_comf_rem(tempcomf_rem):#, dt, dploss):
    
         
    #tempcomf_rem					=	np.memmap(file_next, dtype=dt, mode='c')[0:abs(dploss)]	
    cdef np.ndarray linr				=	np.zeros((len(tempcomf_rem)), dtype=np.int32)
    linr[0:len(tempcomf_rem)]				=	tempcomf_rem['Packet'][0:]	-	tempcomf_rem['Packet'][0]
    cdef np.ndarray tempcomf_rem1			=	np.zeros(((tempcomf_rem['Packet'][-1]-tempcomf_rem['Packet'][0]+10, 1024)), dtype = np.int8)
    tempcomf_rem1[linr]					=	tempcomf_rem['data']
    #tempcomf_rem1					=	tempcomf_rem1.ravel()
    return tempcomf_rem1
        

cpdef pre_gps_align(file_name, file_name1):

    try:
        file_name_1 =       file_name.split('/')[-1]
        file_name1_1=       file_name1.split('/')[-1]
    except:
        file_name_1 =       file_name
        file_name1_1=       file_name1

    cdef str file_000        =       'SAMPLING_INFO/Packet_info_'+file_name_1[:-7]+'000.mbr'
    cdef str file_000_1      =       'SAMPLING_INFO/Packet_info_'+file_name1_1[:-7]+'000.mbr'
    cdef np.ndarray pack     =       np.loadtxt(file_000)
    cdef np.ndarray pack1    =       np.loadtxt(file_000_1)
    Memfactor    			 =       np.array([pack[1], pack1[1]])
    return Memfactor

cpdef comf_read_rem_file(file_next, val):
    
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    try:
        comf_rem    = np.memmap(file_next,  dtype = dt, mode = 'c')
    except:
        print('\n\n\n\n No '+str(file_next)+ ' found in the given location..\n\n\n\n')
        return np.array(['N'])
    len1_fil1    =   comf_rem['Packet'][0]
    len2_fil1    =   comf_rem['Packet'][int(round(val))]

    start_point1 =   len(comf_rem['data'])
    templen      =   len2_fil1  -len1_fil1
    print('value to read..'+str(templen))
    #print('len2_fil1-len1_fil1 '+str(len2_fil1)+'-'+str(len1_fil1)+'='+str(templen))

    cdef np.ndarray  lin1       =   comf_rem['Packet'][0:val]  -len1_fil1
    cdef float       tim_allot  =   time.time()
    cdef np.ndarray tempcomf_rem=   np.zeros((templen, 1024), dtype = np.int8) 
    print('Time for array array allotment...'+str(time.time()-tim_allot))
    print(len(comf_rem['data'][0:val]), len(tempcomf_rem), len(lin1))
    tempcomf_rem[lin1]           =   np.memmap.copy(comf_rem['data'][0:val])
    #tempcomf                     =   compensate_pack_loss_rem(comf_rem, tempcomf_rem, lin1, val)#compensate_pack_loss(comf_rem, comf_rem, tempcomf_rem, tempcomf_rem, lin1, lin1, 0, 0)    
    tempcomf			 =	tempcomf_rem.ravel()    

    return tempcomf



cpdef gps_alingn_files(double ploss1, double ploss2, double pjump1, double pjump2, file_name, file_name1, series):
  
    try:
        file_name_1 =       file_name.split('/')[-1]
        file_name1_1=       file_name1.split('/')[-1]
    except:
        file_name_1 =       file_name
        file_name1_1=       file_name1


    

    print('pjump1, pjump2..'+str(pjump1)+', '+str(pjump2))
    cdef double pjump=(pjump1-pjump2)
    
    cdef double pre_losst1=0 
    cdef double pre_losst2 =0
    cdef double pre_jump =0
    #cdef double pre_jump2 =0
    cdef np.ndarray packet 
    cdef np.ndarray packet1
    cdef double Nxt_file_jump=0 
    if(series!='000'):
        if(pjump>0):
            np.savetxt('SAMPLING_INFO/Packet_info_'+file_name_1, np.array([ploss1, pjump]))
            np.savetxt('SAMPLING_INFO/Packet_info_'+file_name1_1, np.array([ploss2, 0]))
        else:
            np.savetxt('SAMPLING_INFO/Packet_info_'+file_name_1, np.array([ploss1, 0]))
            np.savetxt('SAMPLING_INFO/Packet_info_'+file_name1_1, np.array([ploss2, abs(pjump)]))

        
        
        for i in range(int(series)):
            file_pre    =       'SAMPLING_INFO/Packet_info_'+file_name_1[:-7]+str(int(series)-i-1).rjust(3, '0')+'.mbr'
            file_pre1   =       'SAMPLING_INFO/Packet_info_'+file_name1_1[:-7]+str(int(series)-i-1).rjust(3, '0')+'.mbr'
            #print('Now reading..') 
            #print(file_pre, file_pre1)


            packet          =       np.loadtxt(file_pre)
            packet1         =       np.loadtxt(file_pre1)
            pre_losst1      =       packet[0]+pre_losst1
            pre_losst2      =       packet1[0]+pre_losst2
            #print('packet[0], packet1[0], pre_losst1, pre_losst2')
            #print(packet[0], packet1[0], pre_losst1, pre_losst2)
            #pre_jump1       =       packet[1]
            #pre_jump2       =       packet1[1]

            #if(pre_jump1>0):
            #    pre_jump=pre_jump1
            #else:
            #    pre_jump=pre_jump2
            #Nxt_file_jump   =    Nxt_file_jump+  pre_loss1-pre_loss2
            #print('Nxt_file_jump, pre_loss1, pre_loss2, pjump')
            #print(Nxt_file_jump, pre_loss1, pre_loss2, pjump)
    pre_jump=pre_losst1-pre_losst2
    print('pre_jump..'+str(pre_jump))
    print('pre_losst1, pre_losst2..'+str(pre_losst1)+' , '+str(pre_losst2))
    cdef double time_jump                =    0
    cdef int time_jump_flag              =    1

    if(series=='000'):
        return gps_alingn_000_file(ploss1, ploss2, pjump1, pjump2, file_name, file_name1)
    
    cdef double Nxt_file_jump1           =      0
    if(pjump1 == 0):
        Nxt_file_jump1			 =	pre_jump + abs(pjump)# + (pjump%1)*int(series)#pre_jump
        print('The pjump1 is equal to zero, i.e first file has no file jump..hence adding pjump..')
        print(pre_jump, int(abs(pjump)), (pjump%1)*int(series))
    else:
        print('The pjump2 is equal to zero, i.e second file has no file jump..hence subtracting pjump')
        print(pre_jump, int(abs(pjump)), (pjump%1)*int(series))
        Nxt_file_jump1                   =      pre_jump - abs(pjump)#+(pjump%1)*int(series))
    #Now looking for finer adjustments in the absolute time shifts..in the fparam file..
    #time_jump, time_jump_flag		 =    get_fparam(file_name_1, file_name1_1)
    print('Final jump of the file is..'+str(Nxt_file_jump1))
    print('finer jump is..'+str(time_jump))
    print('Initial GPS jump is...'+str(pjump))
    if(pre_jump<0):
        print('Returning..'+str([(abs(Nxt_file_jump1))*512+time_jump, 0]))
        return np.array([(abs(Nxt_file_jump1))*512, 0])
    else:
        print('Returning..'+str([(abs(Nxt_file_jump1))*512+time_jump, 0]))
        return np.array([0, abs(Nxt_file_jump1)*512])

cpdef get_fparam(file_name_1, file_name1_1):
    cdef double time_jumpX               =    0
    cdef double time_jumpY               =    0
    cdef int time_jump_flag              =    1
    print('In get_fparam') 
    try:
        time_jumpX       =       np.loadtxt('SAMPLING_INFO/ParameterX_'+str(file_name_1[:-7])+'000.mbr_'+str(file_name1_1[:-7])+'000.mbr.fparam')
        time_jumpY       =       np.loadtxt('SAMPLING_INFO/ParameterY_'+str(file_name_1[:-7])+'000.mbr_'+str(file_name1_1[:-7])+'000.mbr.fparam')
        print('Finer Time JumpX from the .fparm file is..'+str(time_jumpX))
        print('Finer Time JumpY from the .fparm file is..'+str(time_jumpY))
        time_jump_flag	=	1
        time_jumpX 	=       time_jump_flag*time_jumpX
        time_jumpY	=	time_jump_flag*time_jumpY
        return time_jumpX, time_jumpY, time_jump_flag
    except:
        print('Did not find the file..'+'SAMPLING_INFO/Parameter<X/Y>_'+str(file_name_1[:-7])+'000.mbr_'+str(file_name1_1[:-7])+'000.mbr.fparam')
        print('Hence looking in reverse direction')
        pass;
    try:
        time_jumpX       =       np.loadtxt('SAMPLING_INFO/ParameterX_'+str(file_name1_1[:-7])+'000.mbr_'+str(file_name_1[:-7])+'000.mbr.fparam')
        time_jumpY       =       np.loadtxt('SAMPLING_INFO/ParameterY_'+str(file_name1_1[:-7])+'000.mbr_'+str(file_name_1[:-7])+'000.mbr.fparam')
        print('Finer Time JumpX from the .fparm file is..'+str(time_jumpX))
        print('Finer Time JumpY from the .fparm file is..'+str(time_jumpY))
        time_jump_flag  =       -1
        time_jumpX       =       time_jump_flag*time_jumpX # If time jump file found in the second file as positive then the second file need to jump in reverse direction
        time_jumpY       =       time_jump_flag*time_jumpY
        return time_jumpX, time_jumpY, time_jump_flag
    except:
        print('Finer measurments of time jump, .fparam file not found!!')
        print('Make sure paths are in order..if there is any doubt about the code..\n\n\n')
        time_jump_flag  =       0
    return time_jumpX, time_jumpY, time_jump_flag

cpdef gps_alingn_000_file(ploss1, ploss2, pjump1, pjump2, file_name, file_name1):
   

    try:
        file_name_1 =       file_name.split('/')[-1]
        file_name1_1=       file_name1.split('/')[-1]
    except:
        file_name_1 =       file_name
        file_name1_1=       file_name1

    print('pjump1 and pjump2 in gps_alingn_000_file is..'+str(pjump1)+','+str(pjump2))
    cdef double pjump            =   (pjump1 - pjump2)#/512.0
    cdef np.ndarray    packet1   =   np.zeros((2), dtype=float)
    cdef np.ndarray    packet2   =   np.zeros((2), dtype=float)
    print('In 000 series file, the inital GPS jump is '+str(pjump))
    print('pjump calculated in the gps_align_000_file is..'+str(pjump)) 
    cdef double time_jump                =      0.0
    cdef int time_jump_flag              =      0
    #time_jump, time_jump_flag            =    get_fparam(file_name_1, file_name1_1)
    print('Time jump for 000 series files from .fparam file is..'+str(time_jump))

    
    #Previously generated files might not be the suitable in all senarios..hence
    #Not using this block..
    #Below commented region is related to that..


    #try:
    #    #Searching for reviously save files..
    #    packet1      =       np.loadtxt('SAMPLING_INFO/Packet_info_'+file_name_1)#, packet1)
    #    packet2      =       np.loadtxt('SAMPLING_INFO/Packet_info_'+file_name1_1)#, packet2)
    #    print('packet1..packet2..')
    #    print(packet1)
    #    print(packet2)
    #    if(packet1[1]>0):
    #        return np.array([packet1[1]*512.0+time_jump, packet2[1]]), time_jump_flag
    #    else:
    #        return np.array([packet1[1], packet2[1]*512.0-time_jump]), time_jump_flag
    #except:
    #    print('\n\nNo previously generated synchronization data exsists..building new one\n\n')
    #    pass;
    

    if(pjump<0):
        #If negative then the second file needs longer jump, hence if -ve then file 2 jump..
        packet2[0]  =   ploss2
        packet2[1]  =   abs(pjump)
        packet1[0]  =   ploss1
        packet1[1]  =   0
        print('In gps 000 series module..Saving packet1, packet2')
        print(packet1, packet2)
        print('pjump..'+str(pjump*512)+' time_jump..'+str(time_jump))

        np.savetxt('SAMPLING_INFO/Packet_info_'+file_name_1, packet1)
        np.savetxt('SAMPLING_INFO/Packet_info_'+file_name1_1, packet2)
        return np.array([0, abs(pjump)*512.0])#, time_jump_flag
    else:
        #If positive then the first file needs longer jump and second zero jump..
        packet1[0]  =   ploss1
        packet1[1]  =   pjump
        packet2[0]  =   ploss2
        packet2[1]  =   0
        print('Saving packet1, packet2')
        print(packet1, packet2)
        print('pjump..'+str(pjump*512)+' time_jump..'+str(time_jump))

        np.savetxt('SAMPLING_INFO/Packet_info_'+file_name_1, packet1)
        np.savetxt('SAMPLING_INFO/Packet_info_'+file_name1_1, packet2)
        return np.array([abs(pjump)*512.0, 0])#, time_jump_flag


cpdef tuple decrypy_file_new_SWAN_without_compensation(file_name, file_name1, ch):


    '''
                ch should be the channel number of first file..

		Takes the read file and sorts the X and Y polarizartion in the file into

		comf_X, comf1_X, comf_Y, comf1_Y.

    '''
    #If available get the earlier data set#
    cdef str series      = file_name[-7:-4]
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    if(series   !=  '000'):
        rem     =   open('SAMPLING_INFO/Info_on_straight_line'+str(file_name.split('/')[-1]), 'r').readlines()#open('SAMPLING_INFO/'+str(file_name.split('/')[-1])+'_rem.data', 'r').readlines()
        print(rem)
        rem_l   =   rem[-2].split(',')
        rem1     =   open('SAMPLING_INFO/Info_on_straight_line'+str(file_name1.split('/')[-1]), 'r').readlines()
        rem_l1   =   rem[-2].split(',')

        if(rem_l[0] =='Y'):
            comf_rem             = np.memmap(rem_l[1], dtype = dt, mode='c')
            comf_rem             = comf_rem[:int(float(rem_l[2][:-1]))]

        elif(rem_l1[0] == 'Y'):
            comf_rem             = np.memmap(rem_l1[1], dtype = dt, mode='c')
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



    if(series != '000' and file_name[-39:-35] == rem_l[1][-39:-35] and rem_l[0] == 'Y'):
        print('\n\n\n\n\n\n\n\nIn 1\n\n\n\n')
        comf    =   np.hstack((comf, comf_rem))
    elif(series != '000' and file_name1[-39:-35] == rem_l[1][-39:-35] and rem_l1[0] == 'Y'):
        print('\n\n\n\n\n\n\n\nIn 2\n\n\n\n')
        comf1   =   np.hstack((comf1, comf_rem))
    
    
    
    print(len(comf), len(comf1))
    LO1                  =   comf['LO'][100]
    LO2                  =   comf1['LO'][100]

    if(LO1!=LO2):
        print('LO1 = ' +str(LO1)+'\n')
        print('LO2 = ' +str(LO2)+'\n')
        #raise RuntimeError ('Both LOs are different!')


    #cdef long long int le       = min(len(comf)/avg, (len(comf1)/avg))
    cdef long long int templen    =   0
    cdef long long int templen1   =   0

    ft      =   np.dtype('>i1')
    #fopen   =   open('START_file_'+str(fil_tag), 'a')
    st      =   max(comf['GPS'][0], comf1['GPS'][0])
    op      =   min(comf['GPS'][-1], comf1['GPS'][-1])
    #fopen.write('\nStarting_GPS\n'+str(st)+'\nEnding GPS = \t'+str(op))
    #fopen.close()

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


    tempcomf    =   np.memmap.copy(comf['data'][10:])#data_temp#np.frombuffer(data_temp, dtype=ft, count=1024)

    print('Time required to decrypt one file..'+str(time.time()-timea))
    tempcomf1   =   np.memmap.copy(comf1['data'][10:])#data_temp#np.frombuffer(data_temp, dtype=ft, count=1024)
    tempcomf    = tempcomf.ravel()
    tempcomf1   = tempcomf1.ravel()


    tempcomf_X  			=   np.array(tempcomf[1::2], order = 'F')
    tempcomf_Y  			=   np.array(tempcomf[0::2],  order = 'F')
    tempcomf1_X 			=   np.array(tempcomf1[1::2], order = 'F')
    tempcomf1_Y 			=   np.array(tempcomf1[0::2], order = 'F')


    print('Time required to decrypt both files..'+str(time.time()-timea))
    #GPS_len =   [op, st]
    return tempcomf_X, tempcomf_Y, tempcomf1_X, tempcomf1_Y, LO1



cpdef fix_time(sec, mi, hr, day):
    cdef int trk =   0
    cdef int fact=   1
    

    hour    =   int(sec/3600.0)+hr
    mint    =   (sec/3600.0%1)*60+mi
    minu    =   int((sec/3600.0%1)*60+mi)
    sece    =   int((mint%1)*60)

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

cpdef shift_geo(np.ndarray comf_X, np.ndarray comf_Y, np.ndarray comf1_X, np.ndarray comf1_Y, dela):
    cdef i                         =   0
    cdef list   del_loc            =   []
    cdef double delay_r1           =    0
    cdef double delay_trc          =    0
    for i in range(len(dela)):

        if(delay_r1 > 1.0):
            del_loc.append([i, int(delay_r1)])
            print(del_loc)
            delay_r1   =   delay_r1-int(delay_r1)
        else:
             delay_r1         =  delay_r1    +    dela[i]
    cdef np.ndarray loc     =   np.linspace(0, len(comf_X)-1, len(del_loc), dtype=int)
    #del_loc =   np.array(del_loc)
    #cdef np.ndarray comf_X_1    =   comf_X#np.zeros((len(comf_X))+int(sum(delay)), dtype = np.int8)
    #cdef np.ndarray comf_Y_1    =   comf_Y#np.zeros((len(comf_X))+int(sum(delay)), dtype = np.int8)
    #cdef np.ndarray comf1_X_1   =   comf1_X#np.zeros((len(comf1_X))+int(sum(delay)), dtype = np.int8)
    #cdef np.ndarray comf1_Y_1   =   comf1_Y#np.zeros((len(comf1_X))+int(sum(delay)), dtype = np.int8)
    
    #cdef double delay_r1           =    0
    #cdef double delay_trc          =    0
    
    #print(dela)
    #cdef np.ndarray delay          =    np.zeros((len(dela)), dtype = np.float64)      
    delay                          =    dela #* 33000000 * 10**-6
    #print(delay)
    for i in range(len(del_loc)):
        #if(delay_r1 > 1.0):
            print('In delay_resv')
            #print(int(delay_r1))
            #delay_trc                 =   loc[i]  +   sum(delay[:i])
            comf_X                    =   np.insert(comf_X, int(del_loc[i][0]) , np.zeros((int(del_loc[i][1])), dtype=np.int8))    
            comf_Y                    =   np.insert(comf_Y, int(del_loc[i][0]) , np.zeros(int((del_loc[i][1])), dtype=np.int8))
            #print(i)
            #print(len(comf_X), len(comf_Y))
            #comf_X_1[int(loc[i]+delay[i])]   =   np.zeros((int(delay_r1)), dtype=np.int8)
            #comf_Y_1[int(loc[i]+delay[i])]   =   np.zeros((int(delay_r1)), dtype=np.int8)
            #delay_r1   =   delay_r1-int(delay_r1)
        #else:
            #delay_r1         =  delay_r1    +    delay[i] 
    return comf_X, comf_Y, comf1_X, comf1_Y



cpdef call_to_read(str file_name, str file_name1, avg, sysargv, Memfactor=[0,0]):#, number, number1):
    print('Sync factor..'+str(sysargv))

    if(sysargv==str(1)):
        comf_X, comf_Y, comf1_X, comf1_Y,LO, file_time, Memfactor, new_Memfactor        =       decrypy_file_new_SWAN_onhold(file_name, file_name1, avg, Memfactor)
        #comf_X, comf_Y, comf1_X, comf1_Y, LO, file_time, Memfactor, new_Memfactor
        print('__________________________________________________')
        print(len(comf_X), len(comf_Y), len(comf1_X), len(comf1_Y))
        print('__________________________________________________')
        #comf_X, comf_Y, comf1_X, comf1_Y        =       shift_geo(comf_X, comf_Y, comf1_X, comf1_Y, delay)
    if(sysargv==str(2)):
        print('MBRDSP without packet compensation..')
        #comf    =   comf.split('MBRDSP')[1:]
        #comf1   =   comf1.split('MBRDSP')[1:]
        #comf_X, comf_Y, comf1_X, comf1_Y        =       decrypy_file_new_packet_comp(comf, comf1, avg)
    if(sysargv==str(0)):
        comf_X, comf_Y, comf1_X, comf1_Y,LO        =       decrypy_file_new_SWAN_without_compensation(file_name, file_name1, avg)
        Memfactor                                  =       -1
    return comf_X, comf_Y, comf1_X, comf1_Y, LO, file_time, Memfactor, new_Memfactor

def gen_RFI_matrix(str file_name, str file_name1):
    '''
        CAUTION: Should be run only one, not for every file..time constraint :/

    '''
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])

    comf_X, comf_Y, comf1_X, comf1_Y,LO        =       decrypy_file_new_SWAN_onhold(file_name, file_name1, 60)
    #Memfactor                                  = _header_Fring_cy.gps_sync(file_name, file_name1, 1)
    creal1, creal2, creal3, creal4creal8, creal9= internal_loop_RFI.external_loop(comf_X_1, comf_Y_1, comf1_X_1, comf1_Y_1, avg, le, 255)
    RFI_X1  =   RFI_Reject(creal1, 2, 256, 60)
    RFI_Y1  =   RFI_Reject(creal2, 2, 256, 60)
    RFI_X2  =   RFI_Reject(creal3, 2, 256, 60)
    RFI_Y2  =   RFI_Reject(creal4, 2, 256, 60)
    
    return RFI_X1, RFI_Y1, RFI_X2, RFI_Y2


cdef np.ndarray rms(a,meana, fftsize):
        cdef np.ndarray rmsa       =       np.zeros(int(fftsize), dtype='double')
        cdef int i              =       0
        for i in range(int(fftsize)):
                #for j in range(len(a[0])):
                 rmsa[i] = np.std(a[i])#np.sqrt(np.mean((a[i] - meana[i])**2))#    rmsa[i] += (a[i][j] - meana[i])**2
                #rmsa[i]=sqrt(rmsa[i]/len(a[0]))
        return rmsa



cpdef np.ndarray RFI_Reject(spec, sig, fftsize, avg):
    '''
    RFI Rejection module,
    input:      spec, sigma deviation, fftsize
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
    efficiency_x                =       SNR#/SNR2
    for i in range(256):
        if(efficiency_x[i] > (np.mean(efficiency_x)+int(sig)*np.std(efficiency_x)) or  efficiency_x[i] < (np.mean(efficiency_x)-int(sig)*np.std(efficiency_x))):
            #spec[i]   =   np.zeros((len(spec[0])))
            RFI_list.append(i)
            FLAGS[i]    =   0   
    
    
    return FLAGS



cpdef phase_compensation(spect, delay, freq):
    '''
        Compensation for intra sample delay.
        Input:  (spectrum, delay, frequency / MHz)
        Output: compensated spectrum, phase compensation spectrum
    '''


    cdef np.ndarray   pha =   spect.copy()
    cdef np.ndarray   comp=   spect.copy()
    cdef np.ndarray   f   =   np.linspace(freq-16.5/2, freq+16.5/2, 256)
    cdef int          i   =   0

    for i in range(256):
        for j in range(len(spect[0])):
            pha[i][j] = (np.exp(complex(0, -2*np.pi*f[i]*delay[j])))

    comp    =   pha*spect

    return comp, pha

cdef genphase(spec, delay):
    cdef np.ndarray pha     =    spec.copy()
    cdef int i, j
    for i in range(256):
        for j in range(len(test[0])):
            pha[i][j] = (exp(complex(0, -2*pi*f[i]*delay[j])))
    return pha

cdef plot_all(file_name, file_name1):
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


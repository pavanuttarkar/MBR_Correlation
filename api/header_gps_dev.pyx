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
import single_file_combined_RFInoRFI as internal_loop_without_RFI_single_file
import internal_loop_RFI
#from __future__ import division

'''
0 1 2 3 4 5
Y X Y X Y X
'''

cpdef sub_gen_finer_shift_parameter(fil1, fil2, comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI):
    cdef double delay_param = 0#delay[-1]
    print('In sub_gen_finer_shift_parameter')
    creal1, creal2, creal3, creal4, creal8, creal9  = internal_loop_RFI.external_loop(comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI[0], RFI[2], RFI[1], RFI[3])
    delayX              =    obsdelay_param(creal8)
    delayY              =    obsdelay_param(creal9)
    delay_paramX        =    delayX/5000.0# The second zeros here is insignificant, added due to wierd Cython only error woth savetxt
    delay_paramY        =    delayY/5000.0
    print('\n\n\n\nNewly generated delayX finer parameter is...'+str(delay_paramX))
    print('\n\n\n\nNewly generated delayY finer parameter is...'+str(delay_paramY))
    #Writing the finer parameter to a file in SAMPLING_INFO folder
    fopenX              =    open('SAMPLING_INFO/ParameterX_'+str(fil1[:-7])+'000.mbr_'+str(fil2[:-7])+'000.mbr.fparam', 'w+')
    fopenY              =    open('SAMPLING_INFO/ParameterY_'+str(fil1[:-7])+'000.mbr_'+str(fil2[:-7])+'000.mbr.fparam', 'w+')
    fopenX.write(str(delay_paramX))#np.savetxt(str('SAMPLING_INFO/Parameter_'+str(fil1[:-7])+'000.mbr_'+str(fil2[:-7])+'000.mbr.fparam'), delay_pa)
    fopenY.write(str(delay_paramY))
    fopenY.close()
    fopenX.close()
    
    return delay_paramX, delay_paramY

cpdef gen_finer_shift_parameter(file_name, file_name1, comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI):

    try:
        fil1    =       file_name.split('/')[-1]
        fil2    =       file_name1.split('/')[-1]
    except:
        fil1    =       file_name
        fil2    =       file_name1
	
    
    cdef str gen_str	=	'SAMPLING_INFO/ParameterX_'+str(fil1[:-7])+'000.mbr_'+str(fil2[:-7])+'000.mbr.fparam'#'SAMPLING_INFO/Parameter_'+fil1+'_'+fil2+'.fparam'
    cdef str gen_str1   =       'SAMPLING_INFO/ParameterX_'+str(fil2[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.fparam'#'SAMPLING_INFO/Parameter_'+fil2+'_'+fil1+'.fparam'
    cdef  double delay_paramanter = 0
    if(len(glob.glob(gen_str)) == 1 or len(glob.glob(gen_str1)) == 1):
        print('Older parameter generated file found, hence not required to find the parameters again.!')
        try:
            return np.loadtxt(gen_str)
        except:
            return np.loadtxt(gen_str1)
    else:
        print('Older parameter generated file not found, hence required to find the parameters again.!')
        return sub_gen_finer_shift_parameter(fil1, fil2, comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI)


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
    cdef np.ndarray xrms        =       np.zeros(256, dtype=float)
    cdef float SNR2             =       avg
    cdef np.ndarray effi	=       np.zeros(256, dtype=float)
    cdef np.ndarray FLAGS       =       np.ones((256), dtype=int)
    cdef float lowlim, highlim
    spec[0]                     =       0
    
    xmean                       =       np.nanmean(spec, axis=1)
    xrms                        =       np.nanstd(spec-np.nanmean(spec), axis=1)#, xmean, fftsize)
    SNR                         =       xmean/xrms
    highlim                     =       np.nanmean(SNR)+sig*np.nanstd(SNR)
    lowlim                      =       np.nanmean(SNR)-sig*np.nanstd(SNR)
    SNR2                        =       np.sqrt(avg)
    effi                        =       SNR/SNR2
    #effi1                       =       np.zeros((15, 256), dtype=float)
    #np.save('Efficiency_latest', efficiency_x) 
    #By default twn itrations.. 
    for i in range(15):
        print('Itration..'+str(i))
        highlim                                 =       np.nanmean(effi)+sig*np.nanstd(effi)
        lowlim                                  =       np.nanmean(effi)-sig*np.nanstd(effi)
        #effi1[i]                =       effi
        effi[effi>highlim]      =       np.nan
        effi[effi<lowlim]       =       np.nan
    #No rejecting 20% upper and lower and data, due to the bandpass effect..
    FLAGS[np.argwhere(np.isnan(effi))]          =       0
    FLAGS[0:10]=0
    FLAGS[246:]=0
    #np.save('Efficiency_latest', effi1)
    return FLAGS

import sys

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


def super_sync_file(fil):
    '''
    Wrapper for super_sync_genrator
    takes file containing paths to the files to be produce syncronization files and synchronization equation
    '''
    fil_list    =   open(fil, 'r').readlines()
    for i in range(len(fil_list)):
        if(i==0):
            temp    =   'super_sync_generator('+str(fil_list[i][:-1])+','
        else:
            temp    =   temp+str(fil_list[i])+','
        temp    =   temp + ')'
    return 0

def super_sync_generator(*args):
    '''
        super_sync_generator:

        This module is use to generate synchroniation files for voltage and stokes mode data sets..using the GPS and the packet
        numbering system of the header. All the files generated will be stored in the SAMPLING_INFO folder (atleast in this machine)
        this can be changes by changing the absolute path hard coded in the 'sub_gen_files' module.

        Usage:

        super_sync_generator(*args)

        Example:

        super_sync_generator('')
    '''
    cdef int i          =   0
    cdef int le         =   len(args)
    dt                  =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    cdef list slp       =   []
    cdef list slp1      =   []
    #cdef list file_name =   []
    cdef list series    =   []

    for i in range(le):
        
        file_name   =   args[i]
        print(file_name)
        series.append(args[i].split('_')[-1])
        slp.append(sub_sync_generator(file_name, series[i], dt))
        slp1.append(str(slp[i].slope) +'*x+'+ str(slp[i].intercept))
        #file_name.append(args[i])
    #sub_gen_files(file_name, series, comf, slp1)



    return slp1

cdef sub_sync_generator(fil_name, series, dt):
    cdef list   comf1         =      []
    cdef list   slp1          =      []     
    cdef list   series1       =      []
    comf1.append(np.memmap(fil_name,  dtype = dt, mode = 'c'))
    cdef list accblipp      =   []
    cdef list accblipg      =   []

    try:
        #Relative or absolute path
        file_name   =   fil_name.split('/')[-1]
    except:
        #Direct path
        file_name   =   fil_name
    print(file_name)
    if(series == '000.mbr'):
        print('Generating for 000.mbr')
        for i in range(10, len(comf1[0])-1):
            if(comf1[0]['GPS'][i+1] - comf1[0]['GPS'][i] == 1 and  comf1[0]['Packet'][i+1] - comf1[0]['Packet'][i] == 1 ):
                accblipp.append(comf1[0]['Packet'][i+1])
                accblipg.append(comf1[0]['GPS'][i+1])
    else:
        print('Generating for non 000.mbr')
        for i in range(len(comf1[0])-1):
            if(comf1[0]['GPS'][i+1] - comf1[0]['GPS'][i] == 1 and  comf1[0]['Packet'][i+1] - comf1[0]['Packet'][i] == 1 ):
                accblipp.append(comf1[0]['Packet'][i+1])
                accblipg.append(comf1[0]['GPS'][i+1])
   

    slp =   linregress(accblipg, accblipp)
    slp1.append(str(slp.slope) +'*x+'+ str(slp.intercept))
    series1.append(series)
    sub_gen_files(file_name, series1, comf1, slp1)
    
    print(slp)
    return slp 

def sub_gen_files(file_name, series, comf, slp):

    cdef int Mem1_index =       0
    cdef int Mem1_fact  =       0
    cdef int Mem1       =       0
    
    #x_st  =   max(gps_st)
    #x_ed  =   max(gps_ed)
    
    try:
        fil1	=	file_name.split('/')[-1]
    except:
        fil1	=	file_name

    for i in range(len(comf)):
        #gps_st[i]   = comf[i]['GPS'][10]
        #gps_ed[i]   = comf[i]['GPS'][-1]
        #packet_st[i]= comf[i]['Packet'][10]
        #packet_ed[i]= comf[i]['Packet'][-1]
        if(series == '000.mbr'):
            x               =               comf[i]['GPS'][0] + 1 #temporarly keeping it 0, for now, will change after changing call_to_read file..
        else:
            x               =               comf[i]['GPS'][0] + 1
        if(series[i] == '000.mbr'):
            print('Writig in loop')
            #If from 000 series jump 10 packets...due to previous buffering issue..
            fileeq1         =               open(str("SAMPLING_INFO/Info_on_straight_line"+str(file_name.split('/')[-1])), 'w+')
            fileeq1.write('First GPS Value == '+str(comf[i]['GPS'][0])+'\n'+str(comf[i]['Packet'][0])+'\n')

        else:
            #If not 000 series do not junp 10 packets, as its just the continuation of the previous file..
            fileeq1         =               open(str("SAMPLING_INFO/Info_on_straight_line"+str(file_name.split('/')[-1])), 'w+')
            fileeq1.write('First GPS Value == '+str(comf[i]['GPS'][0])+'\n'+str(comf[i]['Packet'][0])+'\n')

        fileeq1.write(slp[i]+'\n')


        #For the question of RFI mask generation, we can as well generate this while we are processing the file..this is to avoid the etra proceeing
        #step we will use..so next lines are commented 
        print('Now looking at the series code..'+str(series))
        
        if(series[0] == '000.mbr'):

            #Generating RFI mask if the series is equal to 000.mbr#
            print('Generating RFI mask if the series is equal to 000.mbr')
            comf_X, comf_Y	                 =  call_to_read_for_RFI(comf[i], 10)    #The RFI resolution is fixed at 6000 packets 
            le                                   =  len(comf_X)/(int(int(256)*2)*int(10))
            RFIzzz                               =  np.ones((256), dtype=int)
            creal1, creal2, creal3, creal4       = internal_loop_without_RFI_single_file.external_loop(comf_X, comf_Y, 10, le, RFIzzz, RFIzzz)


            #Getting optimal parameters..#
            RFI1    =   RFI_Reject(creal3, 3.0, 256, 10)
            RFI2    =   RFI_Reject(creal4, 3.0, 256, 10)
            #Writing RFI masks#
            np.savetxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil1))+'_ascii.rfi', RFI1)
            np.savetxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil1))+'_ascii.rfi', RFI2)
    return 0;

cpdef call_to_read_for_RFI(comf_1, avg):
    cdef np.ndarray tempcomf, comf_X, comf_Y
    tempcomf	=	comf_1['data']
    tempcomf	=	tempcomf.ravel()
    comf_X	=	tempcomf[0::2]
    comf_Y	=	tempcomf[1::2]

    return comf_X, comf_Y

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
        if(inter1 >0):
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
    print(slp1)
    print(slp2)
    slp_1   =    linear_fit(accblipg1, accblipp1)
    slp_2   =    linear_fit(accblipg2, accblipp2)
    
    x   =   max(comf['GPS'][0], comf1['GPS'][0])+1
    

    fillist         =               np.loadtxt('DEPEND/files_naming.txt', dtype=str)
    series_code     =               int(np.where(fillist==series)[-1][0])#fillist.index(series)+1
    print('\n\n\n\n\n') 
    print(series_code)

    
    
    #fil_000         =               file_name[:-7]+'000.mbr'
    #fil1_000        =               file_name1[:-7]+'000.mbr'
    


    fileeq1         =               open(str("SAMPLING_INFO/Info_on_straight_line"+str(file_name.split('/')[-1])), 'w+')
    fileeq2         =               open(str("SAMPLING_INFO/Info_on_straight_line"+str(file_name1.split('/')[-1])), 'w+')
    #RFI_mask1       =               np.loadtxt(str("/data/swan/pavan/SAMPLING_INFO/RFI_mask_X_"+str(fil_000.split('/')[-1])))#+fillist[0])
    #RFI_mask2       =               np.loadtxt(str("/data/swan/pavan/SAMPLING_INFO/RFI_mask_Y_"+str(fil_000.split('/')[-1])))
    #RFI_mask3       =               np.loadtxt(str("/data/swan/pavan/SAMPLING_INFO/RFI_mask_X_"+str(fil1_000.split('/')[-1])))#+fillist[0])
    #RFI_mask4       =               np.loadtxt(str("/data/swan/pavan/SAMPLING_INFO/RFI_mask_Y_"+str(fil1_000.split('/')[-1])))

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
    #if(len1 > len2):#if(str(eval(slp_1[-1])) > str(eval(slp_2[-1])) and str(eval(slp_2[-1])) > 0):
    fileeq1.write(str(slp_2[-1])+'\n')
    fileeq2.write(str(slp_2[-1])+'\n')
        
    #    file_to_write       =   file_name.split('_')
    #    file_to_write_fil   =   file_name.split('_')
    #    file_to_write[-1]   =   fillist[int(series_code)+1]
    #    file_to_write_fil[-1]=  fillist[series_code]
    #    file_to_write       =   '_'.join(file_to_write)
    #    print(file_to_write)
    #    file_to_write_fil   =   '_'.join(file_to_write_fil)

    #    filrem          =               open('/data/swan/pavan/SAMPLING_INFO/'+str(file_to_write_fil.split('/')[-1])+'_rem.data', 'w+')
    #    filrem.write(str(rem_len)+'\n')#str(eval(slp_2[-1])-comf1['Packet'][0]-eval(slp_1[-1])+comf['Packet'][0])+'\n')
    #    filrem.write(str(file_to_write)+'\n')

    #    print(fillist[int(series_code)+1])
    fileeq2.write('N,0\n')
    #fileeq1.write('Y,'+str(file_to_write)+','+str(rem_len)+'\n')
    fileeq1.write('N,0\n')

    #    #try:
    #    #    fil_to_rem          =   (file_to_write.split('/')[-1])[0:4]
    #    #    file_to_write       =   (file_to_write.split('/')[-1])[-7:-4]
    #    #except:
    #    #    fil_to_rem          =   file_to_write[0:4]
    #    #    file_to_write       =   file_to_write[-7:-4]
    #    #filrem.write(str(fil_to_rem))
    #    #print(file_to_write)
    #    #print(fil_to_rem)
    #else:
    #    fileeq1.write(str(slp_1[-1])+'\n')
    #    fileeq2.write(str(slp_1[-1])+'\n')


     #   file_to_write       =   file_name1.split('_')
     #   file_to_write_fil   =   file_name1.split('_')
     #   file_to_write[-1]   =   fillist[int(series_code)+1]
     #   file_to_write_fil[-1]=  fillist[series_code]
     #   file_to_write       =   '_'.join(file_to_write)
     #   print('=============')
     #   print(file_to_write)
     #   print('=============')
     #   file_to_write_fil   =   '_'.join(file_to_write_fil)


     #   filrem          =               open('/data/swan/pavan/SAMPLING_INFO/'+str(file_to_write_fil.split('/')[-1])+'_rem.data', 'w+')
     #   filrem.write(str(rem_len)+'\n')#str(eval(slp_1[-1])-comf['Packet'][0]-eval(slp_2[-1])+comf1['Packet'][0])+'\n')
     #   #filrem.write(str(min(eval(slp_1[-1]), eval(slp_2[-1]))-comf['Packet'][0])+'\n')
     #   #filrem(str())
     #   filrem.write(str(file_to_write)+'\n')

     #   #try:
     #   #    fil_to_rem          =   (file_to_write.split('/')[-1])[0:4]
     #   #    file_to_write       =   (file_to_write.split('/')[-1])[-7:-4]
     #   #except:
     #   #    fil_to_rem          =   file_to_write[0:4]
     #   #    file_to_write       =   file_to_write[-7:-4]
     #   #print(file_to_write)
     #   #print(fil_to_rem)
     #   #filrem.write(str(fil_to_rem))

     #   fileeq1.write('N,0\n')
     #   fileeq2.write('Y,'+str(file_to_write)+','+str(rem_len)+'\n')



    fileeq1.write(str(comf['GPS'][0])+','+str(comf['GPS'][-1]))
    fileeq2.write(str(comf1['GPS'][0])+',' + str(comf1['GPS'][-1]))

    num =   comf['Packet'][int(np.where(comf['GPS']==comf['GPS'][-1])[0][0])]
    num1=   comf1['Packet'][int(np.where(comf1['GPS']==comf1['GPS'][-1])[0][0])]

    
    fileeq1.write('\n'+str(comf['Packet'][-1])+'-'+str(num)+' = '+str(comf['Packet'][-1]-num)+','+ str(comf['Packet'][0]))
    fileeq2.write('\n'+str(comf1['Packet'][-1])+'-'+str(num)+' = '+str(comf1['Packet'][-1]-num1)+','+ str(comf1['Packet'][0]))

    print(str(slp1.slope)+'*x'+str(slp1.intercept))

    #Generating RFI mask#
    #if(series == s_mbr):
    comf_X, comf_Y, comf1_X, comf1_Y, LO = call_to_read(file_name, file_name1, 60, '1', np.zeros((30)))
    le      = len(comf_X)/(int(int(256)*2)*int(avg))
    le1     = len(comf1_X)/(int(int(256)*2)*int(avg))
    le      =   min(le, le1)
    creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 = internal_loop2.external_loop(comf_X, comf_Y, comf1_X, comf1_Y, avg, le, 255)
    print('In RFI_Reject')
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

cpdef read_RFI(file_name, file_name1):
    
    try:
        fil         =   file_name.split('/')[-1]
        fil1        =   file_name1.split('/')[-1]
    except:
        fil         =   file_name
        fil1        =   file_name1
    fil		    =		    fil[:-7]+'000.mbr'
    fil1	    =		    fil1[:-7]+'000.mbr'
    print('Reading prior saved RFI data from..'+str(fil)+' and '+str(fil1))
    RFI_mask1       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil))+'_ascii.rfi')#+fillist[0])
    RFI_mask2       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil))+'_ascii.rfi')
    RFI_mask3       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_X_"+str(fil1))+'_ascii.rfi')#+fillist[0])
    RFI_mask4       =               np.loadtxt(str("SAMPLING_INFO/RFI_mask_Y_"+str(fil1))+'_ascii.rfi')


    return RFI_mask1, RFI_mask2, RFI_mask3, RFI_mask4



cpdef tuple decrypy_file_new_SWAN_without_compensation(file_name, file_name1, ch):


    '''
                ch should be the channel number of first file..

		Takes the read file and sorts the X and Y polarizartion in the file into

		comf_X, comf1_X, comf_Y, comf1_Y.

    '''
    #If available get the earlier data set#
    cdef str series      = file_name[-7:-4]
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])
    #if(series   !=  '000'):
    #    rem     =   open('/data/swan/pavan/SAMPLING_INFO/Info_on_straight_line'+str(file_name.split('/')[-1]), 'r').readlines()#open('/data/swan/pavan/SAMPLING_INFO/'+str(file_name.split('/')[-1])+'_rem.data', 'r').readlines()
    #    print(rem)
    #    rem_l   =   rem[-2].split(',')
    #    rem1     =   open('/data/swan/pavan/SAMPLING_INFO/Info_on_straight_line'+str(file_name1.split('/')[-1]), 'r').readlines()
    #    rem_l1   =   rem[-2].split(',')

    #    if(rem_l[0] =='Y'):
    #        comf_rem             = np.memmap(rem_l[1], dtype = dt, mode='c')
    #        comf_rem             = comf_rem[:int(float(rem_l[2][:-1]))]

    #    elif(rem_l1[0] == 'Y'):
    #        comf_rem             = np.memmap(rem_l1[1], dtype = dt, mode='c')
    #        comf_rem             = comf_rem[:int(float(rem_l[2][:-1]))]

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

    #if(series != '000' and file_name[-39:-35] == rem_l[1][-39:-35] and rem_l[0] == 'Y'):
    #    print('\n\n\n\n\n\n\n\nIn 1\n\n\n\n')
    #    comf    =   np.hstack((comf, comf_rem))
    #elif(series != '000' and file_name1[-39:-35] == rem_l[1][-39:-35] and rem_l1[0] == 'Y'):
    #    print('\n\n\n\n\n\n\n\nIn 2\n\n\n\n')
    #    comf1   =   np.hstack((comf1, comf_rem))
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

cpdef call_to_read(str file_name, str file_name1, avg, sysargv, delay):#, number, number1):
    #Binary file structure#
    dt      =    np.dtype([('header', 'S8'), ('Source', 'S10'), ('Attenuator_1', '>u1'),('Attenuator_2', '>u1'), ('Attenuator_3', '>u1'), ('Attenuator_4', '>u1'), ('LO', '>u2'), ('FPGA', '>u2'), ('GPS', '>u2'), ('Packet', '>u4'), ('data', '>i1', 1024)])

    comf_X, comf_Y, comf1_X, comf1_Y,LO        =       decrypy_file_new_SWAN_without_compensation(file_name, file_name1, avg)
    return comf_X, comf_Y, comf1_X, comf1_Y, LO

cpdef obsdelay(creal):
    phase1      =   np.zeros((len(creal[0])), dtype=float)
    phase_total =   []
    ampmax	=   np.zeros((len(creal[0])), dtype=float)
    fact 	=   16/256.0
    fact        =   1/(fact*256.0)   #10**3 is the resampling factor..with 1000 fact come to be around 0.0625ns
    print('10**3 is the resampling factor..with 1000 fact come to be around 0.0625ns')
    print('factor of multiplication is..'+str(fact))
    for i in range(len(creal[0])):                                                                                               
            z            =      np.hstack((creal[:,i], np.zeros((256*(10**3-1)))))                                                                             
            corr         =      np.roll(np.fft.ifft(z)/np.sqrt(len(z)), len(z)/2)                       
            corr1	 =	corr.real                                
            corr2	 =	corr.imag                                
            ampcorr      =   	abs(corr) 
            phase_total.append(ampcorr)	    
            ampmax[i]	 =	(np.argmax(abs(corr)) - len(z)/2)*fact
            #print(argmax(abs(corr)), len(z)/2, ampmax[i])
            phase1[i]	 =	angle(complex(corr1[argmax(ampcorr)],corr2[argmax(ampcorr)]))    
    ampmax	=	ampmax - min(ampmax) 
    return phase1, phase_total, ampmax

cpdef obsdelay_param(creal):
    #ampmax      =   np.zeros((len(creal[0])), dtype=float)
    #for i in range(len(creal[0])):
    z            =      np.hstack((np.mean(creal, axis=1), np.zeros((256*(10**4-1)))))
    corr         =      np.roll(np.fft.ifft(z)/np.sqrt(len(z)), len(z)/2)
    corr1        =      corr.real
    corr2        =      corr.imag
    ampcorr      =      abs(corr)
    #Fitting a parabolic curve to get accurate results
    px, py       =       parabolic(ampcorr, np.argmax(abs(corr)))
    ampmax	 =       (px - len(z)/2)
    print(px, ampmax)
    #np.savetxt('Current_IFFT', corr)
    print('\n\n\n\nampmax is..'+str(ampmax/5000.0))
    return ampmax



  
cpdef parabolic(f, x):
    """Quadratic interpolation for estimating the true position of an
    inter-sample maximum when nearby samples are known.
   
    f is a vector and x is an index for that vector.
   
    Returns (vx, vy), the coordinates of the vertex of a parabola that goes
    through point x and its two neighbors.
   
    Example:
    Defining a vector f with a local maximum at index 3 (= 6), find local
    maximum if points 2, 3, and 4 actually defined a parabola.
   
    In [3]: f = [2, 3, 1, 6, 4, 2, 3, 1]
   
    In [4]: parabolic(f, argmax(f))
    Out[4]: (3.2142857142857144, 6.1607142857142856)
   
    """
    # Requires real division.  Insert float() somewhere to force it?
    xv = 1/2 * (f[x-1] - f[x+1]) / (f[x-1] - 2 * f[x] + f[x+1]) + x
    yv = f[x] - 1/4 * (f[x-1] - f[x+1]) * (xv - x)
    return (xv, yv)

cpdef parabolic_polyfit(f, x, n):
    """Use the built-in polyfit() function to find the peak of a parabola
    f is a vector and x is an index for that vector.
    n is the number of samples of the curve used to fit the parabola.
    """
    a, b, c = polyfit(arange(x-n//2, x+n//2+1), f[x-n//2:x+n//2+1], 2)
    xv = -0.5 * b/a
    yv = a * xv**2 + b * xv + c
    return (xv, yv)


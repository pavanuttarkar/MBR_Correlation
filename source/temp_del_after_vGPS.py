#!/usr/bin/env python3


##Uses Cython header file##
##Name: _header_Fring_cy##
import _header_Fring_cy
import _header_read_cy
import _header_geometric_cy

#import matplotlib
#matplotlib.use('Agg')

print("-------------")
from numpy import *
from pylab import *
import time
import sys
import math
import os
from datetime import datetime
import internal_loop_RFI
import time
import warnings
import gc
#import pstats, cProfile

#import pyximport


'''
    Series 000 Correlation, maximizing correlation value..
'''

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


def Correlation_FORT(comf_X_1, comf_Y_1, comf1_X_1, comf1_Y_1, avg, le, file_name, file_name1, fftsize, external_num, c9, fil, delay, RFI):
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
        delay
        RFI array(int): RFI Matrix

    Returns:
        band
        time1
        c9
        creal8
        creal9
        creal5
        creal6

    """

    time1          = np.zeros((2, le+1), dtype=complex)
    band           = np.zeros((2, fftsize), dtype=complex)
    timin      =       time.time()
    creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 = internal_loop_RFI.external_loop(comf_X_1, comf_Y_1, comf1_X_1, comf1_Y_1, avg, le, RFI[0], RFI[2], RFI[1], RFI[3])#, delay, 60000)
    print('Time for internal loop : ' +str(time.time()-timin))

    #####Generating New Folders, if not available#######################
    savenow     =   time.time()
    os.system('mkdir CORRELATION')
    st          =   'DELAY_Sample_len'+str(external_num)
    fg  =   str(file_name).split('/')
    fg1 =   fg[len(fg)-1]

    fh  =   str(file_name1).split('/')
    fh1 =   fh[len(fh)-1]
    #Now compensating for phase..
    creal8  =    _header_geometric_cy.phase_compensation(array(creal8), delay[0], LO-140)
    creal9  =    _header_geometric_cy.phase_compensation(array(creal9), delay[0], LO-140)

    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX1)+str(chX2)+'_'+str(fg1)+'_'+str(fh1)+'.txt', creal8[0])
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chY1)+str(chY2)+'_'+str(fg1)+'_'+str(fh1)+'.txt', creal9[0])
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X1X1_'+str(fg1)+'_'+str(fh1)+'.txt', creal1)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_Y2Y2_'+str(fg1)+'_'+str(fh1)+'.txt', creal2)


    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX1)+str(chY1)+str(fg1)+'_'+str(fh1)+'.txt', pollkg1)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX2)+str(chY2)+str(fg1)+'_'+str(fh1)+'.txt', pollkg2)

    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX1)+str(chY2)+str(fg1)+'_'+str(fh1)+'.txt', creal5)
    np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_'+str(chX2)+str(chY1)+str(fg1)+'_'+str(fh1)+'.txt', creal6)



    #####Checking Correlation and saving all stuff######################

    c9[int(external_num)]      =       np.nanmean(abs(band[0]))
    #print('Mean  Y1Y2...'+str(np.nanmean(abs(band[0]))))
    #print('Mean  X1X2...'+str(np.nanmean(abs(band[1]))))
    #np.save('CORRELATION/'+str(fil)+'/Correlation_X1X2_Y1Y2_Band'+str(fg1)+'_'+str(fh1)+'.txt',band)
    #np.save('CORRELATION/'+str(fil)+'/Correlation_X1X2_Y1Y2_Time'+str(fg1)+'_'+str(fh1)+'.txt',time1)

    print('Time for python function = '+ str(time.time()-timin))

    return band, time1, c9, creal8, creal9, creal5, creal6

if(len(sys.argv) < 7):
	print('Usage: Program_name <1 for Synchronization  otherwise 0 > <File_name_1> <File_name_2> <NaN> <No of loops to run, e.g 10 delays etc.but check the internal multiplication!> <RA> <Dec> <Avg packets> <Tile 1> <Tile 2> <FFT size>')
	# sys.exit()

file_name		=       sys.argv[2]	
file_name1		=	sys.argv[3]
#delay                   =       load(sys.argv[4])
avg                     =       int(sys.argv[8])
RA                      =      float(sys.argv[6])
Dec                     =      float(sys.argv[7])
T1                      =       int(sys.argv[-2])
T2                      =       int(sys.argv[-3])

#Extracting date and time tags from file name and generating time jump matrix#
series, sec, minu, hour, dat, month, year	=	_header_Fring_cy.extract_time(file_name)
series1, sec1, minu1, hour1, dat1, month1, year1	=	_header_Fring_cy.extract_time(file_name1)
sec =   max(sec, sec1)




#GPS_Synchronization starting#
Memfactor = _header_Fring_cy.gps_sync(file_name, file_name1, sys.argv[1])
slope     = [['64453.125']]#get_slope(file_name, file_name1)

#Getting RFI matrix#
RFI         =   _header_read_cy.read_RFI(file_name, file_name1)

#Generating Delay#
print(sec, minu, hour, dat, month, year)
print(sec1, minu1, hour1, dat1, month1, year1)
sece    =   int(sec)    +   max(list(map(int, Memfactor[2])))+1#    +   30
mint    =   int(minu)
hourt   =   int(hour)

print('Going into fix_time')
sec1, min1, hour1, day1     =   _header_read_cy.fix_time(int(sece), int(mint), int(hourt), int(dat))

lst1    =   _header_geometric_cy.selflst(int(sece), int(minu), int(hour), int(dat), int(month), int(year))

delay                       =       _header_geometric_cy.phase_compensation_toggle(sec1, min1, hour1, day1, month, year, RA,  Dec, avg, 1, 40, T1, T2)#.Cal_time_onhold_v1(float(sec1), float(min1), float(hour1), float(day1), float(month), float(year), float(RA),  float(Dec), float(avg), 1.0, 31, T1, T2)
delay        =   [zeros((40)), 1, zeros((40))]
##################
##CHECK THIS!!!!##
##################
sec1, min1, hour1, day1     =   _header_read_cy.fix_time(int(sec)+min(list(map(int, Memfactor[3]))), mint, hourt, dat)

lst2    =   _header_geometric_cy.selflst(int(sec1), int(min1), int(hour1), int(day1), int(month), int(year))



timea    = time.time()
comf_X, comf_Y, comf1_X, comf1_Y, LO    =  _header_read_cy.call_to_read(file_name, file_name1, file_name[-39:-35], sys.argv[1], delay[2])
print('Time taken to read files and decrypt...' +str(time.time()-timea))



#Warning about swapping#
print(bcolors.WARNING+'Remember!! keep avg below 60000, or swapping will happen!!'+bcolors.ENDC)

if(sys.argv[1] == str(1)):

    comf_X    =   comf_X[int(round(Memfactor[0])):]#-int(sys.argv[5])/2:]#+int(round(float(sys.argv[5])))*1/2:] #This is to ensure the correlation peak is seen at the mid-point.
    comf1_X   =   comf1_X[int(round(Memfactor[1])):]

    comf_Y    =   comf_Y[int(round(Memfactor[0])):]#-int(sys.argv[5])/2:]#+int(round(float(sys.argv[5])))*1/2:]
    comf1_Y   =   comf1_Y[int(round(Memfactor[1])):]


print('\n\n\nNumber of samples jumped at first..'+str(int(round(float(sys.argv[5])))*1/2-15)+'\n\n\n')

print('\n\n\n\n')
print(bcolors.BOLD+'Sample Jump..'+bcolors.ENDC)
print(bcolors.BOLD+str(int(round(Memfactor[0]))-int(round(float(sys.argv[5])))*1/2)+bcolors.ENDC)
print(bcolors.BOLD+str(int(round(Memfactor[1]))-int(round(float(sys.argv[5])))*1/2)+bcolors.ENDC)
print('\n\n\n\n')

#sys.exit()

#Generatinf end point#
le      = len(comf_X)/(int(int(sys.argv[-1])*2)*int(avg))
le1     = len(comf1_X)/(int(int(sys.argv[-1])*2)*int(avg))
le      =   min(le, le1)

'''
diff_len=   abs(len(comf_X)-len(comf1_X))
if(len(comf_X)>len(comf1_X)):
    temp_name   =   file_name.split('_')[-4]
    try:
        temp_fil    =   file_name.split('/')[-1]
    except:
        temp_fil    =   file_name
    save('/data/swan/pavan/tmp/comf_X', comf_X)
    save('/data/swan/pavan/tmp/comf_Y', comf_Y)
elif(len(comf1_X)>len(comf_X)):
    temp_name   =   file_name.split('_')[-4]
    try:
        temp_fil    =   file_name1.split('/')[-1]
    except:
        temp_fil    =   file_name1
    save('/data/swan/pavan/tmp/comf_X', comf_X, header=file_name1)
    save('/data/swan/pavan/tmp/comf_Y', comf_Y, header=file_name1)
'''
#Building file to save#
tim         =   str('{0:02d}'.format((datetime.now().hour)))+str('{0:02d}'.format((datetime.now().minute)))+str('{0:02d}'.format((datetime.now().second)))
dat         =   str('{0:04d}'.format((datetime.now().year)))+str('{0:02d}'.format((datetime.now().month)))+str('{0:02d}'.format((datetime.now().day)))
corrfile    =   dat+'_'+tim+'.corr'

#Range of samples to shift#
#If input number if 10, then the samples are shifted from -10 to 10.

c9          =       zeros((int(round(float(sys.argv[5])))*2+1))
samp_range  =   linspace(0, int(round(float(sys.argv[5]))), int(round(float(sys.argv[5]))), dtype =int)#linspace(-int(round(float(sys.argv[5]))), int(round(float(sys.argv[5]))), int(round(float(sys.argv[5])))*2+1)
count_c9    =   0
#Always only comf will be shifted and not comf1#

fil1    =   file_name.split('_')
fil2    =   file_name1.split('_')
if(fil1[-5] != fil2[-5]):
    warnings.warn('Both file names seems to be different')
chX1      =   'X'+str((fil1[-5])[-1])
chX2      =   'X'+str((fil2[-5])[-1])

chY1      =   'Y'+str((fil1[-5])[-1])
chY2      =   'Y'+str((fil2[-5])[-1])



#Writing Metadata#
fil                         =   fil1[-4]+'_'+fil1[-3]+'_'+fil1[-3]#dat+'_'+tim
os.system('mkdir CORRELATION/'+fil)
point   =   open('CORRELATION/'+str(fil)+'/Metadata_'+str(file_name.split('/')[-1])+'_'+str(file_name1.split('/')[-1])+'.txt', 'w+')
point.write('#SOURCE Obs_Frequency\tRA\tDec\tN-FFT\tIntegration_Time\tStart_Time_LST\tStop_Time_LST\n')
point.write(str(fil)+'\t'+str(LO)+'\t'+str(RA)+'\n'+str(Dec)+'\t'+str(sys.argv[-1])+'\t'+str(avg)+'\t'+str(lst1)+'\t'+str(lst2))
point.close()



print('Time till loop..' +str(time.time()-timea))
for jjj in samp_range:

    comf_X1  =   comf_X[int(jjj)*1:].copy()
    comf_Y1  =   comf_Y[int(jjj)*1:].copy()

    st                          =       str(fil)+"/DELAY_Sample_len"+str(int(jjj))
    os.system("mkdir CORRELATION/"+str(fil)+"/DELAY_Sample_len"+str(int(jjj)))
    #print('Going into Cython')
    YYb, YYt, c9, creal8, creal9, creal5, creal6        =       Correlation_FORT(comf_X1, comf_Y1, comf1_X, comf1_Y, avg, le,  file_name, file_name1, int(sys.argv[-1]), int(jjj), c9, fil, delay, RFI)
    count_c9    =   count_c9+1

savetxt(corrfile+'_vec',c9, header = 'Correlation co-efficient values for complete file')
gc.collect()

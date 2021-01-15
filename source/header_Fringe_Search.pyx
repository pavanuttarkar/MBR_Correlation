#from matplotlib.pyplot import *
import resource
import signal
import numpy as np
cimport numpy as np
#from numpy import *
#from numpy cimport *
from numpy cimport *
import sys
import os
import binascii
import cmath
import time
from os.path import getsize
import sys
#from astropy.io import ascii
import math
import os
from ctypes import c_int8
from scipy.signal import *
from scipy.fftpack import *
#from pyfftw.interfaces.numpy_fft import fft
#from pyfftw import FFTW
#pyfftw.interfaces.scipy_fftpack.fft(a, threads=multiprocessing.cpu_count())
from fractions import Fraction
import multiprocessing
from pyfftw import FFTW
import pyfftw
#import internal_loop2
from libc.stdlib cimport malloc, free
from scipy.interpolate import *

#from ctypes import c_int, byref, cdll
#addlib = cdll.LoadLibrary('./internal_FORTRAN_ctyp1.so')
nthr   =  1 

#Address location
ADDR_LOC = "./SAMPLING_INFO/"

cdef ndarray rms(a,meana, fftsize):
    """
    Function to Calcualte RMS for Efficiency Calculation.

    Args:
        a array(float): Input Array
        meana (float): Mean of the Input Array
        fftsize (int): FFT Size
    
    Returns:
        - **rmsa** array(complex) RMS of Input
    """

    cdef ndarray rmsa	=	np.zeros(int(fftsize), dtype='complex')
    cdef int i		=	0
    for i in range(int(fftsize)):
            #for j in range(len(a[0])):
                rmsa[i] = np.sqrt(np.mean((a[i] - meana[i])**2))#    rmsa[i] += (a[i][j] - meana[i])**2
            #rmsa[i]=sqrt(rmsa[i]/len(a[0]))
    return rmsa

cdef list RFI_Reject(spec, sig, fftsize, avg):
    '''
    RFI Rejection module,

    Args:      
        spec, sigma deviation, fftsize, avg

    Returns:
        - **RFI_list** array(int)
    '''
    cdef int i 			= 0
    cdef ndarray SNR 		=	np.zeros(256, dtype=float)
    cdef ndarray xmean		=       np.zeros(256, dtype=float)
    cdef ndarray arms		=	np.zeros(256, dtype=float)
    cdef float SNR2			=	avg
    cdef ndarray efficiency_x	=       np.zeros(256, dtype=float)
    cdef list RFI_list          =       []       
    xmean				=	np.mean(spec, axis=0)
    xrms            		=       rms(spec, xmean, fftsize)
    SNR             		=       xmean/xrms
    efficiency_x    		=       SNR/SNR2
    for i in range(256):
        if(efficiency_x[i] > (np.mean(efficiency_x)+int(sig)*np.std(efficiency_x)) or  efficiency_x[i] < (np.mean(efficiency_x)-int(sig)*np.std(efficiency_x))):
            #spec[i]   =   np.zeros((len(spec[0])))
            RFI_list.append(i)
    return RFI_list

def gps_sync(file_name, file_name1, fact):
    cdef list SynFile       =   []
    cdef ndarray ext        =   np.zeros((2), dtype=float)
    cdef list f512          =   []
    cdef list args 	    =   []
    cdef int ctr	    =   0
    cdef int numprocess     =   len(ext)
    cdef list NoPacktoSkip  =       []
    cdef list barebinary    =       []
    cdef list GPScount      =       []
    cdef list PackCount     =       []
    #cdef list data          =       np.zeros((2, 512))
    cdef list Memfactor     =       []
    cdef int i		    =	    0
    cdef list GPS_st        =       []
    if(fact == str(1)):
        print('Synchcronization in progress>>>>>>>>')

        ext[0]  =   (os.system(''.join(['ls -lrth ', ADDR_LOC, 'Info_on_straight_line' + str(file_name.split('/')[-1]) ] )))
        ext[1]  =   (os.system(''.join(['ls -lrth ', ADDR_LOC, 'Info_on_straight_line' + str(file_name1.split('/')[-1]) ] )))
        if(ext[0] == 512):
            f512.append(file_name)
        else:
            SynFile.append(open(ADDR_LOC + 'Info_on_straight_line'+str(file_name.split('/')[-1])))

        if(ext[1] == 512):
            f512.append(file_name1)
        else:
            SynFile.append(open(ADDR_LOC + 'Info_on_straight_line'+str(file_name1.split('/')[-1])))


        for i in range(2):
            GPScount.append(int(((SynFile[i].readline().split(" ")[-1]).split('\n')[0]).split('.')[0]))
            PackCount.append(int(SynFile[i].readline().split('.')[0]))
            NoPacktoSkip.append(SynFile[i].readline())
            junk    =   SynFile[i].readline()
            GPS_st.append(SynFile[i].readline().split(','))
            print('GPS count..file no.' +str(i)+'..'+str(GPScount[i]))
        x           =   max(GPScount)+1
        for i in range(2):
            a           =   math.modf(eval(NoPacktoSkip[i]))
            skipint     =   int(a[1])
            skipfloat   =   float(a[0])	
            Memfactor.append(((skipint+skipfloat) - PackCount[i])*512.0)
        Memfactor.append(GPS_st[0])
        Memfactor.append(GPS_st[1])
        #print('Skiping --'+str(Memfactor*1056))

        return Memfactor
    else:
        Memfactor = [0, 0]
        return Memfactor

def optimize(file_name, file_nam1):
    cdef int opt    =   0
    print('<<<<<<<<<<<Optimization in progress>>>>>>>>>>>>>')
        
    ext1  =   os.system('ls -lrth Optimal_paramaeter_'+str(file_name.split('/')[-1])+'_'+str(file_name1.split('/')[-1]))
    ext2  =   os.system('ls -lrth Optimal_paramaeter_'+str(file_name1.split('/')[-1])+'_'+str(file_name.split('/')[-1]))
    
    if(ext == 512):
            f512.append(file_name)
    else:
            SynFile.append(open(ADDR_LOC + 'Info_on_straight_line'+str(file_name.split('/')[-1])))
    

def get_slope(file_name, file_name1):
    cdef list SynFile       =   []
    cdef ndarray ext        =   np.zeros((2), dtype=float)
    cdef list f512          =   []
    cdef list args 		=   []
    cdef int ctr		=   0
    cdef int numprocess  	=   len(ext)
    cdef list NoPacktoSkip  =       []
    cdef list barebinary    =       []
    cdef list GPScount      =       []
    cdef list PackCount     =       []
    #cdef list data          =       np.zeros((2, 512))
    cdef list Memfactor     =       []
    cdef int i		=	0
    if(sys.argv[1] == str(1)):
        print('Synchcronization in progress>>>>>>>>')

        ext[0]  =   (os.system(''.join(['ls -lrth', ADDR_LOC, 'Info_on_straight_line'+str(file_name.split('/')[-1])] )))
        ext[1]  =   (os.system(''.join(['ls -lrth', ADDR_LOC, 'Info_on_straight_line'+str(file_name1.split('/')[-1])] )))
        if(ext[0] == 512):
            f512.append(file_name)
        else:
            SynFile.append(open(ADDR_LOC + 'Info_on_straight_line'+str(file_name.split('/')[-1])))

        if(ext[1] == 512):
            f512.append(file_name1)
        else:
            SynFile.append(open(ADDR_LOC + 'Info_on_straight_line'+str(file_name1.split('/')[-1])))


        for i in range(2):
            GPScount.append(int(((SynFile[i].readline().split(" ")[-1]).split('\n')[0]).split('.')[0]))
            PackCount.append(int(SynFile[i].readline().split('.')[0]))
            NoPacktoSkip.append(SynFile[i].readline().split(" ")[0])
        Memfactor.append(NoPacktoSkip)
        #print('Skiping --'+str(Memfactor*1056))
        return Memfactor





def extract_time(fil):
    """
    Function for extracting time parameter from file name
    
    Args:
        fil (str): File Name

    Returns:
        - **series** (str)
        - **sec** (str)
        - **minu** (str)
        - **hour** (str)
        - **dat** (str)
        - **month** (str)
        - **year** (str)
    """
    cdef str series                  =       (fil[-7:-4])
    cdef str dat                     =       (fil.split("_")[-3])[6:8]
    cdef str month                   =       (fil.split("_")[-3])[4:6]
    cdef str year                    =       (fil.split("_")[-3])[0:4]
    cdef str hour                    =       (fil.split("_")[-2])[0:2]
    cdef str minu                    =       (fil.split("_")[-2])[2:4]
    cdef str sec                     =       (fil.split("_")[-2])[4:6]
    return series, sec, minu, hour, dat, month, year
cpdef tuple decrypy_file_new_packet_comp(comf, comf1, avg):

    '''

		Takes the read file and sorts the X and Y polarizartion in the file into
		Without packet loss compensation..
		comf_X, comf1_X, comf_Y, comf1_Y.

    '''
    le      = min(len(comf)/avg, (len(comf1)/avg))

    #cdef list tempcomf    =   comf#.copy()
    #cdef list tempcomf1   =   comf1#.copy()
    tempcomf    =   comf
    tempcomf1   =   comf1

    comf        =   []
    comf1       =   []
    comf.append('MBRDSP00000000000000000000000000')
    comf1.append('MBRDSP00000000000000000000000000')
    cdef ndarray packet = np.zeros((2))
    global packet1
    global packet2
    global gps1 
    global gps2
    packet1 = np.zeros((len(tempcomf)))
    packet2 = np.zeros((len(tempcomf1)))
    gps1    = np.zeros((len(tempcomf)))	
    gps2    = np.zeros((len(tempcomf1)))
    #sys.exit()
    for i in range(len(tempcomf)):
        gps1[i]     =  int(binascii.hexlify(tempcomf[i][20:22]),16)
        packet1[i]  =  int(binascii.hexlify(tempcomf[i][22:26]),16)   
        packet[i%2] =   packet1[i]
        #if((packet[i%2]-packet[i%2-1] > 1) and i != 0):
        #    for k in range(int(abs(packet[i%2]-packet[i%2-1]))):
        #        comf.append(str('\x00')*1024)
        #        comf.append('MBRDSP00000000000000000000000000')
        comf.append(tempcomf[i][26:])
        comf.append('MBRDSP00000000000000000000000000')
    packet = np.zeros((2))
    for i in range(len(tempcomf1)):
        gps2[i]     =  int(binascii.hexlify(tempcomf1[i][20:22]),16)
        packet2[i]  =  int(binascii.hexlify(tempcomf1[i][22:26]),16)
        packet[i%2] =   packet2[i]
        #if((packet[i%2]-packet[i%2-1] > 1) and i != 0):
        #    for k in range(int(abs(packet[i%2]-packet[i%2-1]))):
        #        comf1.append(str('\x00')*1024)
        #        comf1.append('MBRDSP00000000000000000000000000')
        comf1.append(tempcomf1[i][26:])
        comf1.append('MBRDSP00000000000000000000000000')

    del tempcomf
    del tempcomf1

    comf = ''.join(comf)
    comf1 = ''.join(comf1)

    comf    =   comf.split('MBRDSP00000000000000000000000000')[1:]
    comf1   =   comf1.split('MBRDSP00000000000000000000000000')[1:]
    le      = min(len(comf)/avg, (len(comf1)/avg))

    comf    = ''.join(comf)
    comf1   = ''.join(comf1)
	
    #cdef ndarray comf_X		=	np.zeros((len(comf)/2), dtype=str)
    #cdef ndarray comf_Y     	=       np.zeros((len(comf)/2), dtype=str)
    #cdef ndarray comf1_X     	=       np.zeros((len(comf)/2), dtype=str)
    #cdef ndarray comf1_Y     	=       np.zeros((len(comf)/2), dtype=str)
	
    comf_X  			=   comf[1::2]
    comf_Y  			=   comf[0::2]
    comf1_X 			=   comf1[1::2]
    comf1_Y 			=   comf1[0::2]
    
    ft      =   np.dtype('>i1')
    comf_X  =   np.frombuffer(comf_X,  dtype=ft, count = len(comf_X))
    comf_Y  =   np.frombuffer(comf_Y,  dtype=ft, count = len(comf_Y))
    comf1_X =   np.frombuffer(comf1_X, dtype=ft, count = len(comf1_X))
    comf1_Y =   np.frombuffer(comf1_Y, dtype=ft, count = len(comf1_Y))

    
    
    return comf_X, comf_Y, comf1_X, comf1_Y

cpdef tuple decrypy_file_new_SWAN(comf, comf1, avg):

    '''

		Takes the read file and sorts the X and Y polarizartion in the file into

		comf_X, comf1_X, comf_Y, comf1_Y.

    '''
    le      = min(len(comf)/avg, (len(comf1)/avg))
    #cdef list tempcomf    =   comf#.copy()
    #cdef list tempcomf1   =   comf1#.copy()
    tempcomf    =   comf
    tempcomf1   =   comf1

    comf        =   []
    comf1       =   []
    comf.append('SWAN0000000000000000000000000000')
    comf1.append('SWAN0000000000000000000000000000')
    cdef ndarray packet = np.zeros((2))
    global packet1
    global packet2
    global gps1 
    global gps2
    packet1 = np.zeros((len(tempcomf)))
    packet2 = np.zeros((len(tempcomf1)))
    gps1    = np.zeros((len(tempcomf)))	
    gps2    = np.zeros((len(tempcomf1)))
    #sys.exit()
    for i in range(len(tempcomf)):
        gps1[i]     =  int(binascii.hexlify(tempcomf[i][22:24]),16)
        packet1[i]  =  int(binascii.hexlify(tempcomf[i][24:28]),16)   
        packet[i%2] =   packet1[i]
        #print(packet1[i])
        if((packet[i%2]-packet[i%2-1] > 1) and i != 0):

            #print('Packet loss deteted..')
            for k in range(int(abs(packet[i%2]-packet[i%2-1]))):
                comf.append(str('\x00')*1024)
                comf.append('SWAN0000000000000000000000000000')
        else:
                comf.append(tempcomf[i][28:])
        comf.append('SWAN0000000000000000000000000000')
    packet = np.zeros((2))
    for i in range(len(tempcomf1)):
        gps2[i]     =  int(binascii.hexlify(tempcomf1[i][22:24]),16)
        packet2[i]  =  int(binascii.hexlify(tempcomf1[i][24:28]),16)
        packet[i%2] =   packet2[i]
        #print(packet2[i])
        if((packet[i%2]-packet[i%2-1] > 1) and i != 0):
            #print('Packet loss detected..')
            for k in range(int(abs(packet[i%2]-packet[i%2-1]))):
                comf1.append(str('\x00')*1024)
                comf1.append('SWAN0000000000000000000000000000')
        else:
            comf1.append(tempcomf1[i][28:])
        comf1.append('SWAN0000000000000000000000000000')

    del tempcomf
    del tempcomf1

    comf = ''.join(comf)
    comf1 = ''.join(comf1)

    comf    =   comf.split('SWAN0000000000000000000000000000')[1:]
    comf1   =   comf1.split('SWAN0000000000000000000000000000')[1:]
    le      = min(len(comf)/avg, (len(comf1)/avg))

    comf    = ''.join(comf)
    comf1   = ''.join(comf1)
	
    #cdef ndarray comf_X		=	np.zeros((len(comf)/2), dtype=str)
    #cdef ndarray comf_Y     	=       np.zeros((len(comf)/2), dtype=str)
    #cdef ndarray comf1_X     	=       np.zeros((len(comf)/2), dtype=str)
    #cdef ndarray comf1_Y     	=       np.zeros((len(comf)/2), dtype=str)

    comf_X  			=   comf[1::2]
    comf_Y  			=   comf[0::2]
    comf1_X 			=   comf1[1::2]
    comf1_Y 			=   comf1[0::2]

    ft      =   np.dtype('>i1')
    comf_X  =   np.frombuffer(comf_X,  dtype=ft, count = len(comf_X))
    comf_Y  =   np.frombuffer(comf_Y,  dtype=ft, count = len(comf_Y))
    comf1_X =   np.frombuffer(comf1_X, dtype=ft, count = len(comf1_X))
    comf1_Y =   np.frombuffer(comf1_Y, dtype=ft, count = len(comf1_Y))


    return comf_X, comf_Y, comf1_X, comf1_Y




cpdef tuple internal_loop_exe(np.ndarray[np.int8_t, ndim=1]  chunk_X_1, np.ndarray[np.int8_t, ndim=1]  chunk_Y_1, np.ndarray[np.int8_t, ndim=1]  chunk1_X_1, np.ndarray[np.int8_t, ndim=1]  chunk1_Y_1, unsigned int avg, unsigned int fftsize, unsigned int wis_count):
        
        if(wis_count == 0):
            wisdom  =   pyfftw.export_wisdom()
            np.save('wisdom', wisdom)
        else:
            wisdom  =   np.load('wisdom.npy')
            pyfftw.import_wisdom(wisdom)

        fX1                                     =pyfftw.empty_aligned(512, dtype='complex128')
        fX2                                     =pyfftw.empty_aligned(512, dtype='complex128')
        fY1                                     =pyfftw.empty_aligned(512, dtype='complex128')
        fY2                                     =pyfftw.empty_aligned(512, dtype='complex128')

        iX1                                     =pyfftw.empty_aligned(512, dtype='complex128')
        iX2                                     =pyfftw.empty_aligned(512, dtype='complex128')
        iY1                                     =pyfftw.empty_aligned(512, dtype='complex128')
        iY2                                     =pyfftw.empty_aligned(512, dtype='complex128')
        
        oX1                                     =pyfftw.FFTW(iX1, fX1,flags=('FFTW_MEASURE',), direction='FFTW_FORWARD',threads=nthr)
        oX2                                     =pyfftw.FFTW(iX2, fX2,flags=('FFTW_MEASURE',), direction='FFTW_FORWARD',threads=nthr)
        oY1                                     =pyfftw.FFTW(iY1, fY1,flags=('FFTW_MEASURE',), direction='FFTW_FORWARD',threads=nthr)
        oY2                                     =pyfftw.FFTW(iY2, fY2,flags=('FFTW_MEASURE',), direction='FFTW_FORWARD',threads=nthr)



        cdef    ndarray[np.complex128_t,ndim=2] fx2y2_1                 	=np.zeros((avg,fftsize), dtype=complex) #2D matrix for eaiser mean or wgh. avg.
        cdef    ndarray[np.complex128_t,ndim=2] fx1y1_1   			=np.zeros((avg,fftsize), dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] x33_1			=np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] x44_1			=np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] x99_1			=np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] x100_1			=np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] fx1x2_1	     	        =np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] fy1y2_1	                =np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray[np.complex128_t,ndim=2] fcross1_1			=np.zeros((avg,fftsize),dtype=complex)
        cdef    ndarray fcross2_1					=np.zeros((avg,fftsize),dtype=complex)
        cdef    unsigned	int     jj
        cdef    ndarray[np.complex128_t,ndim=1] X			#=np.zeros((fftsize))
        cdef    ndarray[np.complex128_t,ndim=1] Y			#=np.zeros((fftsize))
        cdef    ndarray[char,ndim=1] X2	     		#=np.zeros((fftsize))
        cdef    ndarray[char,ndim=1] Y2			#=np.zeros((fftsize))
        cdef	ndarray[char,ndim=1]	Temp1			#=np.empty(fftsize*2, dtype =int )
        cdef    ndarray[char,ndim=1] Temp2                  #=np.empty(fftsize*2, dtype =int )
        cdef    ndarray[char,ndim=1] Temp3                  # =np.empty(fftsize*2, dtype =int )
        cdef    ndarray[char,ndim=1] Temp4                  # =np.empty(fftsize*2, dtype =int )
        cdef    unsigned int loop_len       =   min((len(chunk_X_1)/(fftsize*2)), (len(chunk1_X_1)/(fftsize*2))) 
        for jj in range(loop_len):#min(len(chunk), len(chunk1))):
                ##Calculating relative packet time##

                X1 =   chunk_X_1[jj*fftsize*2:(jj+1)*fftsize*2]#T1temp[1::2]
                Y1 =   chunk_Y_1[jj*fftsize*2:(jj+1)*fftsize*2]#T1temp[0::2]
                X2 =   chunk1_X_1[jj*fftsize*2:(jj+1)*fftsize*2]#T2temp[1::2]
                Y2 =   chunk1_Y_1[jj*fftsize*2:(jj+1)*fftsize*2]#T2temp[0::2]
                #X1      =   T1temp1#convfrombuffer(T1temp1)
                #Y1      =   T1temp2#convfrombuffer(T1temp2)
                #X2      =   T2temp1#convfrombuffer(T2temp1)
                #Y2      =   T2temp2#convfrombuffer(T2temp2)
               
                
                oX1(X1)#funX1   =    FFTW(X1, fX1)#fX1      =   (fft(X1)[0:fftsize])#(/fftsize*2))[0:fftsize]#[0:256]
                #funX1()
                #fX1     =    fX1[0:fftsize]

                oX2(X2)#funX2   =    FFTW(X2, fX2)#fX1      =   (fft(X1)[0:fftsize])#(/fftsize*2))[0:fftsize]#[0:256]
                #funX2()
                #fX2     =    fX2[0:fftsize]

                oY1(Y1)#funY1   =    FFTW(Y1, fY1)#fX1      =   (fft(X1)[0:fftsize])#(/fftsize*2))[0:fftsize]#[0:256]
                #funY1()
                #fY1     =    fY1[0:fftsize]


                oY2(Y2)#funY2   =    FFTW(Y2, fY2)#fX1      =   (fft(X1)[0:fftsize])#(/fftsize*2))[0:fftsize]#[0:256]
                #funY2()
                #fY2     =    fY2[0:fftsize]

                fX1      =   (fX1)[0:int(fftsize)] 
                fY1      =   (fY1)[0:int(fftsize)]
                fX2      =   (fX2)[0:int(fftsize)]  
                fY2      =   (fY2)[0:int(fftsize)]
		

                #fX1[39:51]  =   np.zeros((12))
                #fY1[39:51]  =   np.zeros((12))
                #fX2[39:51]  =   np.zeros((12))
                #fY2[39:51]  =   np.zeros((12))


                
                #fX1[207:219]  =   np.zeros((12))
                #fY1[207:219]  =   np.zeros((12))
                #fX2[207:219]  =   np.zeros((12))
                #fY2[207:219]  =   np.zeros((12))

                fX1[147:157]  =   np.zeros((10))
                fY1[147:157]  =   np.zeros((10))
                fX2[147:157]  =   np.zeros((10))
                fY2[147:157]  =   np.zeros((10))


                #fX1[169:172]  =   np.zeros((03))
                #fY1[169:172]  =   np.zeros((03))
                #fX2[169:172]  =   np.zeros((03))
                #fY2[169:172]  =   np.zeros((03))


                #fX1[200:203]  =   np.zeros((03))
                #fY1[200:203]  =   np.zeros((03))
                #fX2[200:203]  =   np.zeros((03))
                #fY2[200:203]  =   np.zeros((03))



                #Cross_Correlation_Spectrum#
                fX1X2  	=   fX1*np.conjugate(fX2)#sqrt(fX1*np.conjugate(fX2))#/fftsize
                fY1Y2	=   fY1*np.conjugate(fY2)#sqrt(fY1*np.conjugate(fY2))#/fftsize

                #Check Polarization Leakage#
                fX1Y1   =   fX1*np.conjugate(fY1)#sqrt(fX1*np.conjugate(fY1))
                fX2Y2   =   fX2*np.conjugate(fY2)#sqrt(fX2*np.conjugate(fY2))
        
                #Auto Correlation filling#
                x33_1[jj]	    = fX1*np.conjugate(fX1)#fX1abs
                x44_1[jj]	    = fY1*np.conjugate(fY1)#abs(fY1)#fY1abs
                x99_1[jj]	    = fX2*np.conjugate(fX2)#abs(fX2)#fX2abs
                x100_1[jj]    = fY2*np.conjugate(fY2)#abs(fY2)#fY2abs
                #Polarization Leakage Test#
                fx1y1_1[jj]   = fX1Y1
                fx2y2_1[jj]   = fX2Y2
                #Cross Corrlation Spectrum#
                fx1x2_1[jj]	    =  fX1X2
                fy1y2_1[jj]	    =  fY1Y2
                #Cross Polarization
                fcross1_1[jj]             =   np.conjugate(fX1)*fY2#fX1Y2
                fcross2_1[jj]             =   np.conjugate(fX2)*fY1#fX2Y1

        return x33_1, x44_1, x99_1, x100_1, fx1y1_1, fx2y2_1, fx1x2_1, fy1y2_1, fcross1_1, fcross2_1

def upsample3 (x, k, workers = None):
  """
  Like upsample, but uses the multi-threaded interp3
  """

  assert k >= 1, 'k must be equal or greater than 1'

  mn = x.shape
  if len(mn) == 2:
    m = mn[0]
    n = mn[1]
  elif len(mn) == 1:
    m = 1
    n = mn[0]
  else:
    raise ValueError ("x is greater than 2D")

  nn = n * k

  xt = np.linspace (1, n, n)
  xp = np.linspace (1, n, nn)

  return interp3 (xp, xt, x, workers)

def interp (xp, xt, x):
  """
  Interpolate the signal to the new points using a sinc kernel
  
  Args:
    xt: time points x is defined on
    x: input signal column vector or matrix, with a signal in each row
    xp: points to evaluate the new signal on
 
  Returns:
    - **y** the interpolated signal at points xp
  """

  mn = x.shape
  if len(mn) == 2:
    m = mn[0]
    n = mn[1]
  elif len(mn) == 1:
    m = 1
    n = mn[0]
  else:
    raise ValueError ("x is greater than 2D")

  nn = len(xp)

  y = np.zeros((m, nn))

  for (pi, p) in enumerate (xp):
    si = np.sinc (xt - p)
    y[:, pi] = np.sum(si * x)

  return y.squeeze ()


default_workers = 6
def interp3 (xp, xt, x, workers = default_workers):
  """
  Interpolate the signal to the new points using a sinc kernel
  Like interp, but splits the signal into domains and calculates them
  separately using multiple threads.

  Args:
    xt: time points x is defined on
    x: input signal column vector or matrix, with a signal in each row
    xp: points to evaluate the new signal on
    workers: number of threaded workers to use (default: 16)
  
  Returns: 
    - **y** array(float) the interpolated signal at points xp
  """

  mn = x.shape
  if len(mn) == 2:
    m = mn[0]
    n = mn[1]
  elif len(mn) == 1:
    m = 1
    n = mn[0]
  else:
    raise ValueError ("x is greater than 2D")

  nn = len(xp)

  y = np.zeros((m, nn))

  # from upsample
  if workers is None: workers = default_workers

  xxp = np.array_split (xp, workers)

  from concurrent.futures import ThreadPoolExecutor
  import concurrent.futures

  def approx (_xp, strt):
    for (pi, p) in enumerate (_xp):
      si = np.tile (np.sinc (xt - p), (m, 1))
      y[:, strt + pi] = np.sum (si * x)

  jobs = []
  with ThreadPoolExecutor (max_workers = workers) as executor:
    strt = 0
    for w in np.arange (0, workers):
      f = executor.submit (approx, xxp[w], strt)
      strt = strt + len (xxp[w])
      jobs.append (f)


  concurrent.futures.wait (jobs)

  return y.squeeze ()


cpdef tuple change_freq_new_sinc(data, data1, current_slope, required_slope):
    if(current_slope == required_slope):
        print('Same frequency retained!')
        print('Returning same vales..')
        return data, data1
    cdef    float   df          =   0
    cdef    float   df1         =   0
    cdef    float   num         =   0
    cdef    float   d_num       =   0
    #cdef    float   num        =   0
    cdef    int     i           =   0
    cdef    int     j           =   0
    cdef    int     data_len    =   len(data)
    cdef    float   delt        =   round(1/(len(data)*512), 7)
    cdef    int     nn          =   len(data)*required_slope/required_slope
    cdef    ndarray xp          =   np.linspace (1, n, nn)






    df      =   (current_slope - required_slope) #*512   # Number of samples to be added / second
    print(df)
    if(df < 0):
        df  =   df*-1
    
    #df      =   df*len(data)/(current_slope*512)
    dft     =   round(1/df, 7)
    df1     =   len(data)/df

    print('\nCurrent Slope  :  '+str(current_slope))
    print('\nRequired Slope :  '+str(required_slope)+'\n df\t:'+str(df)+'\n')
    print('\nOld data length:  '+str(len(data)))
    cdef    ndarray datpoint    =   np.zeros((int((required_slope-current_slope))))
    cdef    ndarray datpoint1   =   np.zeros((int((required_slope-current_slope))))
    cdef    list    insert_table1=   [] 
    cdef    ledf                =   int(df)




   



    #Sinc function genration#
    for i in range(ledf):
        datpoint[i]   =   sum(sinc(df1*i - xp[int(df1*i)])*data)
        datpoint[i]    =   sum(sinc(df1*i - xp[int(df1*i)])*data1)
        print(datapoint[i])
        insert_table1.append(int(round(df1*(i))))
        #print('Inserting at : '+str(i))
        #datpoint     =   sum(np.sinc((i)-df1*(i))*data[(i-1000):(i+1000)])
        #datpoint1    =   sum(np.sinc((i)-df1*(i))*data1[(i-1000):(i+1000)])
    
    
    data            =   np.insert(data, insert_table1, datpoint)
    data1           =   np.insert(data1,insert_table1, datpoint1)
    
    return data, data1
          
            
cpdef tuple change_freq_new(data, data1, current_slope, required_slope):
    if(current_slope == required_slope):
        print('Same frequency retained!')
        print('Returning same vales..')
        return data, data1
    cdef    float   df          =   0
    cdef    int     df1         =   0
    cdef    float   num         =   0
    cdef    float   d_num       =   0
    #cdef    float   num        =   0
    cdef    int     i           =   0
    cdef    int     j           =   0
    cdef    int     data_len    =   len(data)
    cdef    float   delt        =   round(1/(len(data)*512), 7)


    df      =   (current_slope - required_slope)*512    # Number of samples to be added / second
    print(df)
    if(df < 0):
        df  =   df*-1
    df      =   df*len(data)/(current_slope*512)
    dft     =   round(1/df, 7)
    df1     =   len(data)/int(df)                    # time/samples # current_slope*512 Sampling Frequency
   
    print('\nCurrent Slope  :  '+str(current_slope))  
    print('\nRequired Slope :  '+str(required_slope)+'\n df\t:'+str(df)+'\n')
    print('\nOld data length:  '+str(len(data))) 
    cdef    ndarray datpoint    =   np.zeros((int(df1)))
    cdef    ndarray datpoint1   =   np.zeros((int(df1)))
    cdef    ndarray insert_table=   np.zeros((int(df1)))

    #Sinc function genration#
    for i in range(int(df1)):
        insert_table[i]     =   df1*i   
        for j in range(10):
            datpoint[i]     =   datpoint[i]   +   data[df1*i]*np.sin(np.pi*(df1*i/(current_slope*512)-j)/(1/(current_slope*512)))/(df1*i/(current_slope*512)-j)/(1/(current_slope*512))   
            datpoint1[i]    =   datpoint1[i]    +   data1[df1*i]*np.sin(np.pi*(df1*i/(current_slope*512)-j)/(1/(current_slope*512)))/(df1*i/(current_slope*512)-j)/(1/(current_slope*512))
    data            =   np.insert(data, insert_table, datpoint)
    data1           =   np.insert(data1,insert_table, datpoint1)
    print('\nNew data length : '+str(len(data)))
    ''' 
    for i in range(data_len):
        j   =   delt*1/current_slope
        if(j%dft    ==  0):
            data    =   np.insert(data, i, 0)#sinc[i] =   sin(np.pi data_len*(512*num))
            data1   =   np.insert(data1, i, 0)
    '''
    #for i in range(int(num)):
    #    for j in range(data_len):
    #        new_samps[i]    =   new_samps[i]+          
    return data, data1
    
#creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 

#creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 


cpdef loop_external_exe(np.ndarray[np.int8_t, ndim=1] comf_X, np.ndarray[np.int8_t, ndim=1] comf_Y, np.ndarray[np.int8_t, ndim=1] comf1_X, np.ndarray[np.int8_t, ndim=1] comf1_Y, unsigned int avg, unsigned int le, str file_name, str file_name1, unsigned int fftsize, int external_num, np.ndarray[double, ndim=1] c9, str fil):
                now =   time.time()
                cdef   creal1       = np.zeros((fftsize,le), dtype = float)
                cdef   creal2       = np.zeros((fftsize,le), dtype = float)
                cdef   creal3       = np.zeros((fftsize,le), dtype = float)
                cdef   creal4       = np.zeros((fftsize,le), dtype = float)
                cdef np.ndarray[np.complex128_t,ndim=2]  creal5       = np.zeros((fftsize,le), dtype = complex)
                cdef np.ndarray[np.complex128_t,ndim=2]  creal6       = np.zeros((fftsize,le), dtype = complex)
                cdef np.ndarray[np.complex128_t,ndim=2]  creal8       = np.zeros((fftsize,le), dtype=complex)
                cdef np.ndarray[np.complex128_t,ndim=2]  creal9       = np.zeros((fftsize,le), dtype=complex)
                cdef np.ndarray[np.complex128_t,ndim=2] pollkg1       = np.zeros((fftsize,le), dtype=complex)	
                cdef np.ndarray[np.complex128_t,ndim=2] pollkg2       = np.zeros((fftsize,le), dtype=complex)
                cdef ndarray time1          = np.zeros((2, le), dtype=complex)
                cdef ndarray band           = np.zeros((2, fftsize), dtype=complex)

                cdef    ndarray[np.complex128_t,ndim=2] fx1y1                 =np.zeros((avg,fftsize), dtype=complex)
                cdef    ndarray[np.complex128_t,ndim=2] fx2y2                 =np.zeros((avg,fftsize), dtype=complex)
                cdef    x33                   =np.zeros((avg,fftsize),dtype=np.float64)
                cdef    x44                   =np.zeros((avg,fftsize),dtype=np.float64)
                cdef    x99                   =np.zeros((avg,fftsize),dtype=np.float64)
                cdef    x100                  =np.zeros((avg,fftsize),dtype=np.float64)
                cdef    ndarray[np.complex128_t,ndim=2] fx1x2                 =np.zeros((avg,fftsize),dtype=complex)
                cdef    ndarray[np.complex128_t,ndim=2] fy1y2                 =np.zeros((avg,fftsize),dtype=complex)
                cdef    ndarray[np.complex128_t,ndim=2] fcross1               =np.zeros((avg,fftsize),dtype=complex)
                cdef    ndarray[np.complex128_t,ndim=2] fcross2               =np.zeros((avg,fftsize),dtype=complex)

                cdef int p_c		     =	0
                timin                        = time.time()
                
                for p_c in range(le):#while(p_c<le):# and m <= end):

                        chunk_X         = comf_X[avg*p_c*fftsize*2:avg*(p_c+1)*fftsize*2]
                        chunk1_X        = comf1_X[avg*p_c*fftsize*2:avg*(p_c+1)*fftsize*2]
                        chunk_Y         = comf_Y[avg*p_c*fftsize*2:avg*(p_c+1)*fftsize*2]
                        chunk1_Y        = comf1_Y[avg*p_c*fftsize*2:avg*(p_c+1)*fftsize*2]



                        
                        #len_1       =       len(chunk_X)-1
                        #len_2       =       len(chunk1_X)-1
                        fftsiz      =       fftsize-1
                        avg1        =       avg-1
                        #print(chunk_X[0:], chunk_Y[0:], chunk1_X[0:], chunk1_Y[0:], avg-1, fftsize-1)
                        #print(p_c, le)
                        #print('len chunk.. '+str(len(chunk_X)))
                        #print('len chunk1.. '+str(len(chunk1_X)))
                        timg        =       time.time()
                        x33, x44, x99, x100, fx1y1, fx2y2, fx1x2, fy1y2, fcross1, fcross2  =internal_loop2.internal_loop(chunk_X[0:], chunk_Y[0:], chunk1_X[0:], chunk1_Y[0:], avg-1, fftsize-1)
                        #RFI_list1   =  RFI_Reject(x33, 1, fftsize, avg)#RFI_Reject_new(x33, 3, 256, avg) 
                        #RFI_list2   =  RFI_Reject(x44, 1, fftsize, avg)
                        #RFI_list3   =  RFI_Reject(x99, 1, fftsize, avg)
                        #RFI_list4   =  RFI_Reject(x100, 1, fftsize, avg)
                        #RFI_listt   =   max(len(RFI_list1), len(RFI_list2), len(RFI_list3), len(RFI_list4))
                        #for i in range(RFI_listt):
                        #    try:
                        #        if(i==RFI_list1[i] or i == RFI_list3[i]):
                        #            fx1x2[i]    =   np.zeros((len(fx1x2[0])))
                        #        if(i==RFI_list2[i] or i == RFI_list4[i]):
                        #            fy1y2[i]    =   np.zeros((len(fy1y2[0])))
                        #
                        #    except:
                        #        pass
                        
                        fy1y2[:,15]   =   np.zeros((len(fy1y2))) 
                        fy1y2            = np.transpose(fy1y2)
                        fy1y2[39:49]     = np.zeros((49-39, len(fy1y2[0])))   
                        fy1y2[103:106]   = np.zeros((106-103, len(fy1y2[0])))
                        fy1y2[77:80]     = np.zeros((80-77, len(fy1y2[0])))
                        fy1y2[137:158]   = np.zeros((158-137, len(fy1y2[0])))
                        fy1y2[170:176]   = np.zeros((176-170, len(fy1y2[0])))#, fy1y2[197:208]
                        fy1y2            = np.transpose(fy1y2)      
                            
                        #x33, x44, x99, x100, fx1y1, fx2y2, fx1x2, fy1y2, fcross1, fcross2  =internal_loop2.internal_loop(chunk_X[0:].astype('float64'), chunk_Y[0:].astype('float64'), chunk1_X[0:].astype('float64'), chunk1_Y[0:].astype('float64'), avg-1, fftsize-1)

                        #x33, x44, x99, x100, fx1y1, fx2y2, fx1x2, fy1y2, fcross1, fcross2  =internal_loop_exe(chunk_X[0:], chunk_Y[0:], chunk1_X[0:], chunk1_Y[0:], avg, fftsize, p_c)
                        #chunk_X_1, chunk_Y_1, chunk1_X_1, chunk1_Y_1, avg, fftsize, wis_count, len_1, len_2, x33_1,x44_1, x99_1, x100_1, fx1y1_1, fx2y2_1, fx1x2_1, fy1y2_1,fcross1_1, fcross2_1 
                        #print('Time for internal loop --'+str(timg-time.time()))
                        #fx1x2, fy1y1        =   RFI_Reject_new(fx1x2,fy1y2, 3, 256, avg)
			#timenow         =  record_start + p_c*avg*1024/33e6+max(Memfactor)*1024/33e6
                        creal1[:,p_c]     = np.mean(abs(x33), axis   =0)#/(avg*fftsize)
                        creal2[:,p_c]     = np.mean(abs(x44), axis   =0)#/(avg*fftsize)
                        creal3[:,p_c]     = np.mean(abs(x99), axis   =0)#/(avg*fftsize)   
                        creal4[:,p_c]     = np.mean(abs(x100), axis  =0)#/(avg*fftsize)
                        pollkg1[:,p_c]    = np.mean(fx1y1, axis =0)#/(avg*fftsize)
                        pollkg2[:,p_c]    = np.mean(fx2y2, axis =0)#/(avg*fftsize)
                        creal8[:,p_c]     = np.mean(fx1x2, axis =0)#/(avg*fftsize)
                        creal9[:,p_c]     = np.mean(fy1y2, axis =0)#/(avg*fftsize)
                        creal5[:,p_c]     = np.mean(fcross1, axis = 0) 
                        creal6[:,p_c]     = np.mean(fcross2, axis = 0)
               
                        #p_c		+=	1
                        #m = m+1
                        #del chunk_X
                        #del chunk_Y
                        #del chunk1_X
                        #del chunk1_Y
                #    bar.next()
                #bar.finish()
                
                #creal1, creal2, creal3, creal4, creal5, creal6, creal8, creal9, pollkg1, pollkg2 = internal_loop2.external_loop (comf_X, comf_Y, comf1_X, comf1_Y, avg, le-1, 256-1) 
                print('Time for internal loop : ' +str(time.time() - timin))

                creal1[0]   =   np.zeros((len(creal1[0])))
                creal2[0]   =   np.zeros((len(creal1[0])))
                creal3[0]   =   np.zeros((len(creal1[0])))
                creal4[0]   =   np.zeros((len(creal1[0])))
                creal8[0]   =   np.zeros((len(creal1[0])))
                creal9[0]   =   np.zeros((len(creal1[0])))


		#####Generating New Folders, if not available#######################
                savenow     =   time.time()
                os.system('mkdir CORRELATION')
                st          =   'DELAY_Sample_len'+str(external_num)
                #os.system('mkdir '+str(st))
                #exec(st1)
                #os.system("mkdir "+str(st))
                #################For Bound check issue###############################
                fg  =   str(file_name).split('/')
                fg1 =   fg[len(fg)-1]
                
                fh  =   str(file_name1).split('/')
                fh1 =   fh[len(fh)-1]
                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X1X2_'+str(fg1)+'_'+str(fh1)+'.txt', creal8)
                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_Y1Y2_'+str(fg1)+'_'+str(fh1)+'.txt', creal9)
                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X1X1_'+str(fg1)+'_'+str(fh1)+'.txt', creal1)
                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_Y2Y2_'+str(fg1)+'_'+str(fh1)+'.txt', creal2)


                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X2Y2_'+str(fg1)+'_'+str(fh1)+'.txt', pollkg1)
                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X1Y1_'+str(fg1)+'_'+str(fh1)+'.txt', pollkg2)

                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X1Y2_'+str(fg1)+'_'+str(fh1)+'.txt', creal5)
                np.save('CORRELATION/'+str(fil)+'/'+str(st)+'/Correlation_X2Y1_'+str(fg1)+'_'+str(fh1)+'.txt', creal6)
                #####Checking Correlation and saving all stuff######################
                deno240     =   np.sqrt(np.mean(creal2, axis = 0)*np.mean(creal4, axis = 0))
                deno241     =   np.sqrt(np.mean(creal2, axis = 1)*np.mean(creal4, axis = 1))

                time1[0]    =   np.mean(creal9[20:230], axis =0)/deno240
                band[0]     =   np.mean(creal9, axis =1)/deno241

                deno80      =   np.sqrt(np.mean(creal1, axis = 0)*np.mean(creal3, axis = 0))
                deno81      =   np.sqrt(np.mean(creal1, axis = 1)*np.mean(creal3, axis = 1))

                time1[1]    =   np.mean(creal8[30:230], axis =0)/deno80
                band[1]     =   np.mean(creal8, axis =1)/deno81


                c9[int(external_num)]      =       np.nanmean(abs(band[0]))
                print('Mean  '+str(np.nanmean(abs(band[0]))))

                np.save('CORRELATION/'+str(fil)+'/Correlation_X1X2_Y1Y2_Band'+str(fg1)+'_'+str(fh1)+'.txt',band)
                np.save('CORRELATION/'+str(fil)+'/Correlation_X1X2_Y1Y2_Time'+str(fg1)+'_'+str(fh1)+'.txt',time1)

                
                print('Time taken for    saving= '+ str(time.time()-savenow))
                print('Time for external loops = '+ str(time.time()-now))
                
                return band, time1, c9, creal8, creal9


cdef obsdelay(creal9):
    """
    Function to calculate the observed delay
    """
    #cdef    list    phase1      =   []
    #cdef    list    z           =   []
    cdef    ndarray phase_total =   np.zeros((len(creal9), len(creal9[0])))
    cdef    ndarray phase1      =   np.zeros((len(creal9[0])))
    cdef    int i               =   0

    for i in range(len(creal9[0])):
            #z                           =   []
            z                           =   np.append(creal9[:,i], np.zeros((256*9)))#append(creal9[:,p_c], zeros((256*9))) #append(average(xy, axis = 0), zeros((256*9)))
            corr                        = ifft(z)
            corr                        = corr/np.sqrt(len(corr))#[corr[ii]/sqrt(len(corr)) for ii in range(len(corr))]
            corr1                       = corr.real#[j.real for j in corr]
            corr2                       = corr.imag#[j.imag for j in corr]
            ampcorr                     = abs(corr)#roll(abs(corr), len(corr)/2)#[sqrt(corr1[k]**2 + corr2[k]**2) for k in range(2560)]
            ampcorr1                    = np.roll(ampcorr, len(ampcorr)/2)
            phase_total[:,i]            = ampcorr1[1152:1408]# for m in range(1152,1408)]
            phase1[i]                   = np.angle(complex(corr1[np.argmax(ampcorr)],corr2[np.argmax(ampcorr)]))*180/np.pi
    i = 0
    while i < len(phase1)-1:
        if ((phase1[i+1] - phase1[i]) > 180):
                phase1[i+1] = phase1[i+1] - 360
        elif ((phase1[i+1] - phase1[i]) < -180):
                phase1[i+1] = phase1[i+1] + 360
        else:
                i = i + 1
    #phase1 = np.array(phase1)
    phase1 = phase1 - min(phase1)
    return phase1, phase_total

def time_exceeded(signo, frame): 
    print("Time's up !") 
    raise SystemExit(1) 
  
def set_max_runtime(seconds): 
    # setting up the resource limit 
    soft, hard = resource.getrlimit(resource.RLIMIT_CPU) 
    resource.setrlimit(resource.RLIMIT_CPU, (seconds, hard)) 
    signal.signal(signal.SIGXCPU, time_exceeded) 


def limit_memory(maxsize): 
    soft, hard = resource.getrlimit(resource.RLIMIT_AS) 
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard)) 
 

cpdef change_freq_builtin(data, data1, current_slope, required_slope):
    '''
    Uses scipy interpolation technique, of interp1d

    Args: 
        data: Data
        data1: Data1
        currentslope: Current Slope
        required_slope: Required Slope
    '''

    if(required_slope < current_slope):
        print('Running Decimation..')
        data,   data1           =   change_freq_new_deci(data, data1, current_slope, required_slope)
        return data, data1

    if(current_slope == required_slope):
        print('Same frequency retained!')
        print('Returning same vales..')
        return data, data1
    cdef    int     data_len    =   len(data)
    cdef    float   delt        =   round(1/(len(data)*512), 7)
    cdef    int     nn          =   (required_slope-current_slope)#len(data)/33000000*(required_slope-current_slope)

    cdef    ndarray     x       =   np.linspace(0, len(data), num = len(data))
    cdef    ndarray     xnew    =   np.linspace(0, len(data), num= len(data)+nn)

    f                           =   interp1d(x, data)
    f1                          =   interp1d(x, data1)


    return  f(xnew), f1(xnew)

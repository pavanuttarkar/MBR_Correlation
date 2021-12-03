cimport numpy as np
import numpy as np
import time
import mmap
import glob
from scipy.stats import linregress
from sympy import S, symbols
from os.path import getsize
from . import single_file_combined_RFInoRFI as internal_loop_without_RFI_single_file
from . import internal_loop_RFI

'''
0 1 2 3 4 5
Y X Y X Y X
'''

cpdef sub_gen_finer_shift_parameter(fil1, fil2, comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI):
    '''
       Module to generate finer delay, after the inital course delay estimation using the embedded
       timming information. This is done by simple cross-correlation and estimation of the delay.

       Parameters
       ----------
       fil1    : string
           1st file name as string
       fil2    : string 
           2nd file name as string
       comf_X  : 1-D numpy array
           array of 8 bit signed integers, of X-Polarization 1st file
       comf_Y  : 1-D numpy array
           array of 8 bit signed integers, of Y-Polarization 1st file
       comf1_X : 1-D numpy array
           array of 8 bit signed integers, of X-Polarization 2nd file
       comf1_Y : 1-D numpy array
           array of 8 bit signed integers, of Y-Polarization 2nd file
       avg     : int
           average value, supplied by the user
       le      : int
           calculated length of the spectrum, in number of columns
       RFI     : 1-D numpy array
           RFI mask array, calculated earlier

       Returns
       -------
       delay-X : unsigned int
           delay value in samples, for X-Polarization, one sample is 30.30ns
       delay-Y : unsigned int
           delay value in samples, for Y-Polarization, one sample is 30.30ns
    '''
    cdef double delay_param = 0#delay[-1]
    print('In sub_gen_finer_shift_parameter')
    #Calculating the dynamic spectrum
    creal1, creal2, creal3, creal4, creal8, creal9  = internal_loop_RFI.external_loop(comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI[0], RFI[2], RFI[1], RFI[3])
    delayX              =    obsdelay_param(creal8)
    delayY              =    obsdelay_param(creal9)
    delay_paramX        =    delayX/5000.0# The second zeros here is insignificant, added due to wierd Cython only error woth savetxt
    delay_paramY        =    delayY/5000.0
    print('Newly generated delayX finer parameter is...'+str(delay_paramX))
    print('Newly generated delayY finer parameter is...'+str(delay_paramY))
    #Writing the finer parameter to a file in SAMPLING_INFO folder
    fopenX              =    open('SAMPLING_INFO/ParameterX_'+str(fil1[:-7])+'000.mbr_'+str(fil2[:-7])+'000.mbr.fparam', 'w+')
    fopenY              =    open('SAMPLING_INFO/ParameterY_'+str(fil1[:-7])+'000.mbr_'+str(fil2[:-7])+'000.mbr.fparam', 'w+')
    fopenX.write(str(delay_paramX))
    fopenY.write(str(delay_paramY))
    fopenY.close()
    fopenX.close()
    
    return delay_paramX, delay_paramY


cpdef gen_finer_shift_parameter(file_name, file_name1, comf_X, comf_Y, comf1_X, comf1_Y, avg, le, RFI):

    '''
       Wrapper module to generate finer delay, after the inital course delay estimation using the embedded
       timming information. This is done by simple cross-correlation and estimation of the delay, redirects
       if previous save info notavailable to sub_gen_finer_shift_parameter.

       Parameters
       ----------
       fil1    : string
           1st file name as string
       fil2    : string 
           2nd file name as string
       comf_X  : 1-D numpy array
           array of 8 bit signed integers, of X-Polarization 1st file
       comf_Y  : 1-D numpy array
           array of 8 bit signed integers, of Y-Polarization 1st file
       comf1_X : 1-D numpy array
           array of 8 bit signed integers, of X-Polarization 2nd file
       comf1_Y : 1-D numpy array
           array of 8 bit signed integers, of Y-Polarization 2nd file
       avg     : int
           average value, supplied by the user
       le      : int
           calculated length of the spectrum, in number of columns
       RFI     : 1-D numpy array
           RFI mask array, calculated earlier

       Returns
       -------
       delay-X : unsigned int
           delay value in samples, for X-Polarization, one sample is 30.30ns
       delay-Y : unsigned int
           delay value in samples, for Y-Polarization, one sample is 30.30ns
       returned as 1-D numpy array

    '''


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


cpdef np.ndarray RFI_Reject(spec, sig, fftsize, avg):    
    '''
       Module to calculate RFI mask for the given auto-correlation
       spectrum. 
       
       Parameters
       ----------
       spec    : 2-D numpy array
           auto-correlation spectrum, having time in x-axis.
       sig     : int
           threashold level to be considered for calculating the RFI mask.
       ffisize : int
           size of the fft taken, i.e if 512 point, then supply 512/2
       avg     : int
           
           average value, supplied by the user

       Returns
       -------
       RFI FLAG: 1-D numpy array 
           Calulated RFI flag corresponding to each channel of the spectrum, 
           if 0, then, RFI affected
           if 1, then, not affected by RFI 
 
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
    return FLAGS

import sys


def super_sync_file(fil):
    '''
       Wrapper for super_sync_genrator
       takes file containing paths to the files to be produce syncronization files and synchronization equation


       Parameters
       ----------
       file : path as string 
           set of files, to calculate synchronization information.

       Returns
       -------
       ----    : ---- 

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
       Wrapper module is use to generate synchroniation files for voltage and 
       stokes mode data sets..using the GPS and the packet numbering system
       of the header. All the files generated will be stored in the SAMPLING_INFO
       folder (atleast in this machine) this can be changes by changing the
       absolute path hard coded in the 'sub_gen_files' module.
       
       Parameters
       ----------
       kwargs  : string
           set of files, to calculate synchronization information.

       Returns
       -------
       ----    : ---- 

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
    '''
        Lower level module to calculate the syncronization solution to the
        given file name, this checks for the GPS and packet transitions and
        fits a straight line to the curve.

       Parameters
       ----------
       filename  : string
           filename, with path is not local, as string
       series    : string
           series for calculating the RFI mask.
       datatype: numpy custom datatype
           datatype used for decryptying the mbr file

       Returns
       -------
       line param: linregress solution
           array of line parameters 

    '''

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
    '''
        Module to calculate the synchronization information using the embedded GPS information.
        This is mainly used in the Intterferometric mode of operation, this is a lower level
        module, generaly masked to the end user. 
    

        Parameters
        ----------
        filname  : string
            name of the file as with it's path, if not
            local.
        series   : string
            series of the file supplied
        comf     : memmap array
            memory mapped array of the file, to calulate
            the RFI mask, if required.
        slp      : float
            output slope from the linregress function, from
            sub_sync_generator upstream function.
 
        Returns
        -------
        ----    : ---- 

    '''
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
    '''
       Module to decrypt file for RFI measurments
       
       Parameters
       ----------
       comf  : memmap array
           memory mapped array of the file supplied
           by the user.
       avg   : int 
           average value, as supplied by the user.
       Returns
       -------
       comf_X: 1-D numpy array
           array of X-Polarization voltage values.
       comf_Y: 1-D numpy array
           array of Y-Polarization voltage values.
    '''

    cdef np.ndarray tempcomf, comf_X, comf_Y
    tempcomf	=	comf_1['data']
    tempcomf	=	tempcomf.ravel()
    comf_X	=	tempcomf[0::2]
    comf_Y	=	tempcomf[1::2]

    return comf_X, comf_Y


cpdef read_RFI(file_name, file_name1):
    '''
       Module to read calculated RFI flags from SAMPLING_INFO folder
       
       Parameters
       ----------
       fil1    : string
           1st file name as string
       fil2    : string 
           2nd file name as string
       
       Returns
       -------
       X-Pol RFI: 1-D numpy array
           RFI mask array, 1st file, X-Polarization, calculated earlier.
       Y-Pol RFI: 1-D numpy array
           RFI mask array, 1st file, Y-Polarization, calculated earlier
       X-Pol RFI: 1-D numpy array
           RFI mask array, 2nd file, X-Polarization, calculated earlier
       Y-Pol RFI: 1-D numpy array
           RFI mask array, 2nd file, Y-Polarization, calculated earlier

    ''' 
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


cpdef obsdelay(creal):
    '''
       Module to calculate the delay in the correlation spectrum,
       by Hilbert Trabsform.

       Parameters
       ----------
       spec        : 2-D numpy array
           cross correlation 2-D matrix, with time
           in horizontal axis.
       
       Returns
       -------
       phase1      : 1-D numpy array
           delay phase, corresponding to individual
           column, calulated by Hilbert Transform.
       phase_total : 2-D numpy array
           2-D array of the phase delay
       ampmax      : 1-D numpy array
            delay calulated using the relative position
            of the ifft peak. 
                
    '''

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
            phase1[i]	 =	np.angle(complex(corr1[np.argmax(ampcorr)],corr2[np.argmax(ampcorr)]))    
    ampmax	=	ampmax - min(ampmax) 
    return phase1, phase_total, ampmax


cpdef obsdelay_param(creal):
    '''
       Module to calculate the finer delay, using the correlation spectrum,
       by Hilbert Transform.

       Parameters
       ----------
       spec    : 2-D numpy array
           cross correlation 2-D matrix, with time
           in horizontal axis.
       
       Returns
       -------
       ampmax  : float
           calculated relative delay
                
    '''

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
    a, b, c = np.polyfit(np.arange(x-n//2, x+n//2+1), f[x-n//2:x+n//2+1], 2)
    xv = -0.5 * b/a
    yv = a * xv**2 + b * xv + c
    return (xv, yv)

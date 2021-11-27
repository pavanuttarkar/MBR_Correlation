import numpy as np
cimport numpy as np
import _header_geometric_cy
import LatLong
import _header_Fring_cy
import _header_gps_cy
from scipy.stats import linregress
import glob
import os

cpdef genphase_ason_OCT23_2021(spec, delay, f):
    cdef np.ndarray pha     =    spec.copy()
    cdef int i, j
    #Changed to positive sign -- JUN20 2021#
    for i in range(256):
        for j in range(len(spec[0])):
            pha[i][j] = (np.exp(complex(0, -2*np.pi*f[i]*delay[j])))
    return pha


cpdef genphase_phase(spec, delay, f):
    pha      =    np.zeros((256, len(spec[0])), dtype=complex)
    for i in range(256):
        for j in range(len(pha[0])):
            pha[i][j] = (np.exp(complex(0, -2*np.pi*f[i]*delay[j])))
    return pha

cpdef genphase_phase_ason_OCT23_2021(spec, delay, f):
    cdef np.ndarray pha     =    spec.copy()
    cdef int i, j
    #Changed to positive sign -- JUN20 2021#
    for i in range(256):
        for j in range(len(spec[0])):
            pha[i][j] = (np.exp(complex(0, -1*delay[j])))
    return pha


cpdef Equ2local(float RA, float Dec, float phi, float lon,
        float second, float minute, float hour, float day, float month, float year):
    '''
        Module to covert Equatorial coordinates to Local Coordinates 
        Equ2local(RA, Dec, phi, lon, second, minute, hour, day, month, year) 
    '''

    cdef double secu        =   second
    cdef double mint        =   minute
    cdef double hourt       =   hour
    cdef double dayt        =   day
    cdef double sindec      =   np.sin(Dec*np.pi/180.0)
    cdef double cosdec      =   np.cos(Dec*np.pi/180.0)
    cdef double tandec      =   np.tan(Dec*np.pi/180.0)
    cdef double d2r         =   np.pi/180

    #secu, mint, hourt, dayt =   correct_time(second, minute, hour, day, month, year)
    #print('Corrected Time..secu, mint, hourt, dayt')
    #print(secu, mint, hourt, dayt)
    lst                     =                   _header_geometric_cy.selflst(second, minute, hour, day, month, year)
    ha                      =                   (lst-RA)#*15.0*np.pi/180
    if(ha<0):
        ha = ha+24.0
    if(ha>12):
        ha = ha-24.0

    ha  =   ha*15*d2r

    coslat               = np.cos(13.6111*np.pi/180)
    sinlat               = np.sin(13.6111*np.pi/180)#77.451944444*np.pi/180)
   
    alt                     = np.arcsin(sindec*sinlat+cosdec*coslat*np.cos(ha))
    #if(alt<0.0):
    #    alt =   alt+np.pi/2
    #az_num                  = np.sin(ha)#sindec-np.sin(alt)*sinlat
    #az_den                  = np.cos(ha)*sinlat - tandec*coslat
    #az_den                  = np.cos(alt)*coslat
    #az                      =   (np.arctan2(az_den, az_num))
    #az                      =   np.arcsin((-1*np.sin(ha)*cosdec)/np.cos(alt))
    #if(ha< 0.0):#RA < LatLong.Lat_local and az < 0.0):
    az                      =   np.arccos((sindec - sinlat * np.sin(alt))  / (coslat * np.cos(alt)))
    if(ha>0.0): 
        az                      =    np.pi-np.arcsin((-1*np.sin(ha)*cosdec)/np.cos(alt))
        #az                      =   np.arcsin((-1*np.sin(ha)*cosdec)/np.cos(alt)) - 2*np.pi
    #az                          =   az-2*np.pi 
    #az                      =   np.arccos((sindec - sinlat * np.sin(alt))  / (coslat * np.cos(alt)))
    #if(az < 0):
    #    az               = 2*np.pi+az#*180/np.pi

    return az, alt, ha
cpdef  geometric_model_using_setdelay(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, np.ndarray time_array, unsigned int T1, unsigned int T2):
#(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2):

    '''
    float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2
    '''
    cdef double c_per_mus   = 299.792458
    cdef int    i       =       0
    print (year, month, day, hour, minu, sec)
    cdef np.ndarray X_param     =      np.zeros(len(time_array), dtype = float)
    cdef np.ndarray Y_param     =      np.zeros(len(time_array), dtype = float)
    cdef np.ndarray Z_param     =      np.zeros(len(time_array), dtype = float)
    cdef np.ndarray W_geometric =      np.zeros(len(time_array), dtype = float)

    cdef float alt              =       0
    cdef float az               =       0
    cdef list az1               =       []
    cdef list alt1              =       []
    coslat               = np.cos(13.6111*np.pi/180)
    sinlat               = np.sin(13.6111*np.pi/180)#77.451944444*np.pi/180)
    
    #lat=13.6112*u.deg, lon=77.5170*u.deg
    #Gettine ECEF Coordinates of Tiles#
    ecef    =   np.loadtxt('ECEF_from_header_geometric.txt')#np.loadtxt('ENU_v6.txt')
    x_loc   =   ecef[:,0][T1] - ecef[:,0][T2]#np.loadtxt('ECEF_x')
    y_loc   =   ecef[:,1][T1] - ecef[:,1][T2]#np.loadtxt('ECEF_y')
    z_loc   =   ecef[:,2][T1] - ecef[:,2][T2]#np.loadtxt('ECEF_z')
    print(x_loc, y_loc, z_loc)
    
    #Getting LatLong#
    cdef double secu        =   sec
    cdef double mint        =   minu
    cdef double hourt       =   hour
    cdef double dayt        =   day
    cdef double montht      =   montht
    cdef double sindec      =   np.sin(dec*np.pi/180)
    cdef double cosdec      =   np.cos(dec*np.pi/180)
    cdef double d2r         =   np.pi/180
    cdef int dayflag        =   0
    ha1 =   []
    for i in range(len(time_array)):
        hourt_temp   =   int(time_array[i])
        mint_temp    =   (time_array[i]%1)*60
        secu_temp    =   (mint_temp%1)*60

        hourt       =   int(hourt_temp)
        mint        =   int(mint_temp)
        secu        =   secu_temp


        az, alt, ha = Equ2local(RA, dec, LatLong.Lat_local, LatLong.Long_local, secu, mint, hourt, day, month, year)
        #print(secu, mint, hourt, dayt, montht, year, az, alt, i)
        alt1.append(alt)
        az1.append(az)
        ha1.append(ha)
        za                   = np.pi/2-alt
        X_param[i]           = x_loc*np.cos(az1[i])*np.cos(alt1[i])/c_per_mus#*np.cos(dec*np.pi/180)*np.cos(ha*np.pi/180)/c_per_mus;
        Y_param[i]           = y_loc*np.sin(az1[i])*np.cos(alt1[i])/c_per_mus#*np.cos(dec*np.pi/180)*np.sin(ha*np.pi/180)/c_per_mus;
        Z_param[i]           = z_loc*np.sin(alt1[i])/c_per_mus#*np.sin(dec*np.pi/180)/c_per_mus;
        W_geometric[i]       = (x_loc*np.sin(az)+y_loc*np.cos(az))*np.sin(za)/c_per_mus##(X_param[i]+Y_param[i]+Z_param[i])


    return W_geometric, alt1, az1, ha1
'''
cpdef delay_delay_rate_old(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, np.ndarray time_array, unsigned int T1, unsigned int T2, np.ndarray spec, CF, coarse_ran=20, smooth_ran=0.1):
'''
    #Module to calculate the delay delay_rate correlation function
    #INPUT:  np.ndarray spec, np.ndarray delay, coarse_ran=20, smooth_ran=0.1
'''
    cdef BW		 = 16.5
    cdef np.ndarray comp = spec.copy()
    cdef np.ndarray f	 = np.linspace(CF-BW/2, CF+BW/2, 256)
    delay   		 =    geometric_model_using_setdelay(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2) 
    comp_phase		 =    np.zeros((40, 60, len(spec[0])), dtype=float)
    ph_index		 =    np.zeros((40, 60), dtype=float)
    for i in range(-20, 20, 0.001):
    #Coarse delay calculation
        delay1	=    delay[0]+i	
        for j in range(0, 60, 0.5):
            delay2  	=  delay1*i
            ph      	=  genphase(spec, delay2, f) 
            comp        =  ph*spec
            comp_phase[i][j]  =  _header_gps_cy.obsdelay(comp)
            print('Phase std dev and mean for '+str(i)+' delay and '+str(j)+' rate is..'+str(std(comp_phase[0]))+' and '+str(mean(comp_phase[0])))
            ph_index[i][j]= np.mean(comp_phase[i][j]-np.mean(comp_phase[i][j])) 
    #cdef float low_index = np.argmin(ph_index)
    return comp_phase, ph_index
'''
cpdef delay_delay_rate(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2, spec, CF, file_name, file_name1, pol, sf_rate=0.000001, delay_rate=1.0, sf_points=04, dr_points=40):
    '''
            Module to calculate the delay delay_rate correlation function
            INPUT:  np.ndarray spec, np.ndarray delay, coarse_ran=20, smooth_ran=0.1
    '''
    try:
        fil     =       file_name.split('/')[-1]
        fil1    =       file_name1.split('/')[-1]

    except:
        fil     =       file_name
        fil1    =       file_name1



    cdef str series			 = fil.split('_')[-1]
    cdef float BW              		 = 16.5#MHz    
    cdef np.ndarray comp		 = spec.copy()
    cdef np.ndarray f    		 = np.linspace(CF-BW/2, CF+BW/2, 256)
    #delay                		 =    _header_geometric_cy.uvwsim_evaluate_baseline_uvw(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2)#geometric_model_using_setdelay(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2)
    
    #delay                               =    _header_geometric_cy.uvwsim_evaluate_baseline_uvw(file_name, file_name1, RA, dec, avg, time_array, T1, T2)
    #print(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2)
    delay				=	time_array
    delay_1				=     delay-min(delay)
    #cdef np.ndarray delay_1              =    (delay -round(min(delay), 8))#delay[1] - min(delay[1]) # us
    cdef float delay_fact		 =    np.nanmean(delay_1)
    #cdef np.ndarray comp_phase           =    np.zeros((sf_points, dr_points, len(spec[0])), dtype=float)
    cdef np.ndarray ph_index             =    np.zeros((sf_points, dr_points), dtype=float)
    delay_prac				 =    obsdelay(spec, CF)[0]

    #using line regress to fix the mismatch between the practical and
    #estimated delay..
    #if(series=='000.mbr'):
    line_x		=		np.linspace(0, len(delay_1), len(delay_1))
    line_est	=		linregress(line_x, delay_1*2*np.pi*CF)
    line_prac	=		linregress(line_x, delay_prac)
    print(line_est)
    print(line_prac)
    print(line_prac.intercept)
    slp_val		=		line_prac.slope/line_est.slope
    np.savetxt('SAMPLING_INFO/Delay_param_'+str(pol)+'_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.dparam', np.array([line_prac.intercept, slp_val, 0]))
    np.savetxt('SAMPLING_INFO/Delay_param_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.std_dpspec', ph_index)
    #else:
    #    slp_val		=		np.loadtxt('SAMPLING_INFO/Delay_param_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.dparam')[1]
    

    #delay1		=		delay_1*slp_val
    #ph          	=		genphase(spec, delay1, f)
    #comp		=  		ph*spec
    #ri,ci		=		0,0
    #for i in range(-1*sf_points/2, +1*sf_points/2):
    #Coarse delay calculation
    #    print('Now is sf_point no...'+str(i))
    #    delay1  =    delay_1+i*sf_rate
    #    for j in range(0, dr_points):
    #        delay2      =  delay1*j*delay_rate
    #        ph          =  genphase(spec, delay2, f)
    #        comp        =  ph*spec
    #        comp_phase  	=  obsdelay_delay_rate(comp)[0]#angle(mean(comp, axis =0))#obsdelay(comp)[0]
    #        ph_index[i][j]	=np.std(comp_phase)#np.mean(comp_phase[i][j]-np.mean(comp_phase[i][j]))
    #cdef float low_index = np.argmin(ph_index)
    #ri, ci = ph_index.argmin()//ph_index.shape[1], ph_index.argmin()%ph_index.shape[1] 
    #print('The index of ri, ci is..'+str(ri)+' , '+str(ci))
    #ri, ci = (ri-sf_points/2)*sf_rate, ci*delay_rate
    #np.savetxt('SAMPLING_INFO/Delay_param_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.dparam', np.array([0, slp_val, 0]))
    #np.savetxt('SAMPLING_INFO/Delay_param_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.std_dpspec', ph_index)
    
    #return (line_prac.intercept-line_est.intercept), slp_val, delay[1]#ri, ci#np.array([ri, ci])
    #return (-1*line_prac.intercept/line_prac.slope), slp_val, delay[1]#ri, ci#np.array([ri, ci])
    ret_prac		=	line_prac.slope*line_x + line_prac.intercept
    return (line_prac.intercept-line_est.intercept), slp_val, delay[1]
'''
cpdef delay_delay_rate_ason_JUN20_2021(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2, spec, CF, file_name, file_name1, sf_rate=0.000001, delay_rate=1.0, sf_points=4, dr_points=40):
'''
#            Module to calculate the delay delay_rate correlation function
#            INPUT:  np.ndarray spec, np.ndarray delay, coarse_ran=20, smooth_ran=0.1
'''
    try:
        fil     =       file_name.split('/')[-1]
        fil1    =       file_name1.split('/')[-1]

    except:
        fil     =       file_name
        fil1    =       file_name1



    cdef int i, j
    cdef float 	    BW                   = 16.5#MHz
    cdef np.ndarray comp                 = spec.copy()
    cdef np.ndarray f                    = np.linspace(CF-BW/2, CF+BW/2, 256)
    #Swapping frequency axis#
    #f					 = f[::-1]
    delay                                = _header_geometric_cy.uvwsim_evaluate_baseline_uvw(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2)#geometric_model_using_setdelay(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2)
    cdef np.ndarray delay_1              = delay[1] - min(delay[1]) # us
    cdef float delay_fact           	 = np.nanmean(delay_1)
    cdef np.ndarray ph_index             = np.zeros((sf_points, dr_points), dtype=float)
    for i in range(-1*sf_points/2, +1*sf_points/2):
    #Coarse delay calculation 
        print('Now is sf_point no...'+str(i))
        delay1  =    delay_1+i*sf_rate
        for j in range(0, dr_points):
            delay2      =  delay1*j*delay_rate
            ph          =  genphase(spec, delay2, f)
            comp        =  ph*spec
            comp_phase  =  obsdelay_delay_rate (comp, CF)[0]#
            ph_index[i][j]      =np.std(comp_phase)#
    ri, ci = ph_index.argmin()//ph_index.shape[1], ph_index.argmin()%ph_index.shape[1]
    ri, ci = (ri-sf_points/2)*sf_rate, ci*delay_rate
    print('The index of ri, ci is..'+str(ri)+' , '+str(ci))
    np.savetxt('SAMPLING_INFO/Delay_param_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.dparam', np.array([ri, ci, delay_fact]))
    np.savetxt('SAMPLING_INFO/Delay_param_'+str(fil[:-7])+'000.mbr_'+str(fil1[:-7])+'000.mbr.std_dpspec', ph_index)
    return ri, ci
'''

cpdef obsdelay_delay_rate(creal, CF):
    cdef np.ndarray phase1      =   np.zeros((len(creal[0])), dtype=float)
    cdef np.ndarray ampmax      =   np.zeros((len(creal[0])), dtype=float)
    fact        =   16/256.0
    fact        =   1/(fact*256.0)   #10**3 is the resampling factor..with 1000 fact come to be around 0.0625ns
    for i in range(len(creal[0])):
            z            =      np.hstack((creal[:,i], np.zeros((256*(10**3-1)))))
            corr         =      np.roll(np.fft.ifft(z)/np.sqrt(len(z)), len(z)/2)
            ampcorr      =      abs(corr)
            ampmax[i]    =      (np.argmax(ampcorr) - len(z)/2)
            phase1[i]    =      np.angle(complex(corr.real[np.argmax(ampcorr)],corr.imag[np.argmax(ampcorr)]))
    ampmax      =       ampmax - min(ampmax)
    phase1      =       (np.unwrap(phase1)-min(np.unwrap(phase1)))/(2*np.pi*CF)
    return phase1, ampmax


def obsdelay_ason_OCT23_2021(creal):
    cdef np.ndarray phase1      =   np.zeros((len(creal[0])), dtype=float)
    cdef np.ndarray ampmax      =   np.zeros((len(creal[0])), dtype=float)
    fact        =   16/256.0
    fact        =   1/(fact*256.0)   #10**3 is the resampling factor..with 1000 fact come to be around 0.0625ns
    for i in range(len(creal[0])):
            z            =      np.hstack((creal[:,i], np.zeros((256*(10**3-1)))))
            corr         =      np.roll(np.fft.ifft(z)/np.sqrt(len(z)), len(z)/2)
            corr1        =      corr.real
            corr2        =      corr.imag
            ampcorr      =      abs(corr)
            ampmax[i]    =      (np.argmax(abs(corr)) - len(z)/2)
            phase1[i]    =      np.angle(complex(corr1[np.argmax(ampcorr)],corr2[np.argmax(ampcorr)]))
    #ampmax      =       ampmax - min(ampmax)
    #Added unwrap, to avoid the problems with the practical curve fit routine..
    phase1	=	np.unwrap(phase1, discont=2.6)
    return phase1, ampmax



cdef tuple obsdelay(creal, CF):
    cdef np.ndarray phase1      =   np.zeros((len(creal[0])), dtype=float)
    cdef np.ndarray ampmax      =   np.zeros((len(creal[0])), dtype=float)
    fact        =   16/256.0
    fact        =   1/(fact*256.0)   #10**3 is the resampling factor..with 1000 fact come to be around 0.0625ns
    for i in range(len(creal[0])):
            z            =      np.hstack((creal[:,i], np.zeros((256*(10**3-1)))))
            corr         =      np.roll(np.fft.ifft(z)/np.sqrt(len(z)), len(z)/2)
            corr1        =      corr.real
            corr2        =      corr.imag
            ampcorr      =      abs(corr)
            ampmax[i]    =      (np.argmax(abs(corr)) - len(z)/2)
            phase1[i]    =      np.angle(complex(corr1[np.argmax(ampcorr)],corr2[np.argmax(ampcorr)]))
    #ampmax      =       ampmax - min(ampmax)
    #Added unwrap, to avoid the problems with the practical curve fit routine..
    phase1      =       np.unwrap(phase1, discont=2.6)
    #sending delay vaules in seconds 

    return phase1/(2*np.pi*CF), ampmax



cpdef gen_comp_spec(RA, dec, avg, T1, T2, spec, CF, file_name, file_name1, pol, file_time):
    #Get the comp spec out..
    #Generating delay
    #Checking if any previous delay_delay_rate files exsists
    try:
        fil     =       (file_name.split('/')[-1])[:-7]+'000.mbr'
        fil1    =       (file_name1.split('/')[-1])[:-7]+'000.mbr'

    except:
        fil     =       file_name[:-7]+'000.mbr'
        fil1    =       file_name1[:-7]+'000.mbr'
    
    cdef int series		= int(file_name[-7:-4])
    cdef float BW              	= 16.5    
    cdef np.ndarray comp	= spec.copy()
    cdef np.ndarray	f    	= np.linspace(CF-BW/2, CF+BW/2, 256)
    cdef unsigned int sec_junk, minu_junk, day, month, year
    

    #Extracting time..
    test			= _header_Fring_cy.extract_time(file_name) 
    
    day				= test[-3]
    month			= test[-2]
    year			= test[-1]
   
    file_time_1             = max(file_time[:,0])
    file_time_2             = max(file_time[:,1]) + max(file_time[:,0]) 
    
    if(int(test[3]) >= 12):
        print('The time is after 12 noon..')
        file_time_1             = max(file_time[:,0]) + 12*60*60
        file_time_2		= max(file_time[:,1]) + max(file_time[:,0]) + 12*60*60

    cdef long long int hour     = int(file_time_1/3600.0)
    cdef long long int minu     = int(file_time_1/60.0%60)
    cdef float         sec      = file_time_1%60.0

    print('The starting time for the delay calculation is..'+str(file_time_1))
    print('The end time of the delay calculation is...'+str(file_time_2))
    cdef np.ndarray time_array  = np.linspace(file_time_1/3600.0, file_time_2/3600.0, len(spec[0]))
    #Calculating Theoritical Delay..
    delay       		= _header_geometric_cy.uvwsim_evaluate_baseline_uvw(file_name, file_name1, RA, dec, avg, time_array, T1, T2)
    cdef np.ndarray delay_1     = delay# - min(delay)

    #delay_1			= delay_1*2*np.pi*CF 
    
    #Remember all the delay values are in radians.. 

    #Calculating Practical Delay..
    test1			= obsdelay(spec, CF)


    #Fitting a straight line to the practical and theoritical delay curve...
    lin_x                       = np.linspace(0, len(test1[0])-1, len(test1[0]))
    lin_prac			= linregress(lin_x, test1[0])
    #delay_prac                 = test1[0] - lin_prac.intercept
    #lin_prac                    = linregress(lin_x, delay_prac)
    lin_est                     = linregress(lin_x, delay_1)

    #Subtracting the intercept of the fitted straight line to avoid situation where there might be an momentary
    #zero dip, by the delay..
    

    #Comensating inital factor, after fitting a straight line to the practical delay..
    phase_shift			= lin_prac.intercept-lin_est.intercept 
    slp_val             	= lin_prac.slope/lin_est.slope
    
    print('lin_prac, lin_est')
    print(lin_prac, lin_est) 
    #np.save('delay_1', delay)
    #Checking if any previous delay_delay_rate files exsists
   
    #Working
    #delay_1			= delay_1*slp_val+phase_shift#%(2*np.pi)
    delay_1                     = delay_1+phase_shift#%(2*np.pi)
    #delay_1                     = delay_1+phase_shift%(2*np.pi)


    cdef np.ndarray delay_2			=	delay_1#*Delay_rate+Delay_shift#*Delay_rate+Delay_shift#(delay_1+Delay_shift)*Delay_rate
    np.save('delay_2', delay_2)
    cdef np.ndarray ph				=	genphase_phase(spec, delay_2, f)
    return spec, spec*ph, ph
'''
cpdef phase_comp_saved(sec, minu, hour, day, month, year, RA, dec, avg, time_array1, T1, T2, CF, pol):
   
'''
#        This module is used to compensate the fractional delay in the saved data set.
#        Input:
#            sec, minu, hour, day, month, year, RA, dec, avg, time_array1, T1, T2, spec, CF, file_name, file_name1, pol
#        The time_array here is a dummy variable, though the 'sec, minu, hour, day, month, year' should be accurate.
#        Output:
#             Compensated Spectrum, Phase Spectrum used to compensate.
''' 
    try:
        os.system('rm'+ ' Correlation_'+str(pol)+'_spec.npy')
    except:
        print('No older file found..')
        pass;

    files   	  =	glob.glob('Correlation_'+str(pol)+'_*', )
    print(files)
    cdef file_num =	len(files)	
    print('file_num..'+str(file_num))
    files.sort(key=os.path.getmtime)
    
    cdef np.ndarray spec#=   np.zeros((10)) # required for compilation!
    cdef np.ndarray sp  #=   'spec = np.hstack(('
    cdef list  sp1 	=   []
    cdef int   i   	=   0
    cdef float BW 	=   16.5
    for i in range(file_num):
        sp	=	np.load(files[i])
        if(i>0):
            spec	=	np.hstack((spec, sp))
        else:
            spec		=	sp.copy()
    print(spec)
    print(len(spec), len(spec[0]))
    cdef np.ndarray time_array  = 	np.linspace(hour+minu/60.0+sec/3600.0, hour+minu/60.0+sec/3600.0+(30*file_num)/3600.0, len(spec[0]))
    delay			=	_header_geometric_cy.uvwsim_evaluate_baseline_uvw(sec, minu, hour, day, month, year, RA, dec, avg, time_array, T1, T2)
    cdef np.ndarray     f       =	np.linspace(CF-BW/2, CF+BW/2, 256)	
    cdef np.ndarray delay1	=	delay[1] - min(delay[1])
    cdef np.ndarray ph		=	genphase(spec, delay1, f)
    cdef np.ndarray comp	=	ph*spec
    
    return comp, ph
''' 

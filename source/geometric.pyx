cimport numpy as np
import numpy as np
import math
from scipy import fftpack
import scipy

cdef double julian(float dd,float mm, float yr):
    """
    Function for Julian Day Calculation.

    Args:
        dd (float): day
        mm (float): month
        yr (float): year

    Returns:
        - **julday** (float) # Julian days
    """
    cdef double yr1    =   0
    cdef double jy     =   0
    cdef double jm, p, q, r, s 

    if (yr <  0 ):
        yr1 = yr + 1;
        print("Input Year wrong\n")
    if (mm > 2):
        jy = yr1;
        jm = mm+1
    else:
        jy = yr -1;
        jm = mm + 13
    p = jy * 365.25;
    q = jm * 30.6001;
    r = jy * 0.01;
    s = r * 0.25;
    julday = p + q + s + dd - r + 1720995 + 2;
    return(julday);
                                                                        
cpdef float selflst(float second, float minute, float hour, float day, float month, float year):
    """
    Function for Local Sidereal Time Calculation.

    Args:
        second (float): # Second 
        minute (float): # Minute 
        hour (float): # Hour 
        day (float): # Day 
        month (float): # Month 
        year (float): # Year

    Return:
        - **lst** (float) Local Sidereal Time (hour) 
    """

    cdef float tel_long_hr = 77.451944444/15
    cdef float dd = day
    cdef float mm = month
    cdef float yr = year
    cdef float ist          = hour + (minute / 60.0) + (second / 3600.0);
    #print('IST..'+str(ist))
    cdef float ut           = ist - 5.5;
    cdef float ist_ref_sec  = ist*3600.0;
    cdef float julday       = julian(dd,mm,yr);
    cdef float jd           = float(julday - 0.5);        #              /* JD at 0h UT */
    cdef float T            = (jd - 2451545.0) / 36525.0; # Time interval since 2000 Jan 1 12h UT*/
    gmst0        = 24110.54841 + 8640184.812866 * T + 0.093104 * T * T- 0.000006200 * T * T * T;
                                                #    /* Convert Gmst to Hours */
    gmst0 = gmst0 / (86400.0);     # /* gmst0 in days */
    gmst0 = np.modf(gmst0);#/* Get the fraction of a day */
    gmst1 = gmst0[1]
    gmst01 = gmst0[0]
    gmst0  = gmst01
    gmst0 = gmst0 * 24.0; #/* Convert into hours */
    if(gmst0 < -0.0000001):
        gmst0 = gmst0 + 24.0;
    cdef float mst = ut * 1.0027379094;
    if(mst-24.0 > 0.000001):
        mst = mst - 24.0;
#                  /* Greenwich mean sidereal time at required UT */
    cdef float gmst = gmst0 + mst;
    if(gmst-24.0 > 0.000001):
        gmst = gmst - 24.0;
    cdef float lst = gmst + tel_long_hr;
    if(lst-24.0 > 0.000001):
        lst = lst - 24.0;
    if(lst < -0.0000001):
        lst = lst + 24.0;
    cdef float lst_ref_sec = lst*3600.0;
    lsth = lst;
    lstm1 = (lst - lsth) * 60.0;
    lstm = lstm1;
    lsts = (lstm1 - lstm) * 60.0;
    return lst





cpdef np.ndarray Cal_time(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2):
    
    '''
    float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2
    '''
    cdef double c_per_mus   = 299.792458
    cdef int    i       =       0
    print (year, month, day, hour, minu, sec)
    cdef np.ndarray X_param     =      np.zeros(int(time/del_t), dtype = float)
    cdef np.ndarray Y_param     =      np.zeros(int(time/del_t), dtype = float)
    cdef np.ndarray Z_param     =      np.zeros(int(time/del_t), dtype = float)
    cdef np.ndarray W_geometric =      np.zeros(int(time/del_t), dtype = float)
    #Gettine ECEF Coordinates of Tiles#

    coslat               = np.cos(13.6111*np.pi/180)
    sinlat               = np.sin(13.6111*np.pi/180)#77.451944444*np.pi/180)

    ecef    =   np.loadtxt('ENU_v2.txt')
    x_loc   =   (ecef[:,0][T1] - ecef[:,0][T2])*sinlat+(ecef[:,2][T1] - ecef[:,2][T2])*coslat#np.loadtxt('ECEF_x')
    y_loc   =   (ecef[:,1][T1] - ecef[:,1][T2])#np.loadtxt('ECEF_y')
    z_loc   =   (ecef[:,2][T1] - ecef[:,2][T2])*coslat+(ecef[:,2][T1] - ecef[:,2][T2])*sinlat#np.loadtxt('ECEF_z')
    print(x_loc, y_loc, z_loc)
    cdef double secu        =   sec
    cdef double mint        =   minu
    cdef double hourt       =   hour
    cdef double dayt        =   day
    cdef double sindec      =   np.sin(np.pi/180*dec)
    cdef double cosdec      =   np.cos(np.pi/180*dec)

    for i in range(int(time/del_t)):
        #print(selflst(sec+i*del_t, minu, hour, day, month, year))
        secu    =   secu+1
        #secu, mint, hourt, dayt = fix_time(secu, mint, hourt, dayt)
        if( secu > 59):
            mint    =   mint+1
            secu    =   0
        if(mint>59):
            mint    =   0
            hourt   =   hourt+1
        if(hourt>23):
            hourt   =   0
            dayt    =   dayt+1
        #print(str(secu)+':'+str(mint)+':'+str(hourt))
        tim     =   secu/3600.0+mint/60.0+hourt
        #print(tim, (selflst(secu, mint, hourt, day, month, year)))
        ha                      = (selflst(secu, mint, hourt, day, month, year)-RA)*15.0*np.pi/180.0
        print(i, RA, tim, selflst(secu, mint, hourt, day, month, year), (selflst(secu, mint, hourt, day, month, year)-RA)*15.0*np.pi/180, tim)
        X_param[i]           = x_loc*cosdec*np.cos(ha)/c_per_mus;
        Y_param[i]           = -1.0*y_loc*cosdec*np.sin(ha)/c_per_mus;
        Z_param[i]           = z_loc*sindec/c_per_mus;  
        print(X_param[i], Y_param[i], Z_param[i])
        W_geometric[i]       = (X_param[i]+Y_param[i]+Z_param[i])
    return W_geometric


cdef fix_time(sec, mi, hr, day):
    cdef int trk =   0
    cdef int fact=   1

    hour    =   int(sec/3600.0)+hr
    minu    =   (sec/3600.0%1)*60+mi
    sece    =   (minu%1)*60

    if(sece > 59):
        while(fact):
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


cpdef np.ndarray Cal_time_onhold_v1(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2):

    '''
    float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2
    '''
    cdef double c_per_mus   = 299.792458
    cdef int    i       =       0
    print (year, month, day, hour, minu, sec)
    cdef np.ndarray X_param     =      np.zeros(int(time/del_t), dtype = float)
    cdef np.ndarray Y_param     =      np.zeros(int(time/del_t), dtype = float)
    cdef np.ndarray Z_param     =      np.zeros(int(time/del_t), dtype = float)
    cdef np.ndarray W_geometric =      np.zeros(int(time/del_t), dtype = float)
    cdef float alt              =       0
    cdef float az               =       0

    #Gettine ECEF Coordinates of Tiles#
    cdef double coslat               = np.cos(np.pi/180*13.6033333)#np.cos(latlong[T1][0]*np.pi/180)
    cdef double sinlat               = np.sin(np.pi/180*13.6033333)#np.sin(latlong[T1][0]*np.pi/180)
    cdef double cosdec               = np.cos(np.pi/180*dec)
    cdef double sindec               = np.sin(np.pi/180*dec)
    cdef double cosha                = 0
    cdef double sinha                = 0
   

    #Commented
    enu    =   np.loadtxt('ENU_v2.txt')#np.loadtxt('SWAN_ENU.txt')
    #ecef    =   np.loadtxt('ecef.txt')
    ref_east_dist   =   enu[T1][0]
    ref_north_dist  =   enu[T1][1]
    ref_up_dist     =   enu[T1][2]
    east_dist       =   enu[T2][0]
    north_dist      =   enu[T2][1]
    up_dist         =   enu[T2][2]
    ##########     EEEEEEEEEEEE
    #       ##     EEE
    #       ##     EEE
    #       ##     EEE
    #########      EEEEEEEEEEEE
    #       ##     EEE
    #        ##    EEE
    #         ##   EEE
    #          ##  EEEEEEEEEEEE

    #x_loc   =   enu[:,0][T1] - enu[:,0][T2]#np.loadtxt('ECEF_x')
    #y_loc   =   enu[:,1][T1] - enu[:,1][T2]#np.loadtxt('ECEF_y')
    #z_loc   =   enu[:,2][T1] - enu[:,2][T2]#np.loadtxt('ECEF_z')


    x_loc   =    -1*(north_dist-ref_north_dist)*sinlat + (up_dist-ref_up_dist)*coslat;#np.loadtxt('ECEF_x')
    y_loc   =   east_dist - ref_east_dist;#np.loadtxt('ECEF_y')
    z_loc   =   (north_dist-ref_north_dist)*coslat + (up_dist-ref_up_dist)*sinlat#np.loadtxt('ECEF_z')
    print(x_loc, y_loc, z_loc)
    cdef double secu        =   sec
    cdef double mint        =   minu
    cdef double hourt       =   hour
    cdef double dayt        =   day
    for i in range(int(time/del_t)):
        #print(selflst(sec+i*del_t, minu, hour, day, month, year))
        secu    =   secu+1
        if( secu > 59):
            mint    =   mint+1
            secu    =   0
        if(mint>59):
            mint    =   0
            hourt   =   hourt+1
        tim     =   i*del_t/3600.0+minu/60.0+hour
        #ha                   = (selflst(sec+i*del_t, minu, hour, day, month, year)-RA)*15.0
        ha                   = (selflst(secu, mint, hourt, dayt, month, year)-RA)*15.0
        coslat               = np.cos(np.pi/180*13.6033333)#np.cos(latlong[T1][0]*np.pi/180)
        sinlat               = np.sin(np.pi/180*13.6033333)#np.sin(latlong[T1][0]*np.pi/180)
        
        cosha                = np.cos(ha*np.pi/180)
        sinha                = np.sin(ha*np.pi/180)
        alt                  = np.arcsin(np.sin(dec*np.pi/180)*sinlat+np.cos(dec*np.pi/180)*coslat*np.cos(ha*np.pi/180))
        az                   = np.arcsin(-1*np.sin(ha*np.pi/180)*np.cos(dec*np.pi/180)/np.cos(alt))
        alt1                 =  np.pi/2-alt
        
        X_param[i]           = (x_loc*cosdec)/c_per_mus;#x_loc*np.sin(az)*np.cos(alt1)#*np.cos(dec*np.pi/180)*np.cos(ha*np.pi/180)/c_per_mus;
        Y_param[i]           = (-1.0*y_loc*cosdec)/c_per_mus#*np.cos(dec*np.pi/180)*np.sin(ha*np.pi/180)/c_per_mus;
        Z_param[i]           = (z_loc*sindec)/c_per_mus;#z_loc*np.sin(alt)#*np.sin(dec*np.pi/180)/c_per_mus;
        W_geometric[i]       = X_param[i]*cosha + Y_param[i]*sinha + Z_param[i];
    return W_geometric

cpdef phase_compensation(spect, delay, freq):
    '''
    Compensation for intra sample delay.
    
    Args:
        spect: spectrum
        delay: delay
        freq: frequency (MHz)
    
    Returns:
        - **comp** (array) - compensated spectrum
        - **pha** (array) - phase compensation spectrum
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

cpdef phase_compensation_toggle(sec1, min1, hour1, day1, month, year, RA,  Dec, avg, del_t, time, T1, T2):
    """
    ..warning::
        desc required!!
    
    Args:
        sec1
        min1
        hour1
        day1
        month
        year
        RA
        Dec
        avg
        del_t
        time
        T1
        T2

    Returns:
        - **delay** 
        - **fact** 
        - **delay1**
    """ 

    delay   =   Cal_time_onhold_v1(float(sec1), float(min1), float(hour1), float(day1), float(month), float(year), float(RA),  float(Dec), float(avg), del_t, time, T1, T2)
    cdef double fact            =   sum(delay)*10**-6*33000000
    cdef np.ndarray delay1      =   delay-min(delay)
    if(fact<1):
        fact    =   0
        print('Intra sample delay correction required..')
        return delay, fact, delay1
    else:
        fact    =   1
        print('Inter sample delay correction required..')
        return delay, fact, delay1

cpdef obsdelay(creal9):
    cdef np.ndarray phase1          =   np.zeros((len(creal9[0])))
    cdef np.ndarray phase_total     =   np.zeros((len(creal9), len(creal9[0])))
    cdef int i      =   0
    for i in range(len(creal9[0])):
        zemp             = np.append(creal9[:,i], np.zeros((256*9)))
        icorr            = scipy.fftpack.ifft(zemp)/np.sqrt(len(zemp))
        amp              = abs(icorr)
        amp1             = np.roll(amp, len(amp)/2)
        phase_total[:,i] = amp1[1152:1408]
        phase1[i]        = np.angle(complex(icorr[np.argmax(amp)].real,icorr[np.argmax(amp)].imag))*180/np.pi

    return phase1, phase_total
"""
cpdef phasecomp(creal9, delay, cenfreq, fftlen):
    '''
    Used when the delay is less than one sample, i.e intra sample level compensation..

    Args:
        creal9: creal9 (Rename This!!) 
        delay: Delay
        cenfreq: Central Frequency
        fftlen: FFT Length
    '''
    #From shift theorem, x(t-t0)->X(F)*exp(2pifoto)
    cdef np.ndarray ex         =   np.zeros((len(creal9)), dtype=np.float64)
    cdef np.ndarray creal9_1   =   creal9.copy()
    cdef np.ndarray creal9_2   =   creal9.copy()
    cdef np.ndarray freq       =   np.linspace(cenfreq-16.5/2, cenfreq+16.5/2, fftlen, dtype=np.float64)
    cdef fact                  =   np.zeros((len(creal9)), dtype=np.complex128)   
    cdef int i              =   0
    cdef int j              =   0

    for i in range(len(creal9[0])):
            fact    =   np.cos(2*np.pi*freq*delay[i])-1j*np.sin(2*np.pi*freq*delay[i])
            creal9_1[:,i]   =   creal9[:,i]*fact
            creal9_2[:,i]   =   fact
    return creal9_1, creal9_2
"""


"""
def ECEFto_ENU_file(file_path):
    ecef    =   np.loadtxt(file_path)
    enu     =   []
    for i in range(len(ecef)):
        print(i)
        enu.append(EcefToEnu(ecef[i][0], ecef[i][1], ecef[i][2]))
    return enu
"""
def EcefToEnu(x, y, z):
    """
    Convert ECEF to ENU Coordinate System.

    Args:
        x (float): ECEF Coordinate X.
        y (float): ECEF Coordinate Y.
        z (float): ECEF Coordinate Z.


    :rtype: xEast (float), yNorth (float), zUpfloat (float)

    Returns:
        ENU Coordinates (East, North, Up)

    """

    #// Convert to radians in notation consistent with the paper:

    lat0    =   13.6033333
    lon0    =   77.451944444 
    h0      =   483
    a       = 6378137.0;      #   // WGS-84 Earth semimajor axis (m)
    b       = 6356752.314245; #    // Derived Earth semiminor axis (m)
    f       = (a - b) / a;    #       // Ellipsoid Flatness
    f_inv   = 1.0 / f;    #   // Inverse flattening

    a_sq    = a * a;
    b_sq    = b * b;
    e_sq    = f * (2 - f); 
    lambd   = lat0*np.pi/180

    phi     = lon0*np.pi/180
    s       = np.sin(lambd);
    N       = a/np.sqrt(1-e_sq*s*s);

    sin_lambda = np.sin(lambd);
    cos_lambda = np.cos(lambd);
    cos_phi = np.cos(phi);
    sin_phi = np.sin(phi);

    x0 = (h0 + N) * cos_lambda * cos_phi;
    y0 = (h0 + N) * cos_lambda * sin_phi;
    z0 = (h0 + (1 - e_sq) * N) * sin_lambda;

    #xd, yd, zd;
    xd = x - x0;
    yd = y - y0;
    zd = z - z0;

    #// This is the matrix multiplication
    xEast = -sin_phi * xd + cos_phi * yd;
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd;
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd;
    
    return xEast, yNorth, zUp


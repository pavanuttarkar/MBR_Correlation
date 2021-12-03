cimport numpy as np
import numpy as np
import math
from scipy import fftpack
import scipy
import sys

#Including gloabal constants, GBD Lat and Long
import LatLong


import _header_Fring_cy


cpdef double julian (float dd,int mm,int yr):
    """
    Function to  compute julian date(JD) ## Check here MJD or JD

    Parameters
    ----------
        dd: `float`
            Days
        mm: `float`
            Month
        yr: `float`
            Year

    Returns
    -------
        julday: `float`
            Julian day
    """
    cdef double jy,jm,r,s;
    cdef double p,q,ja,julday;
    cdef double A, B, C,E, F
    if (yr <  0 ):
        yr = yr + 1;
        print("Input Year wrong\n")
    if (mm > 2):
        jy = float(yr);
        jm = float(mm)#+1;
    else:
        jy = yr -1;
        jm = mm + 12;

    A = int(jy/100.0)
    B = int(A/4.0)
    C = 2.0-A+B 
    E = int(365.25*(jy+4716))
    F = int(30.6001*(jm+1.0))
    julday= C+dd+E+F-1524.5
    return julday; 


cpdef float  selflst(float second, float minute, float hour, float day, float month, float year):
    """
    Function to compute local sidereal time.

    Parameters
    ----------
        second: `float`
        minute: `float`
        hour: `float`
        day: `float`
        month: `float`
        year: `float`

    Returns
    -------
        lst: `float`
            Local sidereal time
    """
    #cdef float hour     =    int(tim)
    #cdef float minute   =    (tim-int(tim))*60
    #cdef float second   =    math.modf(minute)[0]*3600.0
    cdef float tel_long_hr = 77.451944444/15
    cdef int dd = int(day)
    cdef int mm = int(month)
    cdef int yr = int(year)
    cdef float ist          = hour + (minute / 60.0) + (second / 3600.0);
    #print('IST..'+str(ist))
    cdef float ut           = ist - 5.5;
    cdef float ist_ref_sec  = ist*3600.0;
    cdef double julday       = julian(dd,mm,yr);
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


cpdef double uvwsim_datetime_to_mjd(int year, int month, int day, int hour, int minute, double seconds):
    cdef double day_fraction;
    cdef int a, y, m, jdn;

    #/* Compute Julian Day Number (Note: all integer division). */
    a = (14 - month) / 12;
    y = year + 4800 - a;
    m = month + 12 * a - 3;
    jdn = day + (153 * m + 2) / 5 + (365 * y) + (y / 4) - (y / 100)+ (y / 400) - 32045;

    #/* Compute day fraction. */
    day_fraction = (hour + minute/60. + seconds/3600.) / 24.;

    #/* Compute day fraction. */
    day_fraction -= 0.5;
    return (jdn - 2400000.5) + day_fraction;


#LST equation from uvw_sim
cpdef double uvwsim_convert_mjd_to_gmst(double mjd):

    cdef double d, T, gmst;

    #/* Days from J2000 */
    d = mjd - 51544.5;

    #/* Centuries from J2000 */
    T = d / 36525.0;

    #/* GMST at this time, in radians */
    gmst = 24110.54841 + 8640184.812866*T + 0.093104*T*T - 6.2e-6*T*T*T;
    gmst *= (15.0*np.pi)/(3600.*180.);
    gmst += math.fmod(mjd, 1.0) * 2.0 * np.pi;

    #/* Range check (0, 2pi) */
    T = math.fmod(gmst, 2*np.pi);
    if(T>0.0):
        return T
    else:
        return  T + 2.0 * np.pi


cpdef double uvwsim_equation_of_equinoxes_fast(double mjd):

    cdef double d, omega, L, delta_psi, epsilon, eqeq;

    #/* Days from J2000.0. */
    d = mjd - 51544.5;

    #/* Longitude of ascending node of the Moon. */
    omega = (125.04 - 0.052954 * d) * (np.pi/180.0);

    #/* Mean Longitude of the Sun. */
    L = (280.47 + 0.98565 * d) * (np.pi/180.0);

    #/* eqeq = delta_psi * cos(epsilon). */
    delta_psi = -0.000319 * np.sin(omega) - 0.000024 * np.sin(2.0 * L);
    epsilon = (23.4393 - 0.0000004 * d) * (np.pi/180.0);

    #/* Return equation of equinoxes, in radians. */
    eqeq = delta_psi * np.cos(epsilon) * ((15.0*np.pi)/180.0);
    return eqeq;

cpdef float  gmst(float second, float minute, float hour, float day, float month, float year):
    """
    Function to compute gmst

    Parameters
    ----------
        second: `float`
            Seconds
        minute: `float`
            Minutes
        hour: `float`
            Hours
        day: `float`
            Days
        month: `float`
            Month
        year: `float`
            Year

    Returns
    -------
        gmst: `float`
    """
    #cdef float hour     =    int(tim)
    #cdef float minute   =    (tim-int(tim))*60
    #cdef float second   =    math.modf(minute)[0]*3600.0
    cdef float tel_long_hr = 77.451944444/15
    cdef int   dd = int(day)
    cdef int   mm =int(month)
    cdef int   yr =int(year)
    cdef float ist          = hour + (minute / 60.0) + (second / 3600.0);
    #print('IST..'+str(ist))
    cdef float ut           = ist# - 5.5;
    cdef float ist_ref_sec  = ist*3600.0;
    cdef float julday       = julian(dd,mm,yr);
    cdef float jd           = float(julday - 0.5);        #              /* JD at 0h UT */
    cdef float T            = (jd - 2451545.0) / 36525.0; # Time interval since 2000 Jan 1 12hr UT*/
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
    return gmst


cpdef double uvwsim_convert_mjd_to_gast_fast(double mjd):
    cdef double gmst = uvwsim_convert_mjd_to_gmst(mjd);
    cdef double gast = gmst + uvwsim_equation_of_equinoxes_fast(mjd);
    return gast;


cpdef uvwsim_evaluate_baseline_uvw(str file_name, str file_name1, float RA, float dec, float avg, np.ndarray time_array, unsigned int T1, unsigned int T2):
    '''
    ECEF ENU Transformation..
    E  -  Y
    N  -  Z
    U  -  X 
    '''
    cdef double c_per_mus   = 299.792458
    cdef int    i       =       0
    coslat               = np.cos(13.6111*np.pi/180)
    sinlat               = np.sin(13.6111*np.pi/180)

    cdef np.ndarray ecef    =       np.loadtxt('ECEF_from_header_geometric.txt')#np.loadtxt('ECEF_from_header_geometric_Tile1_Ref.txt')#np.loadtxt('ECEF_v1.txt')
    cdef double x_ecef      =        ecef[:,0][T1]-ecef[:,0][T2]#ecef[:,1][T1]-ecef[:,1][T2]
    cdef double y_ecef      =        ecef[:,1][T1]-ecef[:,1][T2]#ecef[:,2][T1]-ecef[:,2][T2]
    cdef double z_ecef      =        ecef[:,2][T1]-ecef[:,2][T2]#ecef[:,0][T1]-ecef[:,0][T2]

    #cdef double x_ecef   =   -1*(ecef[:,0][T1] - ecef[:,0][T2])*sinlat+(ecef[:,2][T1] - ecef[:,2][T2])*coslat#*sinlat+(ecef[:,2][T1] - ecef[:,2][T2])*coslat#np.loadtxt('ECEF_x')
    #cdef double y_ecef   =   (ecef[:,1][T1] - ecef[:,1][T2])#np.loadtxt('ECEF_y')
    #cdef double z_ecef   =   (ecef[:,0][T1] - ecef[:,0][T2])*coslat+(ecef[:,2][T1] - ecef[:,2][T2])*sinlat#*coslat+(ecef[:,2][T1] - ecef[:,2][T2])*sinlat#np.loadtxt('ECEF_z')
    
    cdef unsigned int day, month, year
    series_junk, sec_junk, minu_junk, hour_junk, day, month, year       =       _header_Fring_cy.extract_time(file_name)

    cdef double gast  
    cdef double ha0         =       0 
    cdef double d2r         =       np.pi/180.0
    cdef np.ndarray W_geometric =      np.zeros(len(time_array), dtype = float)
    cdef list ha1         =       []
    #/* Convert to station uvw */
    cdef double sinha0  
    cdef double cosha0  
    cdef double sindec0 = np.sin(dec*d2r);
    cdef double cosdec0 = np.cos(dec*d2r);

    for i in range(len(time_array)):

        hourt_temp   =   int(time_array[i])
        mint_temp    =   (time_array[i]%1)*60
        secu_temp    =   (mint_temp%1)*60

        hourt       =   int(hourt_temp)
        mint        =   int(mint_temp)
        secu        =   secu_temp


        #print(secu, mint, hourt, day, month, year)        
        #time_mjd    =   uvwsim_datetime_to_mjd(int(year), int(month), int(day), int(hourt), int(mint),   secu)
        #gast        =   uvwsim_convert_mjd_to_gast_fast(time_mjd);

        #time_mjd    =   julian(day, int(month), int(year))#uvwsim_datetime_to_mjd(int(year), int(month), int(day), int(hourt), int(mint),   secu)
        #gast        =   gmst(secu, mint, hourt, day, month, year)#uvwsim_convert_mjd_to_gast_fast(time_mjd);


        lst		=	selflst(secu, mint, hourt, day, month, year)
        ha0		=	lst*15.0*d2r- RA*15.0*d2r
        
        #ha0  = gast*15.0*d2r - RA*15.0*d2r;
        #ha0  = gast - RA*15.0*d2r 
        #if(ha0<0.0):
        #    ha0 =   ha0+2*np.pi
        #if(ha0>np.pi):
        #    ha0 =   ha0-2*np.pi
        #ha0     =   ha0*15.0*np.pi/180.0
        #print(time_array[i], time_mjd, gast, ha0, i)
        #ha1.append(ha0)
        sinha0  = np.sin(ha0);
        cosha0  = np.cos(ha0);
        t = x_ecef * cosha0 - y_ecef*sinha0;
        W_geometric[i] = t * cosdec0 + z_ecef * sindec0;
    return W_geometric/c_per_mus


cpdef Cal_time(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, np.ndarray time_array, unsigned int T1, unsigned int T2):

    '''
    float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, np.array time_array, unsigned int T1, unsigned int T2
    
    .. note:: 
        Hour Angle convention
        -ve while rising from east
        +ve while setting in the west

    Hour angle convention and integrety check -- done on MAY06 2021

    '''

    cdef double c_per_mus   = 299.792458                                       
    cdef int    i       =       0                                              
    print (year, month, day, hour, minu, sec)                                                                      
    cdef np.ndarray X_param     =      np.zeros(len(time_array), dtype = float)                                    
    cdef np.ndarray Y_param     =      np.zeros(len(time_array), dtype = float)                                    
    cdef np.ndarray Z_param     =      np.zeros(len(time_array), dtype = float)                                    
    cdef np.ndarray W_geometric =      np.zeros(len(time_array), dtype = float)                                    
    #Gettine ECEF Coordinates of Tiles#                                                                                                                             
                                                                                                                                                                    
    coslat               = np.cos(13.6111*np.pi/180)                                                                                                                
    sinlat               = np.sin(13.6111*np.pi/180)                                                                                                                
                                                                                                                                                                    
    ecef    =   np.loadtxt('ENU_v2.txt')                                                                                               
    x_loc   =   -1*(ecef[:,0][T1] - ecef[:,0][T2])*sinlat+(ecef[:,2][T1] - ecef[:,2][T2])*coslat#*sinlat+(ecef[:,2][T1] - ecef[:,2][T2])*coslat#np.loadtxt('ECEF_x')
    y_loc   =   (ecef[:,1][T1] - ecef[:,1][T2])#np.loadtxt('ECEF_y')                                                                                             
    z_loc   =   (ecef[:,0][T1] - ecef[:,0][T2])*coslat+(ecef[:,2][T1] - ecef[:,2][T2])*sinlat#*coslat+(ecef[:,2][T1] - ecef[:,2][T2])*sinlat#np.loadtxt('ECEF_z')
                                   
    print(x_loc, y_loc, z_loc)      
    cdef double secu        =   sec 
    cdef double mint        =   minu
    cdef double hourt       =   hour                 
    cdef double dayt        =   day                  
    cdef double sindec      =   np.sin(np.pi/180*dec)
    cdef double cosdec      =   np.cos(np.pi/180*dec)
    ha1           =   []                             
    for i in range(len(time_array)):                 
        hourt_temp   =   int(time_array[i])          
        mint_temp    =   (time_array[i]%1)*60
        secu_temp    =   (mint_temp%1)*60    
                                             
        hourt       =   int(hourt_temp)      
        mint        =   int(mint_temp)   
        secu        =   secu_temp                                                       
                                                                                        
        ha1.append(selflst(secu, mint, hourt, day, month, year)-RA)
        if(ha1[i]<0):                                                                   
               ha1[i] = ha1[i]+24.0                                                     
               print('Negative ha..'+str(ha1[i]))                                       
        if(ha1[i]>12):                           
                ha1[i] = ha1[i]-24.0                   
        #ha1.append(ha)                                       
        ha                      = ha1[i]*15*np.pi/180.0                                                                              
        #print(i, selflst(secu, mint, hourt, day, month, year), RA, ha, str(ha*180/np.pi))                   
        X_param[i]           = x_loc*cosdec*np.cos(ha)/c_per_mus;                                                                    
        Y_param[i]           = -1.0*y_loc*cosdec*np.sin(ha)/c_per_mus;                                                               
        Z_param[i]           = z_loc*sindec/c_per_mus;                                                                               
        W_geometric[i]       = (X_param[i]+Y_param[i]+Z_param[i])                                                                    
    return W_geometric


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
    cdef double montht      =   month
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
        print(secu, mint, hourt, dayt, montht, year, az, alt, i)
        alt1.append(alt)
        az1.append(az)
        ha1.append(ha)
        za                   = np.pi-alt
        X_param[i]           = x_loc*np.cos(az1[i])*np.cos(alt1[i])/c_per_mus#*np.cos(dec*np.pi/180)*np.cos(ha*np.pi/180)/c_per_mus;
        Y_param[i]           = y_loc*np.sin(az1[i])*np.cos(alt1[i])/c_per_mus#*np.cos(dec*np.pi/180)*np.sin(ha*np.pi/180)/c_per_mus;
        Z_param[i]           = z_loc*np.sin(alt1[i])/c_per_mus#*np.sin(dec*np.pi/180)/c_per_mus;
        W_geometric[i]       = (x_loc*np.sin(az)+y_loc*np.cos(az))*np.sin(za)/c_per_mus##(X_param[i]+Y_param[i]+Z_param[i])
    return W_geometric, alt1, az1, ha1


cpdef  Cal_time_onhold(float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, np.ndarray time_array, unsigned int T1, unsigned int T2):
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

        #if(hourt_temp > 24.0):
        #    print('Increasing day..')
        #    dayt    =   day+1
        #    hourt   =   0
        #    minu    =   0
        #    secu    =   0
        #    dayflag =   1

        az, alt, ha = Equ2local(RA, dec, LatLong.Lat_local, LatLong.Long_local, secu, mint, hourt, day, month, year)
        print(secu, mint, hourt, dayt, montht, year, az, alt, i)
        alt1.append(alt)
        az1.append(az)
        ha1.append(ha)

        X_param[i]           = x_loc*np.cos(az1[i])*np.cos(alt1[i])/c_per_mus#*np.cos(dec*np.pi/180)*np.cos(ha*np.pi/180)/c_per_mus;
        Y_param[i]           = y_loc*np.sin(az1[i])*np.cos(alt1[i])/c_per_mus#*np.cos(dec*np.pi/180)*np.sin(ha*np.pi/180)/c_per_mus;
        Z_param[i]           = z_loc*np.sin(alt1[i])/c_per_mus#*np.sin(dec*np.pi/180)/c_per_mus;
        W_geometric[i]       = (X_param[i]+Y_param[i]+Z_param[i])
    return W_geometric, alt1, az1, ha1


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
    enu             =   np.loadtxt('ECEF_from_header_geometric.txt')#enu    =   np.loadtxt('ENU_v2.txt')#np.loadtxt('SWAN_ENU.txt')
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

    x_loc   =   enu[:,0][T1] - enu[:,0][T2]#np.loadtxt('ECEF_x')
    y_loc   =   enu[:,1][T1] - enu[:,1][T2]#np.loadtxt('ECEF_y')
    z_loc   =   enu[:,2][T1] - enu[:,2][T2]#np.loadtxt('ECEF_z')


    #x_loc   =    -1*(north_dist-ref_north_dist)*sinlat + (up_dist-ref_up_dist)*coslat;#np.loadtxt('ECEF_x')
    #y_loc   =   east_dist - ref_east_dist;#np.loadtxt('ECEF_y')
    #z_loc   =   (north_dist-ref_north_dist)*coslat + (up_dist-ref_up_dist)*sinlat#np.loadtxt('ECEF_z')
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


cpdef phase_compensation_dev(spect, delay, freq):
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
            pha[i][j] = (np.exp(complex(0, -2*np.pi*i*delay[j])))

    comp    =   np.flipud(pha)*spect

    return comp, np.flipud(pha)

cpdef phase_compensation_toggle(sec1, min1, hour1, day1, month, year, RA,  Dec, avg, del_t, time, T1, T2):
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
    cdef np.ndarray phase_total     =   np.zeros((len(creal9)*100, len(creal9[0])))
    cdef int i      =   0
    for i in range(len(creal9[0])):
        zemp             = np.append(creal9[:,i], np.zeros((256*99)))
        icorr            = scipy.fftpack.ifft(zemp)/np.sqrt(len(zemp))
        amp              = abs(icorr)
        amp1             = np.roll(amp, len(amp)/2)
        phase_total[:,i] = amp1#[1152:1408]
        phase1[i]        = np.angle(complex(icorr[np.argmax(amp)].real,icorr[np.argmax(amp)].imag))*180/np.pi
        
    #for i in range(len(phase1)-1):
    #    if ((phase1[i+1] - phase1[i]) > 180):
    #        phase1[i+1] = phase1[i+1] - 360
    #    elif ((phase1[i+1] - phase1[i]) < -180):
    #        phase1[i+1] = phase1[i+1] + 360
    #    else:
    #        i = i + 1
    #phase1 = phase1 - min(phase1)
    return phase1, phase_total, icorr


cpdef phasecomp(creal9, delay, cenfreq, fftlen):
    '''
        creal9, delay, cenfreq, fftlen
        Used when the delay is less than one sample, i.e intra sample level compensation..
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
            #fact    =   complex(np.cos(2*np.pi*freq[i]*delay[i]), np.sin(2*np.pi*freq[i]*delay[i]))
            #for j in range(len(creal9[0])):
            fact    =   np.cos(2*np.pi*freq*delay[i])-1j*np.sin(2*np.pi*freq*delay[i])
            creal9_1[:,i]   =   creal9[:,i]*fact
            creal9_2[:,i]   =   fact
    return creal9_1, creal9_2


def ECEFto_ENU_file(file_path):
    ecef    =   np.loadtxt(file_path)
    enu     =   []
    for i in range(len(ecef)):
        print(i)
        enu.append(EcefToEnu(ecef[i][0], ecef[i][1], ecef[i][2]))
    return enu


def EcefToEnu(x, y, z):#, lat0, lon0, h0):
    """
    Function to transform ECEF coordinates to ENU coordinates.
    Returns East, North and Up coordinates.

    Parameters
    ----------
        x: `float`
        y: `float`
        z: `float`
    
    Returns
    -------
        East: `float`
        North: `float`
        Up: `float`
    """
    #// Convert to radians in notation consistent with the paper:
    #Tile 1 as the reference
    lat0    =   13.60246667#13.6033333
    lon0    =   77.42778333#77.451944444
    h0      =   691#483
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

    #xd, yd, zd;# Current reference Tile 1.. 
    xd = x - 1349783.1183164874 #1347182.8954636795#1349783.1183164874#x0;
    yd = y - 6052369.1881498210 #6052720.426106435 #6052369.188149821#y0;
    zd = z - 1490432.8161248143 #1491361.0540407044#1490432.8161248143#x0#z0;

    #// This is the matrix multiplication
    xEast = -sin_phi * xd + cos_phi * yd;
    yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd;
    zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd;
    return xEast, yNorth, zUp


cpdef   WGS842ENU(double lat, double lon, double alt, double lat0, double lon0, double alt0):
    '''
    Function to cover WGS84 to ENU, with lat0, lon0, alt0 as the centre.
    all inputs should be in degrees.
    
    WGS842ENU(double lat, double lon, double alt)
    ENU (X, Y, Z) of lat, lon, alt, lat0, lon0, alt0  inputs.
    
    .. warning::
        Need to check geodetic and geocentric latitude calculation..not sure which latitude was
        measured for the estimation of tile positions..needs to be clarified. Now assuming what
        we have is geodetic latitude.
    
    Parameters
    ----------
        lat: `double`
            Latitute
        lon: `double`
            Longitude
        alt: `double`
            Altitude
        lat0: `double`
            Latitute of the center
        lon0: `double`
            Longitude of the center
        alt0: `double`
            Altitude of the center

    Returns
    -------
        E: `double`
            East
        N: `double`
            North
        U: `double` 
            Up
    '''
    cdef double E=0
    cdef double N=0
    cdef double U=0
    cdef double d2r =   np.pi/180.0

    #Get reference Xr, Yr, Zr..
    Xr, Yr, Zr  =   WGS842ECEF(lat0, lon0, alt0)
    #Get ECEF coordinates of the measurment required..
    Xp, Yp, Zp  =   WGS842ECEF(lat, lon, alt)
    #Get the difference#
    Xd          =   (Xp-Xr)
    Yd          =   (Yp-Yr)
    Zd          =   (Zp-Zr)

    #ECEF2ENU - https://en.wikipedia.org/wiki/Geographic_coordinate_conversion#From_ECEF_to_ENU
    E           =   -1*np.sin(lon*d2r)*Xd+np.cos(lon*d2r)*Yd
    N           =   -1*np.sin(lat*d2r)*np.cos(lon*d2r)*Xd-np.sin(lat*d2r)*np.sin(lon*d2r)*Yd+np.cos(lat*d2r)*Zd
    U           =   np.cos(lat*d2r)*np.cos(lon*d2r)*Xd+np.cos(lat*d2r)*np.sin(lon*d2r)*Yd+np.sin(lat*d2r)*Zd
    
    return E, N, U


cpdef WGS842ECEF(double lat, double lon, double alt):
    '''
    Function to cover WGS84 to ECEF, with lat0, lon0, alt0 as the centre.
    All inputs should be in degrees.
    
    Parameters
    ----------
        lat: `double`
        lon: `double`
        alt: `double`
    
    Returns
    -------
        [X,Y,Z]: `list`
            X, Y, Z (ECEF) coordinates

    Examples
    --------
        WGS842ECEF(double lat, double lon, double alt)
    '''
    cdef double X=0
    cdef double Y=0
    cdef double Z=0
    cdef double d2r =   np.pi/180.0

    cdef double a   =6378137.0
    cdef double b   =6356752.3142
    
    cdef double N_phi=Nphi(lat)

    X   =   (N_phi+alt)*np.cos(lat*d2r)*np.cos(lon*d2r)
    Y   =   (N_phi+alt)*np.cos(lat*d2r)*np.sin(lon*d2r)
    Z   =   (b**2/a**2 * N_phi+ alt)*np.sin(lat*d2r)

    return [X, Y, Z]


cdef Nphi(lat):
    '''
    Function to compute corrected Earth's radius with WSG-84 model.
    
    Earth's equatorial and polar radius from https://en.wikipedia.org/wiki/Earth_radius
    IERS	WGS-84 ellipsoid, semi-major axis (a)	6378137.0
    IERS	WGS-84 ellipsoid, semi-minor axis (b)	6356752.3142	[6]
    IERS	WGS-84 first eccentricity squared (e2)	0.00669437999014

    Parameters
    ----------
        lat: `float`
            Latitude

    Returns
    -------
        `float`
        Returns corrected Earth's radius at latitude `lat`
    '''
    cdef double N   =0
    cdef double a   =6378137.0
    cdef double b   =6356752.3142
    cdef double e2  =0.00669437999014
    N               =a/np.sqrt(1-e2*(np.sin(lat*np.pi/180))**2)
    return N


cpdef AltAz2HARA(float el, float az, float phi, float lon,
                    float second, float minute, float hour, float day, float month, float year):    
    '''
    Module to calculate RA and Dec, HA and Dec from AltAz

    .. code-block:: text

        --------------Convention for HA-----------
                        N
                   +ve     -ve

                W               E

                   +ve     -ve
                        S
        -----------------------------------------
    Parameters
    ----------
        el: `float`
        az: `float`
            Azimuthal angle
        phi: `float`
            Phi angle
        lon: `float`
            Longitude of the observatory
        second: `float`
        minute: `float`
        hour: `float`
        day: `float`
        month: `float`
        year: `float`

    Returns
    -------
        RA and Dec: `~numpy.array`
        HA and Dec: `~numpy.array`

    
    Examples
    --------
    All inputs in degrees

    AltAz2HARA(float elevaltion, float azimuth, float phi/lat, float longitude,
        float sec, float minu, float hour, float day. float month, float year)

    >>> AltAz2HARA(0, 45, 12.97, 77.58, 0, 0, 0, 1, 12, 2021) 
    '''
    
    cdef float sa = np.sin(az*np.pi/180);
    cdef float ca = np.cos(az*np.pi/180);
    cdef float se = np.sin(el*np.pi/180);
    cdef float ce = np.cos(el*np.pi/180);
    cdef float sp = np.sin(phi*np.pi/180);
    cdef float cp = np.cos(phi*np.pi/180);


    x = -ca * ce * sp + se * cp;
    y = -sa * ce;
    z = ca * ce * cp + se * sp;

    r = np.sqrt(x * x + y * y);
    if (r == 0.0):
        ha = 0.;
    else:
        ha = math.atan2(y, x);

    dec =   np.arcsin(sp*se+cp*ce*ca)#math.atan2(z, r);
    cdel=   np.cos(dec) 
    sdel=   np.sin(dec)
    ha  =   np.arccos((se-sp*sdel)/(cp*cdel))   
    
    if(ha < 0):
        ha = ha+24*15*np.pi/180
    cdef float RA       =       selflst(second, minute, hour, day, month, year) - ha/15.0 * 180/np.pi
    return np.array([RA, dec*180/np.pi]), np.array([ha/15.0*180/np.pi, dec*180/np.pi] )


cdef HA2RA(HA, second, minute, hour, day, month, year):
    """
    Function to convert Hour-angle to Right Ascension
    Returns right ascension using the formula RA = LST - HA

    Parameters
    ----------
        HA: `float`
            Hour-angle
        second: `float`
        minute: `float`
        hour: `float`
        day: `float`
        month: `float`
        year: `float`

    Returns
    -------
        RA: `float`
            Right ascension
    """
    LST    =   selflst(second, minute, hour, day, month, year)
    RA     =   LST-HA
    return RA


cpdef map_RA_Dec_to_file(fil, pointing_file, gps, LO, Tile_1, Tile_2):
    '''
    .. code-block:: text

        Pointing file structure:
            chXX_NAME_YYYYMMDD_HHMMSS.mbr,   Az,  za
            .                               .   .
            .                               .   .
            .                               .   .
            .                               .   .
            .                               .   .

    Parameters
    ----------
        fil: `str`
            File name
        pointing_file: `str`
            Pointing file
        gps: `float`
            GPS Count
        Tile_1: `int`
            Tile number 1 (0-N convention)
        Tile_2: `int`
            Tile number 2 (0-N convention)

    Note
    ----
        Without timestamping this could result in inaccurate resuts,
        i.e the gps count has to be from the nearest 12th hour.  
    '''

    cdef float Alt
    cdef float Az
    cdef int noon   =   0
    cdef int fact   =   0
    cdef float hour =   0 
    file_list  =   open(pointing_file, 'r').readlines()
    for i in range(len(file_list)):
        fil_temp    =   file_list[i].split(' ')
        try:
            #Checking if path is absolute
            file_temp1   =   fil_temp[0].split('/')[-1]
        except:
            #Path is direct
            file_temp1   =   fil_temp[0]
        if(file_temp1 == fil):
            if(int((file_temp1.split('_')[-2])[0:2]) >= 12 ):
                hour    =   12
                print((file_temp1.split('_')[-2])[0:2])
            Az  =   float(fil_temp[1])
            Alt =   float(fil_temp[2])
            print('ZA..'+str(Alt)+'\tAzimuth..'+str(Az))
            fact=   1

    if(fact==0):
        print('No file match found..check pointing file')
        return -1
    #cdef list pointing  =   []
    cdef float hour1    =   gps/3600.0
    hour     =   hour+math.modf(hour1)[-1]
    cdef float minu1    =   (hour1%1) * 60
    cdef float minu     =   math.modf(minu1)[-1]
    cdef float sec      =   (minu1%1) * 60
    
    #cdef float dat      =   0
    #cdef float month    =   0
    #cdef float year     =   0
    
    print('Time Derived from GPS..'+str(hour)+'H:'+str(minu)+'M:'+str(sec)+'S') 
    #AltAz2HARA(float el, float az, float phi, float lon, float second, float minute, float hour, float day, float month, float year)
    #float sec, float minu, float hour, float day, float month, float year, float RA, float dec, float avg, float del_t, float time, unsigned int T1, unsigned int T2

    series, sec_junk, minu_junk, hour_junk, dat, month, year               =       _header_Fring_cy.extract_time(file_temp1)
    pointing =  AltAz2HARA((90.0-Alt), Az , LatLong.Lat_local, LatLong.Long_local, sec, minu, hour, float(dat), float(month), float(year))#, float(file_list[i].split(',')[-1]), float(file_list[i].split(',')[-1]))[0]
    delay =  Cal_time_onhold_v1(sec, minu, hour, float(dat), float(month), float(year), pointing[0][0], pointing[0][1], 6000, 0.1, 40, Tile_1, Tile_2)
    return delay, pointing

cpdef correct_time(float second, float minute, float hour, float day, float month, float year):
    """
    Parameters
    ----------
        second: `float`
            Seconds
        minute: `float` 
            Minutes
        hour: `float`
            Hours
        day: `float`
            Days
        month: `float`
            Month(1-12)
        year: `float`
            Year
    
    Returns
    -------
        secu: `float`
        mint: `float`
        hourt: `float`
        dayt: `float`
    """
    cdef double mint    =   minute
    cdef double secu    =   second
    cdef double hourt   =   hour
    cdef double dayt    =   day
    if( secu > 59):
        mint    =   mint+1
        secu    =   0   
    if(mint>59):
        mint    =   0   
        hourt   =   hourt+1
    if(hourt > 23):
        hourt   =   0   
        dayt    =   dayt+1
    return secu, mint, hourt, dayt

cpdef Equ2local(float RA, float Dec, float phi, float lon,
        float second, float minute, float hour, float day, float month, float year):
    '''
    Function to covert Equatorial coordinates to Local Coordinates
    Returns local coordinates, Azimuth, Altitude, Hour-angle

    Parameters
    ----------
        RA: `float`
            Right Ascension
        Dec: `float`
            Declination
        phi: `float`
        lon: `float`
            Longitude
        second: `float`
            Seconds
        minute: `float`
            Minutes
        hour: `float`
            Hours
        day: `float`
            Days
        month: `float`
            Month(1-12)
        year: `float`
            Year

    Returns
    -------
        az: `float`
            Azimuth
        alt: `float`
            Altitude
        ha: `float`
            Hour-angle

    Examples
    --------
        >>> Equ2local(RA, Dec, phi, lon, second, minute, hour, day, month, year) 
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
    lst                     =                   selflst(second, minute, hour, day, month, year)
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

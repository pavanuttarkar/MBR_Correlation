#Cython header#
import _header_gps_cy

import numpy as np
cimport numpy as np
import math
import os


cpdef extract_time(fil):
    '''
    Function to extract time tag from the file name, this is less reliable, as the file names
    are derived from the computer NTP (hopefully!) protocol.

    Parameters
    ----------
    fil : `str`
        file name as string

    Returns
    -------
    series : unsigned int
        series of the file.
    second : unsigned int
        second value form the file
    minute : unsigned int
        minute value form the file
    hour   : unsigned int
        hour value form the file
    day    : unsigned int
        day value form the file
    month  : unsigned int
        month value form the file
    year   : unsigned int
        year value form the file

    '''

    cdef str series = (fil[-7:-4])
    cdef unsigned int  dat = int((fil.split("_")[-3])[6:8])
    cdef unsigned int month = int((fil.split("_")[-3])[4:6])
    cdef unsigned int year = int((fil.split("_")[-3])[0:4])
    cdef unsigned int hour = int((fil.split("_")[-2])[0:2])
    cdef unsigned int minu = int((fil.split("_")[-2])[2:4])
    cdef unsigned int sec = int((fil.split("_")[-2])[4:6])
    return series, sec, minu, hour, dat, month, year

cpdef extract_series(fil):
    '''
    Function to extract series tag from the file name

    Parameters
    ----------
    fil : `str`
        file name as string

    Returns
    -------
    series : `str`
        series of the file
    '''
    cdef str series = (fil[-7:-4])
    return series


cpdef gps_sync(file_name, file_name1, fact):
    '''
    Function to calculate synchronization solution for the given 000 series data set
    if extracts using the synchronization series Info files generated using the
    :meth:`~api._header_gps_cy.super_sync_generator` function. This can oly be used for 000 series file and
    non-000 series files use a different method to extract the synchronization
    solution, using a combination of packet loss extimation and previous sync
    solutions.
            
    Parameters
    ----------
    1st filename : `str`
        file name as string
    2nd filename : `str`
        file name as string
    sync factor  : `int`
        if synchronization is required - 1
        if no synchronization solution - 0

    Returns
    -------
    Memfactor : 1-D array of the order 1X2
        Memfactor[0] - corresponds to the Memfactor jump for 1st filename.
        Memfactor[1] - corresponds to the Memfactor jump for 2st filename. 

    '''
    cdef list SynFile = []
    cdef np.ndarray ext = np.zeros((2), dtype=float)
    cdef list f512 = []
    cdef list args = []
    cdef int ctr = 0
    cdef int numprocess = len(ext)
    cdef list NoPacktoSkip = []
    cdef list barebinary = []
    cdef list GPScount = []
    cdef list PackCount = []
    cdef list Memfactor = []
    cdef int i = 0
    cdef list GPS_st = []
    if(fact == str(1)):
        print('Synchcronization in progress>>>>>>>>')

        ext[0] = (os.system('ls -lrth SAMPLING_INFO/Info_on_straight_line'+str(file_name.split('/')[-1])))
        ext[1] = (os.system('ls -lrth SAMPLING_INFO/Info_on_straight_line'+str(file_name1.split('/')[-1])))
        if(ext[0] == 512):
            _header_gps_cy.super_sync_generator(file_name)
            SynFile.append(open('SAMPLING_INFO/Info_on_straight_line'+str(file_name.split('/')[-1])))
        else:
            SynFile.append(open('SAMPLING_INFO/Info_on_straight_line'+str(file_name.split('/')[-1])))

        if(ext[1] == 512):
            _header_gps_cy.super_sync_generator(file_name1)
            SynFile.append(open('SAMPLING_INFO/Info_on_straight_line'+str(file_name1.split('/')[-1])))
        else:
            SynFile.append(open('SAMPLING_INFO/Info_on_straight_line'+str(file_name1.split('/')[-1])))

        for i in range(2):
            GPScount.append(int(((SynFile[i].readline().split(" ")[-1]).split('\n')[0]).split('.')[0]))
            PackCount.append(int(SynFile[i].readline().split('.')[0]))
            NoPacktoSkip.append(SynFile[i].readline())
            junk    =   SynFile[i].readline()
            GPS_st.append(SynFile[i].readline().split(','))
            print('GPS count..file no.' +str(i)+'..'+str(GPScount[i]))
        x = max(GPScount)+1
        for i in range(2):
            print(NoPacktoSkip[i])
            print(PackCount[i])
            a = math.modf(eval(NoPacktoSkip[i]))
            print(a)
            skipint = int(a[1])
            skipfloat = float(a[0])
            Memfactor.append(((skipint+skipfloat)*512.0))# - PackCount[i])*512.0)
        Memfactor.append(GPS_st[0])
        Memfactor.append(GPS_st[1])
        #print('Skiping --'+str(Memfactor*1056))

        return Memfactor, x
    else:
        Memfactor = [0, 0]
        return Memfactor


cpdef extract_XY_vals(file_name, file_name1):
    '''
    Function to extract polarization and DSP numbers from the file name

    Parameters
    ----------
    1st filename : string
        file name as string

    2nd filename : string
        file name as string

    Returns
    -------
    chX1 : `str`
        DSP number associated with 1st file X-Pol.
    chX2 : `str`
        DSP number associated with 2st file X-Pol. 
    chY1 : `str`
        DSP number associated with 1st file Y-Pol.
    chY2 : `str`
        DSP number associated with 2st file Y-Pol. 
    fil1 : list of `str`
        1st filename string split to get filename parameters
    fil2 : list of `str`
        2nd filename string split to get filename parameters
    '''

    cdef list fil1     =   file_name.split('_')
    cdef list fil2     =   file_name1.split('_')
    cdef str chX1      =   'X'+str((fil1[-5])[-1])
    cdef str chX2      =   'X'+str((fil2[-5])[-1])
    cdef str chY1      =   'Y'+str((fil1[-5])[-1])
    cdef str chY2      =   'Y'+str((fil2[-5])[-1])
    return chX1, chX2, chY1, chY2, fil1, fil2


cpdef gps_slope(file_name):
    '''
    Function to get slope from the 000th file name of the series.

    Prarameters
    -----------
    filename  : string
        file name as string

    Returns
    -------
    slope : `float`
        slope values of the straight line fit, calculated for GPS blip vs packet transitions.  

    .. warning::

        Gives error if there is no SMAPLING_INFO file containing the
        Info_on_straightline* file.
    '''
    cdef str file_name_1
    try:
        file_name_1 = file_name.split('/')[-1]
    except:
        file_name_1 = file_name
    file_name_1 = file_name_1[:-7]+'000.mbr'
    SynFile = open('SAMPLING_INFO/Info_on_straight_line'+str(file_name_1.split('/')[-1]))

    junk = SynFile.readline()
    junk = SynFile.readline()
    NoPacktoSkip = float(SynFile.readline().split('*x')[0])
    return NoPacktoSkip

# data_functions.py
import numpy as np
from spacepy import coordinates as coords
from spacepy.time import Ticktock
import aacgmv2
from datetime import datetime

### check if sliceAB and sliceAB2 can be combined.
def sliceAB(data, subunit, n):
    ''' Slices 1 dimensional data into A and B based on whether subunit is 
        equal to n (goes into A) or not (goes into B). This function is used 
        for slicing PC1 and TIME data for MEPD-A and MEPD-B. '''
    A = []
    B = []
    for i in range(len(data)):
        if subunit[i] == n:
            A.append(data[i])
        else:
            B.append(data[i])
    return A, B

def sliceAB2(data, subunit, n):
    ''' Slices 2 dimensional data into A and B based on whether subunit is 
        equal to n (goes into A) or not (goes into B). This function is used 
        for slicing MEPD detector data for MEPD-A and MEPD-B. '''
    A = []
    B = []
    for i in range(len(data)):
        if subunit[i] == n:
            A.append(data[i, :])
        else:
            B.append(data[i, :])
    return A, B

def B_avg(B):
    ''' Calculates average value (magnitude) for given isss magnetic field 
        data. '''
    avg = []
    for i in range(len(B)):
        avg.append(np.sqrt((B[i, 4])**2 + (B[i, 5])**2 + (B[i, 6])**2))
    return avg

def closest(arr, value):
    ''' Identifies an element in arr which has the closest value to value. '''
    c = 0
    for i in range(len(arr)):
        if (abs(arr[i] - value) < abs(arr[c] - value)):
            c = i
    return c

def geo_lat(alt, start_time, conv_module):
    ''' Calculates the geomagnetic lattitude data at the given alt (altitude) 
        and start_time. This function uses either Spacepy or aacgmv2 to convert 
        between geographic and geomagnetic coordinates. The moduel used to 
        convert coordinates can be selected by setting the conv_module. Make
        sure to use all lowercase for converting module. Keep in mind SpacePy 
        for some reason is not able to use 2020 data. '''
    arr = np.zeros((181, 360))
    geo_lat = np.zeros((5, 360))
    for j in range(360):
        for i in range(181):
            # Altitude is given in meters but Spacepy uses kilometers.
            coordinates = coords.Coords([alt / 1000, i - 90, j - 180], \
                                'GEO', 'sph')
            # For some reason, 2020 data could not be used.
            if conv_module == 'spacepy':
                coordinates.ticks = Ticktock(['2019-07-17T17:51:15'], 'ISO')
                arr[i][j] = coordinates.convert('MAG', 'sph').lati
            elif conv_module == 'aacgmv2':
                arr[i][j] = (np.array(aacgmv2.get_aacgm_coord(i - 90, j - 180,\
                    int(alt / 1000), start_time)))[0]
            else:
                print("Error: coordinate conversion module is invalid.\n\
                    Please choose between:\n  spacepy\n  aacgmv2\n\
                    Please use all lowercase.")
                exit()
    for j in range(360):
        for i in range(5):
            geo_lat[i, j] = closest(arr[:, j], 30 * i - 60) - 90
    return geo_lat

def slice_tel(tel):
    ''' Slices HEPD telescope data into proton and electron data. '''
    tel0, tel1, tel2 = tel[0], tel[1], tel[2]
    proton0 = tel0[:, 17:21]
    proton1 = tel1[:, 17:21]
    proton2 = tel2[:, 17:21]
    proton = [proton0, proton1, proton2]
    electron0 = tel0[:, 2:13]
    electron1 = tel1[:, 2:13]
    electron2 = tel2[:, 2:13]
    electron = [electron0, electron1, electron2]
    return proton, electron

def start_end_time(data):
    ''' Gets the start and end time and converts the UNIX datetime to UTC 
        datetime. '''
    start_time = datetime.fromtimestamp(data.time[0])
    end_time = datetime.fromtimestamp(data.time[-1])
    return start_time, end_time

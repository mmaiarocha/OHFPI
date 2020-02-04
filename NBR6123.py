# -*- coding: utf-8 -*-

import sys
import gzip    as gz
import pickle  as pk
import numpy   as np

import scipy.interpolate as si

#=============================================================================

def profile(Z, Tav=3, Cat=2):
    """ Z:     vector with heights above ground in meters
        Tav:   averaging time in seconds
        Cat:   surface roughness category (1, 2, 3, 4, or 5)
    """

# 1. Table of parameters according to NBR6123

    Ti  =  np.array( 
          (    3,    5,   10,   15,   20,   30,
              45,   60,  120,  300,  600, 3600))

    bi  =  np.array(
          ((1.10, 1.11, 1.12, 1.13, 1.14, 1.15, 
            1.16, 1.17, 1.19, 1.21, 1.23, 1.25),
           (1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 
            1.00, 1.00, 1.00, 1.00, 1.00, 1.00),
           (0.94, 0.94, 0.93, 0.92, 0.92, 0.91, 
            0.90, 0.90, 0.89, 0.87, 0.86, 0.85),
           (0.86, 0.85, 0.84, 0.83, 0.83, 0.82, 
            0.80, 0.79, 0.76, 0.73, 0.71, 0.68),
           (0.74, 0.73, 0.71, 0.70, 0.69, 0.67, 
            0.64, 0.62, 0.58, 0.53, 0.50, 0.44)))

    pi  =  np.array(
          ((0.06, 0.065,0.07, 0.075,0.075,0.08, 
            0.085,0.085,0.09, 0.095,0.095,0.10),
           (0.085,0.09, 0.10, 0.105,0.11, 0.115,
            0.12, 0.125,0.135,0.145,0.15, 0.16),
           (0.10, 0.105,0.115,0.125,0.13, 0.14, 
            0.145,0.15, 0.16, 0.175,0.185,0.20),
           (0.12, 0.125,0.135,0.145,0.15, 0.16, 
            0.17, 0.175,0.195,0.215,0.23, 0.25),
           (0.15, 0.16, 0.175,0.185,0.19, 0.205,
            0.22, 0.23, 0.255,0.285,0.31, 0.35)))

    Fi  =  np.array( 
          ( 1.00, 0.98, 0.95, 0.93, 0.90, 0.87, 
            0.84, 0.82, 0.77, 0.72, 0.69, 0.65))

    z0i =  np.array([ 0.005, 0.070, 0.300, 1.000, 2.500])
    cai =  np.array([ 2.800, 6.500, 10.50, 22.60, 52.70])
    bti =  np.array([ 6.500, 6.000, 5.250, 4.850, 4.000])

# 2. Interpolate parameters according to category and averaging time

    if (Cat < 0.0) and (Cat > 5.0):
        sys.exit('Terrain category out of range.')

    if (Tav < 3) or (Tav > 3600):
        sys.exit('Averaging gust duration out of range.')

    C1 = np.int(Cat//1 - 1)
    C2 = np.int(C1 + 1)
    
    b  = np.interp(Tav,Ti,bi[C1,:]) 
    p  = np.interp(Tav,Ti,pi[C1,:]) 
    Fr = np.interp(Tav,Ti,Fi) 

    z0 = z0i[C1] 
    ca = cai[C1] 
    bt = bti[C1] 

# 3. Calculate the wind profile
    
    Z0 = np.array(Z, dtype=float)
    S2 = np.zeros(Z0.shape) 
    
    if (Cat < 2.0):
        Z0[Z0 > 250] = 250.
  
    S2[Z0 > 0] = b*Fr*(Z0[Z0 > 0]/10)**p 
    
# 4. If category is fractionary, interpolate two categories (by recursion!!!)
    
    if  Cat%C2 > 0.1:
        
        S22, b2, p2, Fr2, z02, bt2, ca2 = profile(Z, Tav, C2+1)
        
        S2 = (S2 + S22)/2
        b  = (b  + b2 )/2
        p  = (p  + p2 )/2
        z0 = (z0 + z02)/2
        bt = (bt + bt2)/2
        ca = (ca + ca2)/2
    
    return S2, b, p, Fr, z0, bt, ca


#=============================================================================

def drag(Bx, By, Hz, case='high'):
    """ 
    Bx:    base dimension parallel to x axis
    By:    base dimension parallel to y axis
    Hz:    building height
    case: 'low' or 'high' turbulence
    """

    BxBy  =  Bx/By if (Bx/By <  4) else  4
    ByBx  =  By/Bx if (By/Bx <  4) else  4
    HzBy  =  Hz/By if (Hz/By < 10) else 10
    HzBx  =  Hz/Bx if (Hz/Bx < 10) else 10

    with gz.GzipFile('NBR6123.gz', 'rb') as target:
        low, high = pk.load(target)
        
    if (case.lower() == 'low'):
        X = low   
    elif (case.lower() == 'high'):
        X = high
    else:
        sys.exit('Turbulence case must be either "low" or "high".')    
    
    Ca =  si.Rbf(X[:,0],X[:,1],X[:,2],function='thin_plate')
    
    Cx =  Ca(ByBx, HzBy)
    Cy =  Ca(BxBy, HzBx)
    
    return Cx, Cy


#=============================================================================



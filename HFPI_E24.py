# -*- coding: utf-8 -*-

import numpy             as np
import pandas            as pd
import matplotlib.pyplot as plt

import URP
import HFPI
import HFPI0
import NBR6123

HFPI.ALIM = 20.     # limit for acceleration plots (mG)
HFPI.DLIM = 50.     # limit for displacement plots (mm)
HFPI.RLIM = 0.5     # limit for rotation plots (degrees)
HFPI.MMAX = 8       # maximum number of modes to be considered

plt.close('all')

#==============================================================================
# 1. Reads all cases to be iterated and define analysis parameters
#------------------------------------------------------------------------------
All_Complete   = pd.read_excel('AllCases.xlsx', sheetname='Complete')
All_Structures = pd.read_excel('AllCases.xlsx', sheetname='Structures')

# overall damping ratio of critical
zt  =  0.01
zts =  '1%'

# basic velocity (also defines time scale)
V0  =  39.2                  

# topographic factor (NBR6123)
S1  =  1.0                  

# statistical factor ("m" years, 63%)
m   =  1
S3  =  0.541*((-(np.log(1 - 0.632))/m)**(-0.157))

# Averaging time for mean and peak analysis
Tav =  600

#==============================================================================
# 2. Iterate over all cases
#------------------------------------------------------------------------------
for kc, c in All_Complete.iterrows():

    if (kc != 54): continue

    dirname   =  c['dirname'] + '/'
    prefix    =  c['prefix']
    structure =  All_Structures.loc[c['strt']] 
    graph     = (dirname, prefix, m, zts)

# length scale
    lL      =  structure['lL']

# Tap offset from wind tunnel reference
    offT    =  structure[['dxT', 'dyT', 'dzT', 'rzT']].values.astype('float')
    offT[3] =  np.pi*offT[3]/180

# TQS offset from wind tunnel reference
    offN    =  structure[['dxN', 'dyN', 'dzN', 'rzN']].values.astype('float')
    offN[3] =  np.pi*offN[3]/180
    
# building dimensions
    Bx     =  structure['Bx']
    By     =  structure['By']
    Hz     =  structure['Hz']
    H0     =  structure['H0']
    H1     =  structure['H1']

# Get Z coordinates to provide static equivalent loads
    CotaZ  =  pd.read_excel(dirname+'E24_FloorCoord.xlsx')
    ZS     =  CotaZ['ZS'].values.astype('float')

#==============================================================================
# 3. Read case batch file
#------------------------------------------------------------------------------
    batch  =  HFPI.read_batch(dirname, prefix)

#==============================================================================
# 4. Define structural model and reset collectors
#------------------------------------------------------------------------------
    coorfile =  dirname + structure['coor']
    massfile =  dirname + structure['mass']
    freqfile =  dirname + structure['freq']
    modefile =  dirname + structure['mode']
    str_file =  dirname + structure['strt']

    NN, Ni, XN, MN, fk, QN, ortho  =  \
    HFPI0.read_TQS(coorfile, massfile, freqfile, modefile, str_file)
              
    XN, QN =  HFPI0.offset_TQS(XN, QN, offN)
        
    NS, Si, XS, MS, fk, QS, KS =  \
    HFPI0.TQS_to_3dof(NN, Ni, XN, MN, fk, QN, ZS)

    MSi = MS.reshape(NS,3)

    FX    =  np.zeros((NS, len(batch)))
    FY    =  np.zeros((NS, len(batch)))
    MZ    =  np.zeros((NS, len(batch)))
     
    u_max =  np.zeros((len(batch),3))
    a_max =  np.zeros((len(batch),3))

    Pt    =  np.zeros((len(batch),3))
    Ce    =  np.zeros((len(batch),3))
    
    fS    =  np.zeros((len(batch),3))
    ST    =  np.zeros((len(batch),3))
    Vc    =  np.zeros((len(batch),3))

#==============================================================================
# 5. Iterate over all wind directions
#------------------------------------------------------------------------------
    for kb, b  in batch.iterrows():

#       if (kb != 8): continue

        wnd, Cat, Vm, dPa, Td =  HFPI.unpack_batch(b)
        
        wgraph = (graph, wnd)
        
        if (Cat > 3): Cat = 4         # Corrects category

        S2, *rest = NBR6123.profile(Hz, Tav=Tav, Cat=Cat)

        Vp  =  S1*S2*S3*V0            # wind velocity at building top
        p0  =  0.613*Vp*Vp            # reference dynamic pressure
        lT  =  lL/(Vp/Vm)             # time scale by Strouhal similarity
        Td  =  lT*Td                  # apply time scale
        
#------------------------------------------------------------------------------
# 5.1 Open pressure data for wind direction and interpolate modal shapes
#------------------------------------------------------------------------------
        master, Cp     =  HFPI.get_taps(b)
        tap, A, XT, CT =  HFPI.unpack_master(master)

        n   =  Cp.shape[0]
        fs  =  Td/n
        Cp  =  p0*Cp

        XT, CT =  HFPI.offset_taps(XT, CT, offT)     ########
        master[['cx','cy','cz']] = CT
    
        if (kb == 0):
            QT =  HFPI0.interp_TQS(XT, XN, QN)

#------------------------------------------------------------------------------
# 5.2 Calculate modal loads
#------------------------------------------------------------------------------
        Pk, mPk, sPk, fPk =  \
        HFPI0.mload_TQS(wgraph, master, Cp, XT, QT, Td, plot=False)

#------------------------------------------------------------------------------
# 5.3 Calculate modal responses (important: modal masses are all unity)
#------------------------------------------------------------------------------
        uk, muk, suk, fuk =  \
        HFPI.modal_response(wgraph, fk, zt, Td, Pk, plot=False)

#------------------------------------------------------------------------------
# 5.4 Find the moment of maximum elastic energy
#------------------------------------------------------------------------------
        Xe =  0.0
        Ye =  0.0
        Ze =  ZS.max()
        
        De =  np.sqrt((XN[:,0] - Xe)**2 +
                      (XN[:,1] - Ye)**2 +
                      (XN[:,2] - Ze)**2)
        
        ne =  np.argmin(De)   # control node for biasing Fx/Fy
        
        i_Emax, Emax, ukmax =  \
        HFPI.maximum_energy(wgraph, fk, uk, fs, Tav, QN, ne, plot=False)

#------------------------------------------------------------------------------
# 5.5 Calculate and diplay displacements and accelerations at a given floor
#------------------------------------------------------------------------------
        Xc =   0. #Bx/2.
        Yc =   0. #By/2.
        Zc =   Ze
        
        Dc =  np.sqrt((XN[:,0] - Xc)**2 +
                      (XN[:,1] - Yc)**2 +
                      (XN[:,2] - Zc)**2)
        
        nc =  np.argmin(Dc)   # control node for displacements to present

# ATTENTION: uT[3,:] and aT[3,:] are the instantaneous resultants
        
        uT, aT =  \
        HFPI0.top_TQS(wgraph, i_Emax, nc, fk, uk, QN, Td, Xd=True, plot=False)

        for i in range(3):
            nmax, Xmax, u_max[kb, i] = URP.splitmax(uT[i,:], fs, Tav)
            nmax, Xmax, a_max[kb, i] = URP.splitmax(aT[i,:], fs, Tav)

#------------------------------------------------------------------------------
# 5.6 Defines equivalent static load and convert to 'tonf'
#------------------------------------------------------------------------------
        XTC = np.array([0.00, 0.00, 0.00], dtype='float')  # Center of torsion

        Pe  = HFPI0.stateq_TQS(fk, ukmax, XN, MN, QN, KS, XTC).reshape(NS,3)

        FX[:,kb] = Pe[:,0]/9810
        FY[:,kb] = Pe[:,1]/9810
        MZ[:,kb] = Pe[:,2]/9810

#------------------------------------------------------------------------------
# 5.7 Total load at Z = 0, and equivalent global coefficients
#------------------------------------------------------------------------------
        Pt[kb,:] = np.sum(Pe, axis=0)

        Ce[kb,0] = Pt[kb,0]/By/Hz/p0
        Ce[kb,1] = Pt[kb,1]/Bx/Hz/p0
        Ce[kb,2] = Pt[kb,2]/By/Bx/Hz/p0
        
#==============================================================================
# 6. Plot a comparison of aerodynamic coefficients
#------------------------------------------------------------------------------
    alpha, Cm, Cs = \
    HFPI.aero_coeffs(graph, batch, structure, Ce, plot=True)

    Cmdf = pd.DataFrame(data   = Cm,
                        index  = alpha,
                        columns= ['Cmx', 'Cmy', 'Cmt'])
    
    Csdf = pd.DataFrame(data   = Cs,
                        index  = alpha,
                        columns= ['Csx', 'Csy', 'Cst'])
    
    Cedf = pd.DataFrame(data   = Ce,
                        index  = alpha,
                        columns= ['Cex', 'Cey', 'Cet'])

    fst = prefix + '_{0:0>3}Y_{1}_AeroCoeffs.xlsx'
    fst = fst.format(m, zts)

    excelCA = pd.ExcelWriter(fst)
    
    Cmdf.to_excel(excelCA,'Cm')
    Csdf.to_excel(excelCA,'Cs')
    Cedf.to_excel(excelCA,'Ce')
    
    excelCA.save()
    
#==============================================================================
# 7. Plot maximum displacements and accelerations at building top
#------------------------------------------------------------------------------
    HFPI.peak_response(graph, alpha, 1000*u_max, 1000*a_max/9.81)
''' 
#==============================================================================
# 8. Save static equivalent loads
#------------------------------------------------------------------------------
    FXdf = pd.DataFrame(data   = FX,
                        index  = ZS,
                        columns= alpha)
                        
    FYdf = pd.DataFrame(data   = FY,
                        index  = ZS,
                        columns= alpha)
                        
    MZdf = pd.DataFrame(data   = MZ,
                        index  = ZS,
                        columns= alpha)
        
    fst = prefix + '_{0:0>3}Y_{1}_StaticEquiv.xlsx'
    fst = fst.format(m, zts)

    excelFE = pd.ExcelWriter(fst)
    
    FXdf.to_excel(excelFE,'FX')
    FYdf.to_excel(excelFE,'FY')
    MZdf.to_excel(excelFE,'MZ')
    
    excelFE.save()

#==============================================================================
# 9. Compare with statically defined loads
#------------------------------------------------------------------------------
    FXst, FYst, MZst = \
    HFPI0.compare_loads(graph, Hz, ZS, alpha, FXdf, FYdf, MZdf)
 
#==============================================================================
'''
# -*- coding: utf-8 -*-

import sys
import gzip    as gz
import pickle  as pk
import numpy   as np
import pandas  as pd

import matplotlib.pyplot as     plt
from   scipy.interpolate import griddata    # for interpolation

import URP
import HFPI

#=============================================================================

def scanivalve_process(dirname, prefix, avg=1, per=64*32):

# 1. Define data type for each column of files to be read           

# 1.1 Batch file
    typ1  = ['scan',   # name of Scanivalve file without extension
             'master', # master table file to be used
             'Cat',    # far field roughness category 
             'block',  # wind tunnel blockage factor
             'Vm',     # wind tunnel mean speed
             'k0']     # wind tunnel k0 factor, referred to model top
    
# 1.2 Processed batch file
    typ2  = ['wnd',    # wind direction
             'file',   # file with processed pressure coefficients
             'Cat',    # far field roughness category 
             'Vm',     # wind tunnel mean speed
             'dPa',    # reference wind tunnel pressure
             'Td']     # total pressure acquisition time

# 1.3 Scanivalve file
    typ3  = ['set',    # ignored
             'frm',    # acquisition frame
             'mdu',    # Scanivalve module (1 to 6)
             'chn',    # Scanivalve channel (1 to 64)
             'prs']    # pressure reading

# 1.4 Master table
    typ4  = ['tap',    # tap number
             'md1',    # Scanivalve first module number
             'ch1',    # first module channel
             'md2',    # Scanivalve second module number
             'ch2',    # second module channel
             'zone',   # pressure zone on model surface
             'phas',   # acquisition phase for multiple takes
             'A',      # surface area assigned to tap
             'X',      # |
             'Y',      #  > tap coordinates in model reference system
             'Z',      # |
             'cx',     # |
             'cy',     #  > direction cosines pointing outwards surface
             'cz']     # |

# 2. Otherwise, preprocess Scanivalve files
        
    print('Processing "{0}" ... '.format(dirname))
    
    batchpickle =  dirname + prefix + '.pk'
    batchfile   =  dirname + prefix + '_Batch.txt'
    batchdata   =  pd.read_table(batchfile, names=typ1, header=None)
    batch       =  pd.DataFrame(columns=typ2)

    for ib, b in batchdata.iterrows():
       
# 2.1 Read Scanivalve file

        rec  =  dirname + b.scan   + '.DAT'
        mst  =  dirname + b.master + '.txt'
        wnd  =  int(rec[-7:-4])               # wind direction from filename
#       wnd  =  int(rec[-11:-8])              # wind direction from filename

# CAARC
#       wnd  =  int(rec[-7:-4]) - 180
#       if wnd < 0: wnd += 360

# New York
#       wnd  =  int(rec[-7:-4]) - 45
#       if wnd < 0: wnd += 360
        
        print('  Direction = {0:0>3}deg: "{1}"... '.format(wnd, rec))

        data =  pd.read_table(rec, sep=' ', names=typ3, header=None)
        
        n    =  1 + int(data.iloc[-1,1])
        Td   =  n*(avg*per*1.e-6)
        
        mdu  =  data.mdu.values              # unicity gets troubled ... 
        chn  =  data.chn.values              # ... with non numpy arrays.
        
        ID   =  np.unique(100*mdu + chn)     # allocate and ...
        nID  =  len(ID)                      # ... get length.
        
        pr0  =  np.zeros((n,nID))            # unordered pressure array

# 2.2 Convert SCANIVALVE file to time series

        for i in range(nID):
            ID[i]    = 100*mdu[i] + chn[i]   # ensure correspondence
            pr0[:,i] = data.prs[i::nID]

# 2.3 Read master table and find indexing of pressure taps within records             

        table  = pd.read_table(mst, names=typ4, index_col=None, header=None)
        master = pd.DataFrame(columns=typ4)

        ID1 = []
        ID2 = []
        
        for it, t in table.iterrows():
            
            md1 = t.md1.astype('i')
            ch1 = t.ch1.astype('i')      
            md2 = t.md2.astype('i')
            ch2 = t.ch2.astype('i')
            
#           tap = t.tap.astype('i')
#           print(tap,md1,ch1,md2,ch2)
            
# CAARC
#           t.cx = -t.cx
#           t.cy = -t.cy
#           t.cz = -t.cz

            try:
                if t.zone not in [-1, -2]:
                    
                    ID1.append(list(ID).index(100*md1 + ch1))
                    ID2.append(list(ID).index(100*md2 + ch2))
                    
                    master = master.append(t, ignore_index=True)

                elif t.zone == -1:
                    IA1 = list(ID).index(100*md1 + ch1)
            
                elif t.zone == -2:
                    IA2 = list(ID).index(100*md1 + ch1)

            except ValueError:
                sys.exit('Channel not in Scanivalve file!')

# 2.4 Convert relative pressure to pressure coefficients

        A1  =  np.mean(pr0[:,IA1])        
        A2  =  np.mean(pr0[:,IA2])     
        dPA =  b.k0*(A1 - A2)
        Cp  =  pd.DataFrame()   # reset pressure coefficients
        
        for id1, id2, m in zip(ID1, ID2, master.index):
            Cp[m] = (pr0[:,id1] + pr0[:,id2])/2/dPA/b.block

# 2.5 Dump pressure data to a zipped pickle file

        swn = '{0:0=3}'.format(wnd)
        pck =  dirname + prefix + '/' + prefix + '_' + swn + '.gz'

        with gz.GzipFile(pck, 'wb') as target:
            pk.dump((master, Cp), target)

# 2.6 Append data to batch dataframe

        new_b = pd.Series([wnd, pck, b.Cat, b.Vm, dPA, Td], index=typ2)
        batch = batch.append(new_b, ignore_index=True)

# 3. Dump and return batch dataframe

    batch.wnd = batch.wnd.astype('int')
    
    with open(batchpickle, 'wb') as target:
        pk.dump(batch, target)

    print('... done!!!')
    return batch


#=============================================================================

def read_3dof(massfile, freqfile, modefile, str_file):
 
    try:
        with gz.GzipFile(str_file+'.gz', 'rb') as target:
            return pk.load(target)
    except:
        pass

# 1. Define data type for each column of files to be read           

# 1.1 Height and lumped masses
    typ1  = [('Si','i'),        # storey ID
             ('ZS','f'),        # vertical coordinate of storey center of mass
             ('MS','f'),        # storey mass 
             ('IS','f')]        # storey moment of inertia (vertical axis)

# 1.2 Frequencies
    typ2  = [('K','i'),         # mode ID
             ('fk','f')]        # modal frequency

# 1.3 Modal shapes
    typ3  = [('Si','i'),        # storey ID
             ('qX','f'),        # modal displacement in X
             ('qY','f'),        # modal displacement in Y
             ('qR','f')]        # modal displacement in Z (rotation)

# 2. Build up lumped mass matrices

    data =  np.genfromtxt(massfile+'.txt', dtype=typ1)

    Si   =  data['Si']
    NS   =  len(Si)
    ZS   =  data['ZS']
    X0   =  np.zeros(ZS.shape)

    XS   =  np.vstack((      X0,         X0,         ZS  )).T
    MS   =  np.vstack((data['MS'], data['MS'], data['IS'])).T.reshape((3*NS,))

# 3. Read modal frequencies (must be already sorted)

    data =  np.genfromtxt(freqfile+'.txt', dtype=typ2)
    fk   =  data['fk']
    fk   =  fk[0:HFPI.MMAX]

# 4. Read modal shapes (must be in the same order as fk)

    data =  np.genfromtxt(modefile+'.txt', dtype=typ3)

    QS   =  np.zeros((3*NS,HFPI.MMAX))
    Si0  =  data['Si'][0]            # to mark beginning of new modes
    k    = -1

    for s, qx, qy, qr in data:

        if (s  ==  Si0):
            k  =  k + 1
            
        i  = list(Si).index(s)     # locate store index by store ID

        QS[3*i+0, k] = qx
        QS[3*i+1, k] = qy
        QS[3*i+2, k] = qr

# 5. Normalize modal shapes for unitary modal masses
    
    for k in range(HFPI.MMAX):
        Mk = np.sum(MS*QS[:,k]*QS[:,k])
        QS[:,k] = QS[:,k] / np.sqrt(Mk)

# 6. Check orthogonality of modal shapes
    
    ortho = np.zeros((HFPI.MMAX,HFPI.MMAX))    

    for ii in range(HFPI.MMAX):
        for jj in range(HFPI.MMAX):
            ortho[ii,jj] = np.dot(np.dot(QS[:,ii],np.diag(MS)),QS[:,jj])

# 7. Save pickle file with all structural data
           
    with gz.GzipFile(str_file+'.gz', 'wb') as target:
        pk.dump((NS, Si, XS, MS, fk, QS, ortho), target)

    return NS, Si, XS, MS, fk, QS, ortho


#=============================================================================

def read_TQS(coorfile, massfile, freqfile, modefile, str_file):

    try:
        with gz.GzipFile(str_file+'.gz', 'rb') as target:
            return pk.load(target)
    except:
        pass
    
# 1. Define data type for each column of files to be read           

# 1.1 Nodal coordinates
    typ0  = [('Ni','i'),        # node ID
             ('XN','f'),        #
             ('YN','f'),        # node coordinates in designer reference system
             ('ZN','f')]        #

# 1.2 Lumped nodal masses
    typ1  = [('Ni','i'),        # node ID
             ('MN','f')]        # node mass 
  
# 1.3 Frequencies
    typ2  = [('K', 'i'),        # mode ID (not used)
             ('fk','f')]        # modal frequency (must be sorted!!!)

# 1.4 Modal shapes
    typ3  = [('Ni','i'),        # node ID
             ('qX','f'),        # modal displacement in X
             ('qY','f'),        # modal displacement in Y
             ('qZ','f')]        # modal displacement in Z

# 2. Read nodal coordinates

    dat0 =  np.genfromtxt(coorfile+'.txt', dtype=typ0)

    Ni   =  dat0['Ni']
    NN   =  len(Ni)
    XN   =  np.vstack((dat0['XN'], dat0['YN'], dat0['ZN'])).T

# 3. Build up lumped mass matrices and lumped storey mass

    dat1 =  np.genfromtxt(massfile+'.txt', dtype=typ1)
    
    MNi  =  9810.*dat1['MN']
    MN   =  np.vstack((MNi, MNi, MNi)).T.reshape((3*NN,))

# 4. Read modal frequencies (must be sorted)

    dat2 =  np.genfromtxt(freqfile+'.txt', dtype=typ2)

    fk   =  1/dat2['fk']
    fk   =  fk[0:HFPI.MMAX]
    
# 5. Read modal shapes (must be ascending sorted by fk)

    dat3 =  np.genfromtxt(modefile+'.txt', dtype=typ3)
    print('Start reading modes: ')

    NiN  =  list(Ni)
    NiQ  =  list(dat3['Ni'][0:NN])   # assume all modes with same labels/order
    ni   =  np.empty((NN,), dtype='int')
    QN   =  np.empty((3*NN,HFPI.MMAX))
    
    for ki, n in enumerate(NiN):     # use labels/order from coordinates file
        ni[ki] = NiQ.index(n)
    
    if len(ni) != len(np.unique(ni)):
        sys.exit('Node label not unique in modal shapes file!')
        
    for k in range(HFPI.MMAX):
        print('  Reading mode {0}... '.format(k+1))

        qX = dat3['qX'][k*NN:(k+1)*NN]
        qY = dat3['qY'][k*NN:(k+1)*NN]
        qZ = dat3['qZ'][k*NN:(k+1)*NN]

        QN[:,k] = np.vstack((qX[ni], qY[ni], qZ[ni])).T.reshape((3*NN,))       
        
    print('... done!')
        
# 6. Normalize modal shapes for unitary modal masses

    for k in range(HFPI.MMAX):
        QN[:,k] = QN[:,k] / np.sqrt(np.sum(MN*QN[:,k]*QN[:,k]))

# 7. Check orthogonality of modal shapes
    
    ortho = np.zeros((HFPI.MMAX,HFPI.MMAX))    

    for ii in range(HFPI.MMAX):
        for jj in range(HFPI.MMAX):
            ortho[ii,jj] = np.dot(QN[:,ii],QN[:,jj])

# 8. Save pickle file with all structural data       

    with gz.GzipFile(str_file+'.gz', 'wb') as target:
        pk.dump((NN, Ni, XN, MN, fk, QN, ortho), target)

    return NN, Ni, XN, MN, fk, QN, ortho


#=============================================================================

def offset_TQS(X, Q, offset):

# 1. Prepare translation vector and rotation matrix

    D  = offset[0:3].reshape(3,1)
    rz = offset[3]
    
    R  = np.array([[ np.cos(rz), np.sin(rz),  0.0],
                   [-np.sin(rz), np.cos(rz),  0.0],
                   [  0.0,        0.0,        1.0]], dtype='float')

# 2. Apply transformation to coordinates
                      
    Xoff = (np.dot(R,(X.T  + D))).T

# 3. Apply transformation also to modal shapes (if not a placeholder)

    if len(Q) > 0:

        N3, NQN  =  Q.shape
        Qoff     =  np.empty(Q.shape)
    
        for k in range(NQN):
        
            qk   =  Q[:,k].reshape(N3//3, 3)  
            qoff = (np.dot(R, qk.T)).T
        
            Qoff[:,k] = qoff.ravel()
    
        return Xoff, Qoff
    else:
        return Xoff


#=============================================================================

def TQS_to_3dof(NN, Ni, XN, MN, fk, QN, ZS):

# Contribution by each node

    NN  =  XN.shape[0]
    RNi =  np.sqrt(XN[:,0]**2 + XN[:,1]**2)
    MNi =  MN.reshape(NN,3)[:,0]
    INi =  RNi*RNi*MNi

# Boundaries for collecting nodes from each floor

    ZI = ZS[ 0] + (ZS[ 0] - ZS[ 1])
    ZF = ZS[-1] - (ZS[-2] - ZS[-1])
    
    Zdiff = np.diff(np.hstack((ZI,ZS,ZF)))/2

    if (ZS[0] > ZS[1]):    
        Z2 = ZS - Zdiff[:-1]
        Z1 = ZS + Zdiff[1: ]
    else:    
        Z1 = ZS - Zdiff[:-1]
        Z2 = ZS + Zdiff[1: ]

# Reset collectors

    NS  =  len(ZS)
    Si  =  np.arange(NS)
    MSi =  np.zeros(NS)
    ISi =  np.zeros(NS)
    XS  =  np.zeros((NS,3))
    KS  =  np.zeros((NS,NN), dtype='bool')
    
# Collect

    for s in range(NS):

        kS = (XN[:,2] > Z1[s]) & (XN[:,2] <= Z2[s])
        
        if any(kS):
            KS[s,:] = kS
            MSi[s]  = np.sum(MNi[kS])
            ISi[s]  = np.sum(INi[kS])
    
# Prepare output

    MS      =  np.vstack((MSi, MSi, ISi)).T.reshape((3*NS,))
    XS[:,2] =  ZS
    Nk      =  len(fk)
    QS      =  np.zeros((3*NS, Nk))

    return NS, Si, XS, MS, fk, QS, KS


#=============================================================================

def interp_TQS(XT, XN, QN):

# 1. Get taps coordinates (with offset over XT and XN already applyed)

    NT = XT.shape[0]
    NN = XN.shape[0]
    
    QT = np.empty((3*NT,HFPI.MMAX))

# 2. Loop over all available modes
    
    for k in range(HFPI.MMAX):
        
        qTk = np.zeros((NT, 3))
        qNk = QN[:,k].reshape(NN, 3)

# 4. Perform mode interpolation at each tap

        qTk[:,0] = griddata(XN, qNk[:,0], XT, method='nearest')
        qTk[:,1] = griddata(XN, qNk[:,1], XT, method='nearest')
        qTk[:,2] = griddata(XN, qNk[:,2], XT, method='nearest')
        
        QT[:,k]  = qTk.reshape(3*NT,)

    return QT


#=============================================================================

def view_TQS(H0, Bx, By, Hz, XT, XN, ZS, fk, QT, QN, scale):

# 1. Loop over all modes

    NT  =  XT.shape[0]
    NN  =  XN.shape[0]
    kN  =  XN[:,2] > H0
    
    Z0  =  ZS.min()
    Z1  =  ZS.max()
    DZ  =  Z1 - Z0
    Z0  =  Z0 - 0.1*DZ
    Z1  =  Z1 + 0.1*DZ
    
#   ticks = [-100, -50, 0, 50, 100]
#   ticks = [ -50, -25, 0, 25,  50]
    ticks = [ -30, -15, 0, 15,  30]
    
    for k in range(HFPI.MMAX):
              
        text = r'MODE {0}: $f_k$ = {1:5.2f}Hz'.format(k+1, fk[k])
        
        plt.figure(figsize=(15, 8))
        plt.suptitle(text, fontsize=18)
        
        ax1 = plt.subplot(131)
        ax1.set_xlim([ticks[0], ticks[3]])
        ax1.set_ylim([ticks[0], ticks[3]])
        ax1.set_xticks(ticks)
        ax1.set_yticks(ticks)
        ax1.set_aspect('equal', adjustable='box')
        ax1.grid(True)

        ax2 = plt.subplot(132)
        ax2.grid(True)
        ax2.set_xlim([ticks[0],ticks[3]])
        ax2.set_ylim([Z0, Z1])
        ax2.set_xticks(ticks)
        ax2.set_aspect('equal', adjustable='box')
        
        ax3 = plt.subplot(133)
        ax3.grid(True)
        ax3.set_xlim([ticks[0],ticks[3]])
        ax3.set_ylim([Z0, Z1])
        ax3.set_xticks(ticks)
        ax3.set_aspect('equal', adjustable='box')

# 2. Plot structural nodes

        qk   =  QN[:,k].reshape(NN,3)
        sc   =  scale/np.abs(qk).max()
        qk   =  sc*qk
        Xoff =  XN + qk
        
        ax1.plot(XN[kN,0], XN[kN,1], 'g.', markersize=1)
        ax2.plot(XN[kN,0], XN[kN,2], 'g.', markersize=1)
        ax3.plot(XN[kN,1], XN[kN,2], 'g.', markersize=1)
        
        ax1.plot(Xoff[kN,0], Xoff[kN,1], 'b.', markersize=1)
        ax2.plot(Xoff[kN,0], Xoff[kN,2], 'b.', markersize=1)
        ax3.plot(Xoff[kN,1], Xoff[kN,2], 'b.', markersize=1)

# 3. Plot taps

        if len(XT) == 0: continue

        qk   =  sc*QT[:,k].reshape(NT,3)
        Xoff =  XT + qk

        ax1.plot(Xoff[:,0], Xoff[:,1], 'ro', markersize=4)
        ax2.plot(Xoff[:,0], Xoff[:,2], 'ro', markersize=4)
        ax3.plot(Xoff[:,1], Xoff[:,2], 'ro', markersize=4)

# 4. Plot floors
#       for ks in range(len(ZS)):
#           ax2.plot([-Bx/1.8, Bx/1.8], [ZS[ks], ZS[ks]], 'k')
#           ax3.plot([-By/1.8, By/1.8], [ZS[ks], ZS[ks]], 'k')

        plt.savefig('MODE_{0:0>2}.png'.format(k+1), 
                    format='png', dpi=300, bbox_inches='tight')

    return

#==============================================================================

def mload_TQS(wgraph, master, Cp, XT, QT, Td, plot=False):

    graph,  wnd              =  wgraph        
    dirname, prefix, m,  zts =   graph

# 1. Calculate modal loads for all modes

    CAx =  master[['A', 'cx']].product(axis=1)
    CAy =  master[['A', 'cy']].product(axis=1)
    CAz =  master[['A', 'cz']].product(axis=1)

    n   =  Cp.shape[0]
    fs  =  n/Td
    N3  =  QT.shape[0]
    Pt  =  np.empty((N3,n))
    fPk =  np.zeros(3)
    
    Pt[0::3,:]  = -Cp.multiply(CAx, axis=1).values.T
    Pt[1::3,:]  = -Cp.multiply(CAy, axis=1).values.T
    Pt[2::3,:]  = -Cp.multiply(CAz, axis=1).values.T
        
    Pk  =  np.dot(QT.T,Pt)
    mPk =  np.mean(Pk, axis=1)
    sPk =  np.std (Pk, axis=1)

    if plot:
        
        t =  np.linspace(0,Td,n)
        plt.figure(figsize=(15, 8))        

    for k in range(3):
        
        Pk0    = (Pk[k,1:] - mPk[k])/sPk[k]   # normalize process
        f, SPk =  URP.periodogram(Pk0, fs)    # get spectrum
        SPk    =  URP.mov_average(SPk, 21)    # smooth spectrum

# 3. Plot all modal loads and spectra if required
   
        if plot:
       
            plt.subplot(3,2,2*k+1)
            plt.plot(t,Pk[k,:])
            plt.grid(True)
            plt.ylabel('Mode {0}'.format(k+1))
            
            if k == 0:
                plt.axis([0,Td,-3500,3500])
                plt.text(0.75*Td,2800,
                         r'$\mu_k =$ {0:6.1f}'.format(mPk[k]), fontsize=14)

            if k == 1:
                plt.axis([0,Td,-1500,1500])
                plt.text(0.75*Td,1200,
                         r'$\mu_k =$ {0:6.1f}'.format(mPk[k]), fontsize=14)

            if k == 2:
                plt.axis([0,Td, -1000, 1000])
                plt.xlabel('Time (s)')
                plt.text(0.75*Td, 800,
                         r'$\mu_k =$ {0:6.1f}'.format(mPk[k]), fontsize=14)
                 
            plt.subplot(3,2,2*k+2)
            plt.loglog(f,SPk,'b')#, f,np.exp(SS),'r:')
            plt.grid(True)
            plt.axis([1e-3,1,0.001,100])
            plt.text(0.2,20,
                     r'$\sigma_k =$ {0:6.1f}'.format(sPk[k]), fontsize=14)
            
            if k == 2:
                plt.xlabel('Frequency (Hz)')
            
            if fPk[k] > 0.1:
                plt.loglog([fPk[k], fPk[k]], [0.001, 100], color='r', lw=2)
                plt.text(fPk[k], 10, '{0:.3f}Hz'.format(fPk[k]))
                
# 4. Save plot

    if plot:
        
        tst = prefix + r': Modal loads for '
        tst = tst + r' m = {0:0>2} years, $\alpha$ = {1} deg'
        tst = tst.format(m, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1:0>3}D_Mloads.png'
        fst = fst.format(m, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

    return Pk, mPk, sPk, fPk

    
#=============================================================================

def top_TQS(wgraph, i_Emax, nc, fk, uk, QN, Td, Xd=False, plot=False):
    
    graph,  wnd              =  wgraph        
    dirname, prefix, m,  zts =   graph
    
# 2. Sum up the modal contribution to top displacement
    
    n   =  uk.shape[1]
    NN  =  QN.shape[0]//3
    uTc =  np.zeros((3,n))
    
    fs  =  n/Td
    f0  =  np.max([1./Td, fk[ 0]/10.])
    f1  =  np.min([   fs, fk[-1]*2.0])

    for k in range(HFPI.MMAX):
        qk  =  QN[:,k].reshape(NN,3) 
        uTc[0,:] +=  qk[nc,0]*uk[k,:]
        uTc[1,:] +=  qk[nc,1]*uk[k,:]
    
    uTc[2,:] = np.sqrt(uTc[0,:]**2 + uTc[1,:]**2)
    
# 3. Derivate to acceleration
    
    aTc  =  np.empty(uTc.shape)

    aTc[0,:] = URP.differentiate(uTc[0,:], fs, [f0, f1])
    aTc[0,:] = URP.differentiate(aTc[0,:], fs, [f0, f1])
    aTc[1,:] = URP.differentiate(uTc[1,:], fs, [f0, f1])
    aTc[1,:] = URP.differentiate(aTc[1,:], fs, [f0, f1])
    
    aTc[2,:] = np.sqrt(aTc[0,:]**2 + aTc[1,:]**2)

# 4. Plot top displacements and spectra if required

    if plot:
        t    =  np.linspace(0,Td,n)
        
        '''
# 4.1 Displacement

    muTc =  np.mean(uTc, axis=1)
    suTc =  np.std (uTc, axis=1)

        plt.figure(figsize=(15, 8))
        
        for i in range(3):
       
            plt.subplot(3,2,2*i+1)
            plt.plot(t,uTc[i,:]*1000)
            plt.grid(True)

            if i == 0:
                plt.axis([0,Td,-HFPI.DLIM,HFPI.DLIM])
                plt.ylabel('X displ. (mm)')              
            if i == 1:
                plt.axis([0,Td,-HFPI.DLIM,HFPI.DLIM])
                plt.ylabel('Y displ. (mm)')  
            if i == 2:
                plt.axis([0,Td,0,HFPI.DLIM])
                plt.ylabel('Resultant (mm)')
                plt.xlabel('Time (s)')

            uTc0    = (uTc[i,1:] - muTc[i])/suTc[i]
            f, SuTc =  URP.periodogram(uTc0, fs)
            SuTc    =  URP.mov_average(SuTc, 5)

            plt.subplot(3,2,2*i+2)
            plt.loglog(f,SuTc)
            plt.grid(True)
            plt.axis([1e-2,1,1e-6,1000])
            if i == 2: plt.xlabel('Frequency (Hz)')

        tst = prefix + r': Top displacements for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}, $\alpha$ = {2} deg'
        tst = tst.format(m, zts, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1}_{2:0>3}D_TopDispl.png'
        fst = fst.format(m, zts, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

# 4.2 Acceleration
    
    maTc =  np.mean(aTc, axis=1)
    saTc =  np.std (aTc, axis=1)

        plt.figure(figsize=(15, 8))

        for i in range(3):

            plt.subplot(3,2,2*i+1)
            plt.plot(t,aTc[i,:]*1000/9.81)
            plt.grid(True)

            if i == 0:
                plt.axis([0,Td,-HFPI.ALIM,HFPI.ALIM])
                plt.ylabel('X accel. (mG)')
                
            if i == 1:
                plt.axis([0,Td,-HFPI.ALIM,HFPI.ALIM])
                plt.ylabel('Y accel. (mG)')
               
            if i == 2:                
                plt.axis([0,Td,0,HFPI.ALIM])
                plt.ylabel('Resultant (mG)')
                plt.xlabel('Time (s)')

            aTc0    = (aTc[i,1:] - maTc[i])/saTc[i]
            f, SaTc =  URP.periodogram(aTc0, fs)
            SaTc    =  URP.mov_average(SaTc, 5)
        
            plt.subplot(3,2,2*i+2)
            plt.loglog(f,SaTc)
            plt.grid(True)
            plt.axis([1e-2,1,1e-6,1000])
            if i == 2: plt.xlabel('Frequency (Hz)')

        tst = prefix + r': Top accelerations for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}, $\alpha$ = {2} deg'
        tst = tst.format(m, zts, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1}_{2:0>3}D_TopAccel.png'
        fst = fst.format(m, zts, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')
        '''
            
# 4.5 Top horizontal walk

        uxm = 1000*np.mean(uTc[0,:])
        uym = 1000*np.mean(uTc[1,:])        
        uxp = 1000*uTc[0,i_Emax]
        uyp = 1000*uTc[1,i_Emax]

        if (Xd):
            uxw = HFPI.DLIM*np.cos((wnd + 180)*np.pi/180)
            uyw = HFPI.DLIM*np.sin((wnd + 180)*np.pi/180)
        else:
            uxw = HFPI.DLIM*np.cos( wnd       *np.pi/180)
            uyw = HFPI.DLIM*np.sin( wnd       *np.pi/180)
        
        plt.figure(figsize=(15, 8))
        
        plt.plot(uTc[0,:]*1000,uTc[1,:]*1000,lw=0.2)
        
        plt.plot([0,  uxw],[0,   uyw],'k',lw=1)
        plt.plot([0,  uxm],[0,   uym],'g',lw=1)
        plt.plot([uxm,uxp],[uym, uyp],'r',lw=1)
        
        plt.plot(0.0, 0.0, 'k', marker='o')
        plt.plot(uxm, uym, 'g', marker='o')
        plt.plot(uxp, uyp, 'r', marker='o')
        
        plt.axis([-HFPI.DLIM,HFPI.DLIM,-HFPI.DLIM,HFPI.DLIM])
        plt.axes().set_aspect('equal')
        plt.xlabel('X displ. (mm)')
        plt.ylabel('Y displ. (mm)')
        plt.grid(True)
        
        trst = 't_ref = {0:4}s'.format(np.int(t[i_Emax]))
        plt.text(-3*HFPI.DLIM/4, 3*HFPI.DLIM/4, trst, fontsize=14)
        
        tst = prefix + r': Top walk for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}, $\alpha$ = {2} deg'
        tst = tst.format(m, zts, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1}_{2:0>3}D_TopWalk.png'
        fst = fst.format(m, zts, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

    return uTc, aTc

    
#=============================================================================

def stateq_TQS(fk, ukmax, XN, MN, QN, KS, XTC):

# 1. Sum up the modal contribution to top displacement
    
    NS  =  KS.shape[0]
    NN  =  XN.shape[0]
    
    Peq =  np.zeros((NS,3))
    
    for s  in range(NS):
        Xs  = XN[KS[s,:],:] - XTC
        Ms  = MN.reshape(NN,3)[KS[s,:],:]
        Fs  = np.zeros(Ms.shape)

        for k in range(HFPI.MMAX):
            qk  = QN[:,k].reshape(NN,3)[KS[s,:],:]
            Fs += ukmax[k]*((2*np.pi*fk[k])**2)*Ms*qk

        Fxs = np.sum(Fs[:,0])
        Fys = np.sum(Fs[:,1])
        Mzs = np.sum(Fs[:,1]*Xs[:,0] - Fs[:,0]*Xs[:,1])

        Peq[s,:] = [Fxs, Fys, Mzs]
            
    return Peq.ravel()

    
#=============================================================================

def compare_loads(graph, Hz, ZS, alpha, FXdf, FYdf, MZdf):
    
    dirname, prefix, m,  zts =   graph
    static_loads = dirname + prefix + '_Static.xlsx'
    
    Z0  =  ZS.min()
    Z1  =  ZS.max()
    DZ  =  Z1 - Z0
    Z0  =  Z0 - 0.1*DZ
    Z1  =  Z1 + 0.1*DZ
    
# 1. Force in X direction

    FXst = pd.read_excel(static_loads, sheetname='FX')
    
    FXmin = min((FXst.values.min(),FXdf.values.min()))
    FXmax = max((FXst.values.max(),FXdf.values.max()))    
    
    plt.figure(figsize=(15, 8))
    for k in range(12):
        
        plt.subplot(3,4,k+1)
        plt.plot(FXst[15*k].values, ZS,'r')
        plt.plot(FXdf[15*k].values, ZS,'b')
        plt.axis([FXmin, FXmax, Z0, Z1])
        plt.text(FXmax/2, Z1/2, r'$\alpha = $' + str(alpha[k]))
        plt.grid(True)
        
    pst = prefix + '_FX_Comparison_A.png'
    plt.suptitle(pst, fontsize=18)
    plt.savefig(pst,format='png', dpi=300, bbox_inches='tight')    
    
    plt.figure(figsize=(15, 8))
    for k in range(12):
        
        plt.subplot(3,4,k+1)
        plt.plot(FXst[15*(k+12)].values, ZS,'r')
        plt.plot(FXdf[15*(k+12)].values, ZS,'b')
        plt.axis([FXmin, FXmax, Z0, Z1])
        plt.text(FXmax/2, Z1/2, r'$\alpha = $' + str(alpha[k+12]))
        plt.grid(True)
        
    pst = prefix + '_FX_Comparison_B.png'
    plt.suptitle(pst, fontsize=18)
    plt.savefig(pst,format='png', dpi=300, bbox_inches='tight')    
        
# 2. Force in Y direction

    FYst = pd.read_excel(static_loads, sheetname='FY')
 
    FYmin = min((FYst.values.min(),FYdf.values.min()))
    FYmax = max((FYst.values.max(),FYdf.values.max()))    
    
    plt.figure(figsize=(15, 8))
    for k in range(12):
        
        plt.subplot(3,4,k+1)
        plt.plot(FYst[15*k].values, ZS,'r')
        plt.plot(FYdf[15*k].values, ZS,'b')
        plt.axis([FYmin, FYmax, Z0, Z1])
        plt.text(FYmax/2, Z1/2, r'$\alpha = $' + str(alpha[k]))
        plt.grid(True)
        
    pst = prefix + '_FY_Comparison_A.png'
    plt.suptitle(pst, fontsize=18)
    plt.savefig(pst,format='png', dpi=300, bbox_inches='tight')    
    
    plt.figure(figsize=(15, 8))
    for k in range(12):
        
        plt.subplot(3,4,k+1)
        plt.plot(FYst[15*(k+12)].values, ZS,'r')
        plt.plot(FYdf[15*(k+12)].values, ZS,'b')
        plt.axis([FYmin, FYmax, Z0, Z1])
        plt.text(FYmax/2, Z1/2, r'$\alpha = $' + str(alpha[k+12]))
        plt.grid(True)
        
    pst = prefix + '_FY_Comparison_B.png'
    plt.suptitle(pst, fontsize=18)
    plt.savefig(pst,format='png', dpi=300, bbox_inches='tight')    
    
# 4. Moment in Z direction

    MZst = pd.read_excel(static_loads, sheetname='MZ')

    MZmin = min((MZst.values.min(),MZdf.values.min()))
    MZmax = max((MZst.values.max(),MZdf.values.max()))    
    
    plt.figure(figsize=(15, 8))
    for k in range(12):
        
        plt.subplot(3,4,k+1)
        plt.plot(MZst[15*k].values, ZS,'r')
        plt.plot(MZdf[15*k].values, ZS,'b')
        plt.axis([MZmin, MZmax, Z0, Z1])
        plt.text(MZmax/2, Z1/2, r'$\alpha = $' + str(alpha[k]))
        plt.grid(True)

    pst = prefix + '_MZ_Comparison_A.png'
    plt.suptitle(pst, fontsize=18)
    plt.savefig(pst,format='png', dpi=300, bbox_inches='tight')    
    
    plt.figure(figsize=(15, 8))
    for k in range(12):
        
        plt.subplot(3,4,k+1)
        plt.plot(MZst[15*(k+12)].values, ZS,'r')
        plt.plot(MZdf[15*(k+12)].values, ZS,'b')
        plt.axis([MZmin, MZmax, Z0, Z1])
        plt.text(MZmax/2, Z1/2, r'$\alpha = $' + str(alpha[k+12]))
        plt.grid(True)

    pst = prefix + '_MZ_Comparison_B.png'
    plt.suptitle(pst, fontsize=18)
    plt.savefig(pst,format='png', dpi=300, bbox_inches='tight')    

    return FXst, FYst, MZst
    
#=============================================================================

# -*- coding: utf-8 -*-

import sys
import gzip   as gz
import pickle as pk
import pandas as pd
import numpy  as np

import matplotlib.pyplot as     plt
from   scipy.interpolate import griddata

import URP
import NBR6123

#=============================================================================

ALIM = 40.     # limit for acceleration plots (mG)
DLIM = 500.    # limit for displacement plots (mm)
RLIM = 1.      # limit for rotation plots (degrees)
ELIM = 1000    # limit for energy plots (kJ)
MMAX = 12      # maximum number of modes to be considered

#=============================================================================

def read_batch(dirname, prefix):

    batchpickle = dirname + prefix + '.pk'

    try:
        with open(batchpickle, 'rb') as target:
            print('Loading "{0}" ... '.format(dirname))
            return pk.load(target)
            
    except:
        sys.exit('Could not read batch file {0}!'.format(batchpickle))


#=============================================================================

def get_taps(b):
    
    print('  Direction = {0:0>3}deg: "{1}"... '.format(b.wnd, b.file))
    
    with gz.GzipFile(b.file, 'rb') as target:
        return pk.load(target)


#=============================================================================

def unpack_batch(batch):

    wnd = batch['wnd']
    Cat = batch['Cat']
    Vm  = batch['Vm']
    dPa = batch['dPa']
    Td  = batch['Td']
        
    if batch.ndim == 2:        
        wnd = wnd.values.astype('int')
        Cat = Cat.values.astype('int')
        Vm  =  Vm.values.astype('float')
        dPa = dPa.values.astype('float')
        Td  =  Td.values.astype('float')

    return wnd, Cat, Vm, dPa, Td


#=============================================================================

def unpack_master(master):
    
    tap  = master['tap'].astype('int')
    A    = master['A'].astype('float')
    XT   = master[['X', 'Y', 'Z' ]].values.astype('float')
    CT   = master[['cx','cy','cz']].values.astype('float')

    if master.ndim == 2:
        tap = tap.values.astype('int')
        A   = A.values.astype('float')

    return tap, A, XT, CT


#=============================================================================

def build_3dof(Bx, By, Hz, rho, gamma, NS):

# 1. Define storey coordinates and masses          

    dz   =  Hz/NS
    Si   =  np.linspace(1,  NS, NS)
    ZS   =  np.linspace(dz, Hz, NS)
    X0   =  np.zeros(ZS.shape)
    Mz   = (Bx*By*dz*rho)*np.ones(NS)
    Iz   = (Bx**2 + By**2)*Mz/12

    XS   =  np.vstack(( X0, X0, ZS)).T        # stifness center for each floor
    MS   =  np.vstack(( Mz, Mz, Iz)).T.reshape((3*NS,))   # masses and inertia

# 2. Prepare the 3 modal shapes

    NQ   =  3
    QS   =  np.zeros((3*NS,NQ))

    for k in range(NQ):

        QS[k:3*NS:3, k] = (ZS/Hz)**gamma[k]

# 5. Normalize modal shapes for unitary modal masses

        Mk = np.sum(MS*QS[:,k]*QS[:,k])
        QS[:,k] = QS[:,k] / np.sqrt(Mk)

    return Si, XS, MS, QS


#=============================================================================

def interp_3dof(XT, XS, QS):

# 1. Get taps coordinates (with offset over XT and ZS already applyed)
#    Note: verify if extrapolation will be necessary beyond ZS.

    NT = XT.shape[0]
    NS = XS.shape[0]    
    
    ZT = np.array(XT[:,2], dtype=np.float)
    ZS = np.array(XS[:,2], dtype=np.float)

    top = False
    if (ZS.max() < ZT.max()):
            ZS  = np.append(ZS, ZT.max()+0.1)
            top = True
            
    bot = False
    if (ZS.min() > ZT.min()):
            ZS  = np.append(ZT.min()-0.1, ZS)
            bot = True

    QT = np.empty((3*NT,MMAX))

# 2. Loop over valid modes
    
    for k in range(MMAX):
        
        qTk = np.zeros((NT, 3))
        qSk = QS[:,k].reshape(NS, 3)

# 3. If necessary, extrapolate modal shape as constant
        
        if top:
            qSk = np.vstack((qSk, qSk[-1,:]))
            
        if bot:
            qSk = np.vstack((qSk[ 0,:], qSk))

# 4. Perform mode interpolation at each tap
        
        r        = griddata(ZS, qSk[:,2], ZT, method='linear')
        qTk[:,0] = griddata(ZS, qSk[:,0], ZT, method='linear') - r*XT[:,1]
        qTk[:,1] = griddata(ZS, qSk[:,1], ZT, method='linear') + r*XT[:,0]

        QT[:,k]  = qTk.reshape(3*NT,)

    return QT


#=============================================================================

def view_3dof(Bx, By, Hz, fk, XT, QT, XS, QS, scale):

# 1. Loop over all modes

    NT  =  XT.shape[0]
    NS  =  XS.shape[0]
    
    BX  =  np.vstack((np.array([-Bx/2., -Bx/2.,  Bx/2.,  Bx/2., -Bx/2.]),
                      np.array([-By/2.,  By/2.,  By/2., -By/2., -By/2.]),
                      np.array([ 0.00 ,  0.00,   0.00,   0.00,   0.00 ]))).T

    gray = [0.8,0.8,0.8]
    blue = [0.0,0.0,1.0]
    
    for k in range(MMAX):
              
        text = r'MODE {0}: $f_k$ = {1:5.2f}Hz'.format(k+1, fk[k])
        
        plt.figure(figsize=(15, 8))
        plt.suptitle(text, fontsize=18)
        
        ax1 = plt.subplot(131)
        ax1.set_xlim([-1.0*Bx, 1.0*Bx])
        ax1.set_ylim([-1.0*Bx, 1.0*Bx])
        ax1.set_aspect('equal', adjustable='box')
        ax1.grid(True)

        ax2 = plt.subplot(132)
        ax2.grid(True)
        ax2.set_xlim([-1.0*Bx, 1.0*Bx])
        ax2.set_ylim([-0.05*Hz, 1.1*Hz])
        ax2.set_xticks([-40, -20, 0, 20, 40])
        ax2.set_aspect('equal', adjustable='box')
        
        ax3 = plt.subplot(133)
        ax3.grid(True)
        ax3.set_xlim([-1.0*Bx, 1.0*Bx])
        ax3.set_ylim([-0.05*Hz, 1.1*Hz])
        ax3.set_xticks([-40, -20, 0, 20, 40])
        ax3.set_aspect('equal', adjustable='box')

# 2. Plot frames

        qkT =  QT[:,k].reshape(NT,3)
        sc  =  scale/np.abs([qkT[:,0], qkT[:,1]]).max()
        qk  =  sc*QS[:,k].reshape(NS,3)
        
        for xs, qs  in zip(XS, qk):
            
            off =  np.array([qs[0], qs[1], 0.0, -qs[2]])
            BS  =  offset_taps(BX, [], off)
            ZS  =  xs[2]*np.ones((5,))
        
            ax1.plot(BX[:,0], BX[:,1], color=gray, linewidth=2)
            ax1.plot(BS[:,0], BS[:,1], color=blue, linewidth=2)
            
            ax2.plot(BX[:,0], ZS,      color=gray, linewidth=2)
            ax2.plot(BS[:,0], ZS,      color=blue, linewidth=2)
            
            ax3.plot(BX[:,1], ZS,      color=gray, linewidth=2)
            ax3.plot(BS[:,1], ZS,      color=blue, linewidth=2)

# 3. Plot taps

        if len(XT) == 0: continue

        qk   =  sc*QT[:,k].reshape(NT,3)
        Xoff =  XT + qk

        ax1.plot(Xoff[:,0], Xoff[:,1], 'ro')
        ax2.plot(Xoff[:,0], Xoff[:,2], 'ro')
        ax3.plot(Xoff[:,1], Xoff[:,2], 'ro')

        plt.savefig('MODE_{0:0>2}.png'.format(k+1), 
                    format='png', dpi=300, bbox_inches='tight')

    return
    

#=============================================================================

def offset_taps(XT, CT, offset):

# 1. Prepare translation vector and rotation matrix

    D  = offset[0:3].reshape(3,1)
    rz = offset[3]
    
    R  = np.array([[ np.cos(rz), np.sin(rz),  0.0],
                   [-np.sin(rz), np.cos(rz),  0.0],
                   [  0.0,        0.0,        1.0]], dtype='float')

# 2. Apply transformation to coordinates

    if len(CT):
        return (np.dot(R,(XT.T + D))).T, (np.dot(R, CT.T)).T
    else:
        return (np.dot(R,(XT.T + D))).T


#=============================================================================

def aero_coeffs(graph, batch, structure, C_eqv=[], plot=False):
        
    dirname, prefix, m, zts = graph
    
# 1. Preliminaries, and get NBR6123 drag coefficients

    stpref =  structure[ 'prefix']
    Bx     =  structure[ 'Bx']
    By     =  structure[ 'By']
    Hz     =  structure[ 'Hz']
    H0     =  structure[ 'H0']
    H1     =  structure[ 'H1']
    
    offT    = structure[['dxT', 'dyT', 'dzT', 'rzT']].values.astype('float')
    offT[3] = np.pi*offT[3]/180
    
    alpha  =  batch['wnd'].values.astype('int64') 
    
    Cm     =  np.empty((len(batch),3))
    Cs     =  np.empty((len(batch),3))
    
    Cx, Cy =  NBR6123.drag(Bx, By, Hz, case='low')
    
    if Bx > By: Ct = 0.075*Cy*Bx/By
    else:       Ct = 0.075*Cx*By/Bx

# 2. Iterate over all wind directions

    for ib, b  in batch.iterrows():

        wnd, Cat, Vm, dPa, Td =  unpack_batch(b)
     
# 3. Open pressure data for wind direction

        master, Cp     =  get_taps(b)
        tap, A, XT, CT =  unpack_master(master)
        
        XT, CT  =  offset_taps(XT, CT, offT)      ##########
        
        A[(XT[:,2] < H0) | (XT[:,2] > H1)] = 0

# 4. Time series for total forces and moments at ground level

        CAx =  A*CT[:,0]
        CAy =  A*CT[:,1]
        CAt =  XT[:,0]*CAy - XT[:,1]*CAx

        Cx0 = -Cp.dot(CAx)/(By*(H1-H0))
        Cy0 = -Cp.dot(CAy)/(Bx*(H1-H0))
        Ct0 = -Cp.dot(CAt)/(By*Bx*(H1-H0))

# 5. Reset mean and rms collection of coefficients

        Cm[ib,0] =  Cx0.mean()
        Cm[ib,1] =  Cy0.mean()
        Cm[ib,2] =  Ct0.mean()
        
        Cs[ib,0] =  Cx0.std()
        Cs[ib,1] =  Cy0.std()
        Cs[ib,2] =  Ct0.std()

# 6. Plot results if required

    if plot:
        xticks =  30*np.arange(13)
        grid   = {'linestyle':'-','linewidth':'0.2'}
        leg2   = {'loc':4, 'fancybox':True, 'ncol':2, 'fontsize':10}
        leg3   = {'loc':4, 'fancybox':True, 'ncol':3, 'fontsize':10}
                   
# 6.1 Plot mean coefficients
            
        plt.figure(figsize=(14, 8))

        yticks = [-2.0, -1.5, -1.0, -0.5,  0.0, 0.5, 1.0, 1.5, 2.0]
        tticks = [-0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5]
    
        ax1 = plt.subplot(221)
        ax1.plot(alpha,      Cm[:,0],  color='b', marker='o')
        
        if len(C_eqv) > 0:
            ax1.plot(alpha, C_eqv[:,0], color='k', ls=':', lw=2)

        ax1.plot([0, 360], [ Cx,  Cx], color='r', ls='--', lw=2)
        ax1.plot([0, 360], [-Cx, -Cx], color='r', ls='--', lw=2)
        
        ax1.grid(**grid)
        ax1.set_xlabel('Wind direction (deg)')
        ax1.set_xlim([0, 360])
        ax1.set_xticks(xticks)
        ax1.set_ylabel(r'mean $C_x$', fontsize=14)
        ax1.set_yticks(yticks)
    
        ax2 = plt.subplot(222)
        ax2.plot(alpha, Cm[:,1], color='b', marker='o')

        if len(C_eqv) > 0:
            ax2.plot(alpha, C_eqv[:,1], color='k', ls=':', lw=2)

        ax2.plot([0, 360], [ Cy,  Cy], color='r', ls='--', lw=2)
        ax2.plot([0, 360], [-Cy, -Cy], color='r', ls='--', lw=2)

        ax2.grid(**grid)
        ax2.set_xlabel('Wind direction (deg)')
        ax2.set_xlim([0, 360])
        ax2.set_xticks(xticks)
        ax2.set_ylabel(r'mean $C_y$', fontsize=14)
        ax2.set_yticks(yticks)
    
        ax3 = plt.subplot(223)
        ax3.plot(alpha, Cm[:,2], color='b', marker='o')
        
        if len(C_eqv) > 0:
            ax3.plot(alpha, C_eqv[:,2], color='k', ls=':', lw=2)

        ax3.plot([0, 360], [ Ct,  Ct], color='r', ls='--', lw=2)
        ax3.plot([0, 360], [-Ct, -Ct], color='r', ls='--', lw=2)
        
        ax3.grid(**grid)
        ax3.set_xlabel('Wind direction (deg)')
        ax3.set_xlim([0, 360])
        ax3.set_xticks(xticks)
        ax3.set_ylabel(r'mean $C_t$', fontsize=14)
        ax3.set_yticks(tticks)

        if len(C_eqv) > 0:
            
            ax1.legend(('Mean', 'HFPI', 'NBR'), **leg3)
            ax2.legend(('Mean', 'HFPI', 'NBR'), **leg3)
            ax3.legend(('Mean', 'HFPI', 'NBR'), **leg3)

        else:
            ax1.legend(('Mean', 'NBR'), **leg2)
            ax2.legend(('Mean', 'NBR'), **leg2)
            ax3.legend(('Mean', 'NBR'), **leg2)

        img = plt.imread('Sections/' + stpref + '_plan.png')
        
        ax4 = plt.subplot(224)
        ax4.imshow(img)
        plt.gca().axison = False
        
        tst = prefix + ': Mean coefficients for'
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}'
        tst = tst.format(m, zts)
        
        fst = prefix + '_{0:0>3}Y_{1}_AeroCoeffs.png'
        fst = fst.format(m, zts)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

        '''
# 6.2 Plot rms coefficients

        plt.figure(figsize=(14, 8))

        yticks = [ 0.00, 0.05, 0.10, 0.15]
        tticks = [ 0.00, 0.01, 0.02, 0.03]
    
        ax1 = plt.subplot(221)
        ax1.plot(alpha, Cs[:,0], color='b', marker='o')
        ax1.grid(**grid)
        ax1.set_xlabel('Wind direction (deg)')
        ax1.set_xlim([0, 360])
        ax1.set_xticks(xticks)
        ax1.set_ylabel(r'rms $C_x$', fontsize=14)
        ax1.set_yticks(yticks)
    
        ax2 = plt.subplot(222)
        ax2.plot(alpha, Cs[:,1], color='b', marker='o')
        ax2.grid(**grid)
        ax2.set_xlabel('Wind direction (deg)')
        ax2.set_xlim([0, 360])
        ax2.set_xticks(xticks)
        ax2.set_ylabel(r'rms $C_y$', fontsize=14)
        ax2.set_yticks(yticks)
    
        ax3 = plt.subplot(223)
        ax3.plot(alpha, Cs[:,2], color='b', marker='o')
        ax3.grid(**grid)
        ax3.set_xlabel('Wind direction (deg)')
        ax3.set_xlim([0, 360])
        ax3.set_xticks(xticks)
        ax3.set_ylabel(r'rms $C_t$', fontsize=14)
        ax3.set_yticks(tticks)

        img = plt.imread('Sections/' + stpref + '_plan.png')
        
        ax4 = plt.subplot(224)
        ax4.imshow(img)
        plt.gca().axison = False
        
        tst = prefix + r': RMS coefficients for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}'
        tst = tst.format(m, zts)
        
        fst = prefix + '_{0:0>3}Y_{1}_RMSCoeffs.png'
        fst = fst.format(m, zts)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')
        '''

    return alpha, Cm, Cs


#==============================================================================

def mload_3dof(wgraph, master, Cp, XT, QT, Td, plot=False):

    graph,  wnd              =  wgraph        
    dirname, prefix, m,  zts =   graph

# 1. Calculate modal loads for all modes

    CAx =  master[['A', 'cx']].product(axis=1)
    CAy =  master[['A', 'cy']].product(axis=1)
    CAt =  XT[:,0]*CAy - XT[:,1]*CAx

    n   =  Cp.shape[0]
    fs  =  n/Td
    N3  =  QT.shape[0]
    Pt  =  np.empty((N3,n))
    fPk =  np.zeros(3)
    
    Pt[0::3,:]  = -Cp.multiply(CAx, axis=1).values.T
    Pt[1::3,:]  = -Cp.multiply(CAy, axis=1).values.T
    Pt[2::3,:]  = -Cp.multiply(CAt, axis=1).values.T
        
    Pk  =  np.dot(QT.T,Pt)
    mPk =  np.mean(Pk, axis=1)
    sPk =  np.std (Pk, axis=1)

    if plot:
        
        t =  np.linspace(0,Td,n)
        plt.figure(figsize=(15, 8))        

# 2. Search for vortex shedding frequency

    for k in range(3):
        
        Pk0    = (Pk[k,1:] - mPk[k])/sPk[k]   # normalize process
        f, SPk =  URP.periodogram(Pk0, fs)    # get spectrum
        SPk    =  URP.mov_average(SPk, 21)    # smooth spectrum
  
        '''
        i0 =  np.argmin(np.abs(f -  0.15))
        f0 =  np.log(f[i0])
        S0 =  np.log(SPk[i0])

        i1 =  np.argmin(np.abs(f -  0.80))
        f1 =  np.log(f[i1])
        S1 =  np.log(SPk[i1])

        SS     =  S0 + (S1 - S0)*(np.log(f) - f0)/(f1 - f0)
        bump   =  np.log(SPk) - SS
        ik     =  np.argmax(bump)
        
        fPk[k] =  f[ik]
        if (fPk[k] < 0.15) | (bump[ik] < 1.6): fPk[k] = 0.
        '''

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

def modal_response(wgraph, fk, zt, Td, Pk, plot=False):

    graph,  wnd              =  wgraph        
    dirname, prefix, m,  zts =   graph

# 1. Calculate modal response for all modes

    n   =  Pk.shape[1]
    fs  =  n/Td
    uk  =  np.zeros((MMAX,n))
    fuk =  np.zeros(MMAX)

    for k in range(MMAX):
        t, uk[k,:] = URP.EQ_Fourier(2*np.pi*fk[k], zt, fs, Pk[k,:])

    muk =  np.mean(uk, axis=1)
    suk =  np.std (uk, axis=1)

    for k in range(MMAX):
        
        uk0    = (uk[k,1:] - muk[k])/suk[k]
        f, Suk =  URP.periodogram(uk0, fs)
        Suk    =  URP.mov_average(Suk, 11)
        fuk[k] =  f[np.argmax(Suk)]
   
# 2. Plot all modal responses and spectra if required

    if plot:
        plt.figure(figsize=(15, 8))
        
        for k in range(3):
       
            plt.subplot(3,2,2*k+1)
            plt.plot(t,uk[k,:])
            plt.grid(True)
            plt.ylabel('Mode {0}'.format(k+1))
            
            if k == 0:
                plt.axis([0,Td,-2500,2500])
                plt.text(0.8*Td,2000,
                         r'$\mu_k =$ {0:6.1f}'.format(muk[k]), fontsize=14)
                
            if k == 1:
                plt.axis([0,Td,-2500,2500])
                plt.text(0.8*Td,2000,
                         r'$\mu_k =$ {0:6.1f}'.format(muk[k]), fontsize=14)

            if k == 2:
                plt.axis([0,Td,-1000,1000])
                plt.text(0.8*Td, 800,
                         r'$\mu_k =$ {0:6.1f}'.format(muk[k]), fontsize=14)
                plt.xlabel('Time (s)')

            uk0    = (uk[k,1:] - muk[k])/suk[k]
            f, Suk =  URP.periodogram(uk0, fs)
            Suk    =  URP.mov_average(Suk, 11)
        
            plt.subplot(3,2,2*k+2)
            plt.loglog(f,Suk)
            plt.grid(True)
            plt.axis([1e-3,1,1e-6,1000])
            plt.text(0.1,100,
                     r'$\sigma_k =$ {0:6.1f}'.format(suk[k]), fontsize=14)

            if k == 2:
                plt.xlabel('Frequency (Hz)')
        
        tst = prefix + r': Modal displacements for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}, $\alpha$ = {2} deg'
        tst = tst.format(m, zts, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1}_{2:0>3}D_Mdispl.png'
        fst = fst.format(m, zts, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')
   
    return uk, muk, suk, fuk


#=============================================================================

def maximum_energy(wgraph, fk, uk, fs, Tav, QN, ne, plot=True):

    graph,  wnd              =  wgraph        
    dirname, prefix, m,  zts =   graph
    
# 1. Sum up the potential elastic energy for each mode and get
#    total displacement at a reference node
    
    n   =  uk.shape[1]
    NN  =  QN.shape[0]//3
    Td  =  n*fs
    
    Ep  =  np.zeros(n)
    uX  =  np.zeros(n)
    uXm =  0.0
    uY  =  np.zeros(n)
    uYm =  0.0
    
    for k in range(MMAX):
        
        qk   =  QN[:,k].reshape(NN,3)
        
        uXm +=  qk[ne,0]*(uk[k,:].mean())
        uX  +=  qk[ne,0]* uk[k,:]

        uYm +=  qk[ne,1]*(uk[k,:].mean())
        uY  +=  qk[ne,1]* uk[k,:]

        Ep  += (((2*np.pi*fk[k])**2) * uk[k,:]**2)/2

    mdl =  np.sqrt(uXm*uXm + uYm*uYm)

    uXm =  uXm/mdl
    uYm =  uYm/mdl      
    
    Ep  =  Ep/1000

# 2. Define statistical peak for observation time Tav < Td

    nmax, Xmax, Eth =  URP.splitmax(Ep, fs, Tav)
    isort           =  np.argsort(np.abs(Ep - Eth))

    for ii in range(n):
        
        i_Emax =  isort[ii]

        uXi    =  uX[i_Emax]
        uYi    =  uY[i_Emax]

        mdl =  np.sqrt(uXi*uXi + uYi*uYi)
        
        uXi =  uXi/mdl
        uYi =  uYi/mdl      
        
        prj = np.abs(np.sum(uXi*uXm + uYi*uYm))
        
        if prj > 0.95: break
    
    Emax   =  Ep[i_Emax]
    sE     =  np.std(Ep)
    mE     =  np.mean(Ep)
    
#   print(prj, Emax/Eth)
    
# 3. Plot potential elastic energy if required

    if plot:
        leg  = {'loc':8, 'fancybox':True, 'ncol':4, 'fontsize':12}
        
        t    =  np.linspace(0, Td, n)
        tmax =  t[i_Emax]
        
        plt.figure(figsize=(15, 8))
        plt.plot(t,Ep)
        plt.plot([0, Td],[mE,  mE ], 'k', lw=2)       
        plt.plot([0, Td],[sE,  sE ], 'g', lw=2, ls='--')       
        plt.plot([0, Td],[Eth, Eth], 'r', lw=2, ls='--')       
        plt.plot( tmax,   Emax,      'r', marker='o')
        
        plt.legend(('instantaneous','mean','rms','600s peak'),**leg)
        
        plt.text(1.01*Td,Eth, '{0:4}kJ'.format(np.int(Eth)), fontsize=14)
        plt.text(1.01*Td, mE, '{0:4}kJ'.format(np.int(mE )), fontsize=14)
        plt.text(1.01*Td, sE, '{0:4}kJ'.format(np.int(sE )), fontsize=14)
        
        trst = 't_ref = {0:4}s'.format(np.int(tmax))
        plt.text(0.85*tmax, 1.05*Emax, trst, fontsize=14)
        
        plt.xlabel('Time (s)')
        plt.ylabel('Energy (kJ)')
        plt.grid(True)
        plt.axis([0,Td,0,ELIM])
        
        tst = prefix + r': Elastic Energy for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}, $\alpha$ = {2} deg'
        tst = tst.format(m, zts, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1}_{2:0>3}D_Energy.png'
        fst = fst.format(m, zts, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

    return i_Emax, Emax, uk[:,i_Emax]


#=============================================================================

def top_3dof(wgraph, i_Emax, floor, edge, fk, uk, QS, Td, plot=False):

    graph,  wnd              =  wgraph        
    dirname, prefix, m,  zts =   graph

# 2. Sum up the modal contribution to top displacement
    
    n   =  uk.shape[1]
    uT  =  np.zeros((3,n))
    
    fs  =  n/Td
    f0  =  np.max([1./Td, fk[ 0]/10.])
    f1  =  np.min([   fs, fk[-1]*2.0])
    
    for k in range(MMAX):
        for i in range(3):
            uT[i,:] +=  QS[3*(floor-1)+i,k]*uk[k,:]
   
    muT =  np.mean(uT, axis=1)
    suT =  np.std (uT, axis=1)
    
# 3. Integrate to acceleration
    
    aT  =  np.empty(uT.shape)

    for i in range(3):
        aT[i,:] = URP.differentiate(uT[i,:], fs, [f0, f1])
        aT[i,:] = URP.differentiate(aT[i,:], fs, [f0, f1])

    aT[2,:] = np.sqrt((aT[0,:] - aT[2,:]*edge[1])**2 + 
                      (aT[1,:] + aT[2,:]*edge[0])**2)
    
    maT =  np.mean(aT, axis=1)
    saT =  np.std (aT, axis=1)

# 4. Plot top displacements and spectra if required

    if plot:
        t    =  np.linspace(0,Td,n)
        
# 4.1 Displacement

        plt.figure(figsize=(15, 8))
        
        for i in range(3):
       
            plt.subplot(3,2,2*i+1)

            if i == 0:
                plt.plot(t,uT[i,:]*1000)
                plt.axis([0,Td,-DLIM,DLIM])
                plt.ylabel('X displ. (mm)')
                
            if i == 1:
                plt.plot(t,uT[i,:]*1000)
                plt.axis([0,Td,-DLIM,DLIM])
                plt.ylabel('Y displ. (mm)')
     
            if i == 2:
                plt.plot(t,uT[i,:]*180/np.pi)
                plt.axis([0,Td,-RLIM,RLIM])
                plt.ylabel('Rotation (deg)')
                plt.xlabel('Time (s)')

            uT0    = (uT[i,1:] - muT[i])/suT[i]
            f, SuT =  URP.periodogram(uT0, fs)
            SuT    =  URP.mov_average(SuT, 5)
        
            plt.subplot(3,2,2*i+2)
            plt.loglog(f,SuT)
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

        plt.figure(figsize=(15, 8))

        for i in range(3):

            plt.subplot(3,2,2*i+1)
            plt.plot(t,aT[i,:]*1000/9.81)

            if i == 0:
                plt.axis([0,Td,-ALIM,ALIM])
                plt.ylabel('X accel. (mG)')
                
            if i == 1:
                plt.axis([0,Td,-ALIM,ALIM])
                plt.ylabel('Y accel. (mG)')
               
            if i == 2:                
                plt.axis([0,Td,0,ALIM])
                plt.ylabel('Edge Accel. (mG)')
                plt.xlabel('Time (s)')

            aT0    = (aT[i,1:] - maT[i])/saT[i]
            f, SaT =  URP.periodogram(aT0, fs)
            SaT    =  URP.mov_average(SaT, 5)
        
            plt.subplot(3,2,2*i+2)
            plt.loglog(f,SaT)
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

# 4.3 Top horizontal walk

        plt.figure(figsize=(15, 8))
        plt.plot(uT[0,:]*1000,uT[1,:]*1000)
        plt.plot(uT[0,i_Emax]*1000,uT[1,i_Emax]*1000, 'r', marker='o')
        plt.axis([-1500,1500,-1500,1500])
        plt.axes().set_aspect('equal')
        plt.xlabel('X displ. (mm)')
        plt.ylabel('Y displ. (mm)')
        plt.grid(True)
        
        trst = 't_ref = {0:4}s'.format(np.int(t[i_Emax]))
        plt.text(-450, 1300, trst, fontsize=14)
        
        tst = prefix + r': Top walk for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}, $\alpha$ = {2} deg'
        tst = tst.format(m, zts, wnd)
        
        fst = prefix + '_{0:0>3}Y_{1}_{2:0>3}D_TopWalk.png'
        fst = fst.format(m, zts, wnd)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

    return uT, aT

   
#=============================================================================

def peak_response(graph, alpha, u_max, a_max):

    dirname, prefix, m,  zts =   graph

    plt.figure(figsize=(15, 8))
    xticks =  30*np.arange(13)
    
# Displacements

    for i in range(3):

        plt.subplot(3,1,i+1)
        plt.plot(alpha, u_max[:,i], color='b', marker='o')
        plt.axis([0,360,0,DLIM])

        if i == 0:
            plt.ylabel('X displ. (mm)')
        if i == 1:
            plt.ylabel('Y displ. (mm)')
        if i == 2:                
            plt.ylabel('Resultant (mm)')
            plt.xlabel('Wind direction (deg)')

        plt.xticks(xticks)                
        plt.grid(True)
    
        tst = prefix + r': Maximum top displacements for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}'
        tst = tst.format(m, zts)
        
        fst = prefix + '_{0:0>3}Y_{1}_PeakDispl.png'
        fst = fst.format(m, zts)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')
                     
# Accelerations

    plt.figure(figsize=(15, 8))

    for i in range(3):

        plt.subplot(3,1,i+1)
        plt.plot(alpha, a_max[:,i], color='b', marker='o')
        plt.axis([0, 360, 0, ALIM])

        if i == 0: 
            plt.ylabel('X accel. (mG)')
        if i == 1: 
            plt.ylabel('Y accel. (mG)')
        if i == 2:                
            plt.ylabel('Resultant (mG)')
            plt.xlabel('Wind direction (deg)')
                
        plt.xticks(xticks)
        plt.grid(True)
    
        tst = prefix + r': Maximum top accelerations for '
        tst = tst + r' m = {0:0>2} years, $\zeta$ = {1}'
        tst = tst.format(m, zts)
        
        fst = prefix + '_{0:0>3}Y_{1}_PeakAccel.png'
        fst = fst.format(m, zts)
        
        plt.suptitle(tst, fontsize=18)
        plt.savefig( fst, format='png', dpi=300, bbox_inches='tight')

    return

    
#=============================================================================

def stateq_3dof(fk, ukmax, MS, QS):

# 1. Sum up the modal contribution to top displacement
    
    N3  =  QS.shape[0]
    Peq =  np.zeros(N3)

    for k in range(3):
        Peq += ukmax[k]*((2*np.pi*fk[k])**2)*MS*QS[:,k]
 
    return Peq
    
   
#=============================================================================

def save_Scruton(prefix, zt, fk, QS, Bx, By):

    Sc      =  np.zeros((len(fk),4))
    Sc[:,0] =  1/np.sum(QS**2, axis=0)
    
    Sc[:,1] = (2*np.pi*zt)*Sc[:,0]/0.613/By/By
    Sc[:,2] = (2*np.pi*zt)*Sc[:,0]/0.613/Bx/Bx
    Sc[:,3] = (2*np.pi*zt)*Sc[:,0]/0.613/Bx/By

    Scrouton = pd.DataFrame(data   = Sc,
                            index  = range(len(fk)),
                            columns=['M_equiv',
                                     'bref^2 = By^2',
                                     'bref^2 = Bx^2',
                                     'bref^2 = Bx.By'])
        
    fst = prefix + '_Scrouton.xlsx'
    
    excelSc  = pd.ExcelWriter(fst)
    
    Scrouton.to_excel(excelSc,'Scruton')
    
    excelSc.save()
    
    return
    
    
#=============================================================================

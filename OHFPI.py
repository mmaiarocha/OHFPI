# -*- coding: utf-8 -*-

import sys
import pickle  as pk
import gzip    as gz
import numpy   as np
import pandas  as pd
import NBR6123

#==============================================================================
# 1. Loading table with all wind tunnel tests available
#==============================================================================

try:
    AllCases = pd.read_excel('AllCases.xlsx', sheet_name='AllCases')
except:
    sys.exit('File "AllCases.xlsx" not available in current folder!')

#==============================================================================
# 2. Preliminary test data conversion (from raw Scanivalve files)
#==============================================================================

def scanivalve_process(path, prefix, dwnd=0, avg=1, per=64*32):

# 1. Define data type for each column of files to be read           

# 1.1 Batch file
    typ1  = ['scan',   # name of Scanivalve file without extension
             'master', # master table file to be used
             'categ',  # far field roughness category 
             'block',  # wind tunnel blockage factor
             'Vm',     # wind tunnel mean speed
             'k0']     # wind tunnel k0 factor, referred to model top
    
# 1.2 Processed batch file
    typ2  = ['swnd',   # wind direction
             'categ',  # far field roughness category 
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

# 2. Preprocess all available Scanivalve files
        
    print('Processing "{0}" ... '.format(path))
    
    batchpickle =  path + prefix + '_batch.pk'
    batchfile   =  path + prefix + '_batch.txt'
    
    batchdata   =  pd.read_table(batchfile, names=typ1, header=None)
    batch       =  pd.DataFrame(columns=typ2)

    for ib, b in batchdata.iterrows():

# 2.1 Read Scanivalve file

        rec  =  path + prefix + '/' + b.scan   + '.DAT'
        mst  =  path +                b.master + '.txt'    
        
        wnd  =  int(rec[-7:-4]) - int(dwnd)  # wind direction from filename...
        if (wnd < 0): wnd = wnd + 360        # ... corrected and re-converted.
        swnd = '{0:0>3}'.format(wnd)   
        
        print('  Reading {0} deg from: "{1}"... '.format(swnd, rec))
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

            try:
                if t.zone not in [-1, -2]:
                    
                    ID1.append(list(ID).index(100*t.md1 + t.ch1))
                    ID2.append(list(ID).index(100*t.md2 + t.ch2))
                    
                    master = master.append(t, ignore_index=True)

                elif t.zone == -1:
                    IA1 = list(ID).index(100*t.md1 + t.ch1)
            
                elif t.zone == -2:
                    IA2 = list(ID).index(100*t.md1 + t.ch1)

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

        file =  path + prefix + '/' + prefix + '_' + swnd + '.gz'

        with gz.GzipFile(file, 'wb') as target:
            pk.dump((master, Cp), target)

# 2.6 Append data to batch dataframe

        new_b = pd.Series([swnd, b.categ, b.Vm, dPA, Td], index=typ2)
        batch = batch.append(new_b, ignore_index=True)

# 3. Dump and return batch dataframe
    
    with open(batchpickle, 'wb') as target:
        pk.dump(batch, target)

    print('... done!!!')
    return batch

#==============================================================================
# 3. OHFPI Class constructor
#==============================================================================

class OHFPI:
   
    def __init__(self, path, prefix, stype='3dof'):

        try:
            case  =  AllCases.loc[prefix]
        except:
            sys.exit('Entry not available in "AllCases.xlsx" file!')

        batchpickle = path + prefix + '_batch.pk'

        try:
            with open(batchpickle, 'rb') as target:
                self.batch = pk.load(target)
                
            self.path   = path
            self.prefix = prefix
            self.scale  = case['scale']
            self.Bx     = case['Bx']
            self.By     = case['By']
            self.Hz     = case['Hz']
            self.H0     = case['H0']
            self.H1     = case['H1']
            
        except:
            sys.exit('Could not read batch file {0}!'.format(batchpickle))

        self.coorfile = self.path + self.prefix + '_Coordinates.txt'
        self.massfile = self.path + self.prefix + '_Masses.txt'
        self.freqfile = self.path + self.prefix + '_Frequencies.txt'
        self.modefile = self.path + self.prefix + '_Modes.txt'

        try:
            if   stype.lower() == '3dof':
                self.read_3dof()
                
            elif stype.lower() == 'tqs':
                self.read_TQS()
                
            else:
                sys.exit('Type of structural system not recognized!')
            
        except:
            sys.exit('Could not read structural system files!')

        return


#==============================================================================
# 4. OHFPI Class methods
#==============================================================================
# 4.1. Read structure defined as 3 d.o.f. per floor
#==============================================================================

    def read_3dof(self):

# 1. Define data type for each column of files to be read           

# 1.1 Height and lumped masses
        typ1  = [('Si','i'),        # floor ID
                 ('ZS','f'),        # vertical coordinate of storey
                 ('MS','f'),        # storey mass 
                 ('IS','f')]        # storey moment of inertia (vertical axis)

# 1.2 Frequencies
        typ2  = [('Ki','i'),        # mode ID
                 ('fk','f')]        # modal frequency

# 1.3 Modal shapes
        typ3  = [('Si','i'),        # floor ID
                 ('qX','f'),        # modal displacement in X (translation)
                 ('qY','f'),        # modal displacement in Y (translation)
                 ('qR','f')]        # modal displacement in Z (rotation)

# 2. Build up lumped mass matrices

        dat1    =  np.genfromtxt(self.massfile, dtype=typ1)

        self.Si =  dat1['Si']
        self.NS =  len(self.Si)
        self.ZS =  dat1['ZS']

        X0 =  np.zeros(self.ZS.shape)

        self.XS =  np.vstack((      X0,         X0,    self.ZS  )).T
        self.MS =  np.vstack((dat1['MS'], dat1['MS'], dat1['IS'])).T.reshape((3*self.NS,))

# 3. Read modal frequencies (must be already sorted)

        dat2    =  np.genfromtxt(self.freqfile, dtype=typ2)
    
        self.fk =  dat2['fk']
        self.NQ =  len(self.fk)
        self.QS =  np.zeros((3*self.NS,self.NQ))

# 4. Read modal shapes (must be in the same order as fk)

        dat3 =  np.genfromtxt(self.modefile, dtype=typ3)
        
        Si0  =  dat3['Si'][0]            # to mark beginning of new modes
        k    = -1

        for s, qx, qy, qr in dat3:

            if (s  ==  Si0):
                k  =  k + 1
            
            i  = list(self.Si).index(s)  # locate store index by store ID

            self.QS[3*i+0, k] = qx
            self.QS[3*i+1, k] = qy
            self.QS[3*i+2, k] = qr

# 5. Normalize modal shapes for unitary modal masses
    
        for k in range(self.NQ):
            Mk = np.sum(self.MS*self.QS[:,k]*self.QS[:,k])
            self.QS[:,k] = self.QS[:,k] / np.sqrt(Mk)
        
# 6. Check orthogonality of modal shapes
    
        self.ortho = np.zeros((self.NQ,self.NQ))    

        for ii in range(self.NQ):
            for jj in range(self.NQ):
                self.ortho[ii,jj] = \
                np.dot(np.dot(self.QS[:,ii],np.diag(self.MS)),self.QS[:,jj])

        return

#==============================================================================
# 4.2. Read structure defined through TQS output files
#==============================================================================

    def read_TQS(self):
    
# 1. Define data type for each column of files to be read           

# 1.1 Nodal coordinates
        typ0  = [('Ni','i'),        # node ID
                 ('XN','f'),        #
                 ('YN','f'),        # node coordinates (designer ref. system)
                 ('ZN','f')]        #

# 1.2 Lumped nodal masses
        typ1  = [('Ni','i'),        # node ID
                 ('MN','f')]        # node mass 
    
# 1.3 Frequencies
        typ2  = [('Ki','i'),        # mode ID (not used)
                 ('fk','f')]        # modal frequency (must be sorted!!!)

# 1.4 Modal shapes
        typ3  = [('Ni','i'),        # node ID
                 ('qX','f'),        # modal displacement in X  (translation)
                 ('qY','f'),        # modal displacement in Y  (translation)
                 ('qZ','f')]        # modal displacement in Z  (translation)

# 2. Read nodal coordinates

        dat0    =  np.genfromtxt(self.coorfile, dtype=typ0)

        self.Ni =  dat0['Ni']
        self.NN =  len(self.Ni)
        self.XN =  np.vstack((dat0['XN'], dat0['YN'], dat0['ZN'])).T

# 3. Build up lumped mass matrices and lumped storey mass

        dat1    =  np.genfromtxt(self.massfile, dtype=typ1)
    
        self.Mi =  9810.*dat1['MN']
        self.MN =  np.vstack((self.Mi, self.Mi, self.Mi)).T.reshape((3*self.NN,))

# 4. Read modal frequencies (must be sorted)

        dat2    =  np.genfromtxt(self.freqfile, dtype=typ2)

        self.fk =  1/dat2['fk']
        self.NQ =  len(self.fk)
        self.QN =  np.zeros((3*self.NN,self.NQ))
    
# 5. Read modal shapes (must be ascending sorted by fk)

        dat3 =  np.genfromtxt(self.modefile, dtype=typ3)

        NiN  =  list(self.Ni)
        NiQ  =  list(dat3['Ni'][0:self.NN])    # all modes labels/order
        ni   =  np.empty((self.NN,), dtype='int')
    
        for ki, n in enumerate(NiN):           # labels/order from coord. file
            ni[ki] = NiQ.index(n)
    
        if len(ni) != len(np.unique(ni)):
            sys.exit('Node label not unique in modal shapes file!')
        
        for k in range(self.NQ):

            qX = dat3['qX'][k*self.NN:(k+1)*self.NN]
            qY = dat3['qY'][k*self.NN:(k+1)*self.NN]
            qZ = dat3['qZ'][k*self.NN:(k+1)*self.NN]

            self.QN[:,k] = np.vstack((qX[ni], qY[ni], qZ[ni])).T.reshape((3*self.NN,))       
        
# 6. Normalize modal shapes for unitary modal masses

        for k in range(self.NQ):
            Mk = np.sum(self.MN*self.QN[:,k]*self.QN[:,k])
            self.QN[:,k] = self.QN[:,k] / np.sqrt(Mk)

# 7. Check orthogonality of modal shapes
    
        ortho = np.zeros((self.NQ,self.NQ))    

        for ii in range(self.NQ):
            for jj in range(self.NQ):
                ortho[ii,jj] = np.dot(self.QN[:,ii],self.QN[:,jj])

        return

#==============================================================================
# 4.3. Object iterator over all wind directions
#==============================================================================

    def batchIter(self):

        for kb, b in self.batch.iterrows():

            file =  self.path + self.prefix + '/' + \
                                self.prefix + '_' + b.swnd + '.gz'  
                                
            try:
                with gz.GzipFile(file, 'rb') as target:
                    master, Cp = pk.load(target)
                    yield kb, b, master, Cp
            
            except:
                sys.exit('Could not read data file {0}!'.format(file))

#==============================================================================
# 4.4. Aerodynamic coefficients
#==============================================================================
    def aerocoeffs(self):
    
        Cm     =  np.empty((3,len(self.batch)))
        Cs     =  np.empty((3,len(self.batch)))
    
        Cx, Cy =  NBR6123.drag(self.Bx, self.By, self.Hz, case='low')
    
        if self.Bx > self.By: Ct = 0.075*Cy*self.Bx/self.By
        else:                 Ct = 0.075*Cx*self.By/self.Bx

        for kb, b, master, Cp in self.batchIter():

            A    = master['A'].values
            XT   = master[['X', 'Y', 'Z' ]].values
            CT   = master[['cx','cy','cz']].values

            CAx =  A*CT[:,0]
            CAy =  A*CT[:,1]
            CAt =  XT[:,0]*CAy - XT[:,1]*CAx

            Cx0 = -Cp.dot(CAx)/(self.By*        (self.Hz-self.H0))
            Cy0 = -Cp.dot(CAy)/(self.Bx*        (self.Hz-self.H0))
            Ct0 = -Cp.dot(CAt)/(self.By*self.Bx*(self.Hz-self.H0))

            Cm[0,kb] = Cx0.mean()
            Cm[1,kb] = Cy0.mean()
            Cm[2,kb] = Ct0.mean()
        
            Cs[0,kb] = Cx0.std()
            Cs[1,kb] = Cy0.std()
            Cs[2,kb] = Ct0.std()
            
        return Cx, Cy, Ct, Cm, Cs

#==============================================================================


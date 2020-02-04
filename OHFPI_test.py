# -*- coding: utf-8 -*-

import OHFPI

# E01: CAARC
'''
path   = 'Data/E01/'
prefix = 'E01_IS023'
OHFPI.scanivalve_process(path, prefix)
E01    = OHFPI.HFPI(path, prefix, stype='3dof')
'''

# E24: Campanário da Basílica do Divino Pai Eterno

path   = 'Data/E24/'
prefix = 'E24_Conf1'
OHFPI.scanivalve_process(path, prefix)

E24    = OHFPI.HFPI(path, prefix, stype='tqs')



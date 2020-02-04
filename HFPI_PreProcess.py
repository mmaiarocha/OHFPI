# -*- coding: utf-8 -*-

import HFPI0

#==============================================================================
# Process all Scanivalve files
#==============================================================================
'''
dirname = 'Data/E01/' # CAARC

prefix  = 'E01_IS023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E02/' # Rossi Curitiba
 
prefix  = 'E02_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E02_CV034'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E02_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E02_SV034'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E03/' # Torre Caelum

prefix  = 'E03_CV011'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E03_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E03_SV011'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E03_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E04/' # W-Torre Morumbi

prefix  = 'E04a_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E04b_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E05/' # Fotografia

prefix  = 'E05a_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E05a_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E05b_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E05b_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E06/' # Carlos Steinen

prefix  = 'E06_CM023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E06_IS023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E07/' # Yuni Santos

prefix  = 'E07_CV011'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E07_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E07_CV034'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E07_SV011'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E07_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E07_SV034'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E08/' # Horizonte

prefix  = 'E08_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E08_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E09/' # Multiplan Office

prefix  = 'E09_Mixed'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E10/' # Multiplan Homestay

prefix  = 'E10_Mixed'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E11/' # Autentico Belem

prefix  = 'E11_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E12/' # Legend

prefix  = 'E12a_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E12a_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E12b_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E12b_SV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E13/' # FG Infinity Coast
prefix  = 'E13_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E13_Conf2'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E14/' # VR Resort Caxias

prefix  = 'E14_CV023'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E15/' # VR Faria Lima

prefix  = 'E15_CM023_22Hz'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E15_CM023_30Hz'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E16/' # Setin Rio Branco

prefix  = 'E16_CV034_200'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E16_CV034_400'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E16_IS034_200'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E16_IS034_400'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E17/' # Pininfarina Yachthouse

prefix  = 'E17a_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E17a_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E17b_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E17b_Conf2'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E18/' # Embraed New York Tower

prefix  = 'E18_Mixed'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
 
#==============================================================================

dirname = 'Data/E19/' # Melnick Parque do Pontal

prefix  = 'E19_Mixed'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E20/' # Tishman Speyer Pátio da Marítima

prefix  = 'E20_F1Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E20_F1Conf2'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E20_F2Conf2'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E21/' # Multiplan Vila Hípica

prefix  = 'E21_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
prefix  = 'E21_Conf2'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E22/' # Winter Gomes Portinho

prefix  = 'E22_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

dirname = 'Data/E23/' # GPL-HFA Goiânia

prefix  = 'E23_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)
'''
#==============================================================================

dirname = 'Data/E24/' # Campanário da Basílica do Divino Pai Eterno

prefix  = 'E24_Conf1'
batch   =  HFPI0.scanivalve_process(dirname, prefix)

#==============================================================================

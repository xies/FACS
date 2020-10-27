#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 12:12:06 2020

@author: xies
"""

import numpy as np
import seaborn as sb
import pandas as pd
import matplotlib.path as mplPath
import matplotlib.pyplot as plt
from FlowCytometryTools import FCMeasurement
from roipoly import roipoly
from scipy import stats

sb.set()


#%%

# HMEC
filename = '/Users/xies/Box/Others/EZ_HMEC_facs/RB-Flow-Mimi/HMEC-16Aug17/specimen_HMEC-wt+II_Abs_only-AF488+AF594+DAPI_011.fcs'
hmec_bg = FCMeasurement(ID='Test Sample', datafile=filename)

dapi = hmec_bg['DAPI-A']
fsc = hmec_bg['FSC-A']
plt.scatter(dapi,fsc,alpha=0.01)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((hmec_bg['DAPI-A'],hmec_bg['FSC-A'])).T )
bg_g1 = hmec_bg[g1_gateI]
size_bins_g1 = np.linspace(bg_g1['FSC-A'].min(),bg_g1['FSC-A'].max(),num = 26)

plt.scatter(dapi,fsc,alpha=0.01)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((hmec_bg['DAPI-A'],hmec_bg['FSC-A'])).T )
bg_g2 = hmec_bg[g2_gateI]
size_bins_g2 = np.linspace(bg_g2['FSC-A'].min(),bg_g2['FSC-A'].max(),num = 26)

bg_rb_g1,_ = get_bin_means(bg_g1['FSC-A'],bg_g1['Alexa Fluor 488-A'],size_bins_g1)
bg_rb_g2,_ = get_bin_means(bg_g2['FSC-A'],bg_g2['Alexa Fluor 488-A'],size_bins_g2)
bg_actin_g1,_ = get_bin_means(bg_g1['FSC-A'],bg_g1['PE-Texas Red-A'],size_bins_g1)
bg_actin_g2,_ = get_bin_means(bg_g2['FSC-A'],bg_g2['PE-Texas Red-A'],size_bins_g2)



filename = '/Users/xies/Box/Others/EZ_HMEC_facs/RB-Flow-Mimi/HMEC-16Aug17/specimen_HMEC-wt+RB-AF488+bActin-AF594+DAPI_002.fcs'
hmec = FCMeasurement(ID='Test Sample', datafile=filename)

dapi = hmec['DAPI-A']
fsc = hmec['FSC-A']
plt.scatter(dapi,fsc,alpha=0.01)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((hmec['DAPI-A'],hmec['FSC-A'])).T )
hmec_g1 = hmec[g1_gateI]

plt.scatter(dapi,fsc,alpha=0.01)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((hmec['DAPI-A'],hmec['FSC-A'])).T )
hmec_g2 = hmec[g2_gateI]

rb_g1,_ = get_bin_means(hmec_g1['FSC-A'],hmec_g1['Alexa Fluor 488-A'],size_bins_g1)
actin_g1,_ = get_bin_means(hmec_g1['FSC-A'],hmec_g1['PE-Texas Red-A'],size_bins_g1)
rb_g2,_ = get_bin_means(hmec_g2['FSC-A'],hmec_g2['Alexa Fluor 488-A'],size_bins_g2)
actin_g2,_ = get_bin_means(hmec_g2['FSC-A'],hmec_g2['PE-Texas Red-A'],size_bins_g2)
sizes_g1 = np.array( [ (x+y)/2 for x,y in zip( size_bins_g1[0:-1],size_bins_g1[1:] ) ] )
sizes_g2 = np.array( [ (x+y)/2 for x,y in zip( size_bins_g2[0:-1],size_bins_g2[1:] ) ] )


plt.plot(standardize(sizes_g1),standardize(rb_g1 - bg_rb_g1)), plt.xlim([0,2]); plt.ylim([0,2])
plt.plot(standardize(sizes_g2),standardize(rb_g2 - bg_rb_g2)), plt.xlim([0,2]); plt.ylim([0,2])



plt.plot(standardize(sizes_g1),standardize(actin_g1 - bg_actin_g1)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.plot(standardize(sizes_g2),standardize(actin_g2 - bg_actin_g2)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.xlabel('FSC')
plt.ylabel('Protein total')

plt.legend(['Rb - G1','Rb - G2','Actin - G1','Actin - G2'])


#%%

# HMEC
filename = '/Users/xies/Box/Others/EZ_HMEC_facs/RB-Flow-Mimi/HLF-25Mar16/specimen_HLF-wt-fix-II_Abs_only-Mouse-TexRed+Rab-AF488+Hoechst_009.fcs'
hlf_bg = FCMeasurement(ID='Test Sample', datafile=filename)

dapi = hlf_bg['DAPI-A']
fsc = hlf_bg['FSC-A']
plt.scatter(dapi,fsc,alpha=0.1)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((hlf_bg['DAPI-A'],hlf_bg['FSC-A'])).T )
bg_g1 = hlf_bg[g1_gateI]
size_bins_g1 = np.linspace(bg_g1['FSC-A'].min(),bg_g1['FSC-A'].max(),num = 26)

plt.scatter(dapi,fsc,alpha=0.1)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((hlf_bg['DAPI-A'],hlf_bg['FSC-A'])).T )
bg_g2 = hlf_bg[g2_gateI]
size_bins_g2 = np.linspace(bg_g2['FSC-A'].min(),bg_g2['FSC-A'].max(),num = 26)

bg_rb_g1,_ = get_bin_means(bg_g1['FSC-A'],bg_g1['Alexa Fluor 488-A'],size_bins_g1,minimum_n=2)
bg_rb_g2,_ = get_bin_means(bg_g2['FSC-A'],bg_g2['Alexa Fluor 488-A'],size_bins_g2,minimum_n=2)
bg_actin_g1,_ = get_bin_means(bg_g1['FSC-A'],bg_g1['PE-Texas Red-A'],size_bins_g1,minimum_n=2)
bg_actin_g2,_ = get_bin_means(bg_g2['FSC-A'],bg_g2['PE-Texas Red-A'],size_bins_g2,minimum_n=2)



filename = '/Users/xies/Box/Others/EZ_HMEC_facs/RB-Flow-Mimi/HLF-25Mar16/specimen_HLF-wt-fix-RB-TexRed+Actin-AF488+Hoechst_007.fcs'
hlf = FCMeasurement(ID='Test Sample', datafile=filename)

dapi = hlf['DAPI-A']
fsc = hlf['FSC-A']
plt.scatter(dapi,fsc,alpha=0.1)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((hlf['DAPI-A'],hlf['FSC-A'])).T )
hlf_g1 = hlf[g1_gateI]

plt.scatter(dapi,fsc,alpha=0.1)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((hlf['DAPI-A'],hlf['FSC-A'])).T )
hlf_g2 = hlf[g2_gateI]

rb_g1,_ = get_bin_means(hlf_g1['FSC-A'],hlf_g1['Alexa Fluor 488-A'],size_bins_g1,minimum_n=3)
actin_g1,_ = get_bin_means(hlf_g1['FSC-A'],hlf_g1['PE-Texas Red-A'],size_bins_g1,minimum_n=3)
rb_g2,_ = get_bin_means(hlf_g2['FSC-A'],hlf_g2['Alexa Fluor 488-A'],size_bins_g2,minimum_n=3)
actin_g2,_ = get_bin_means(hlf_g2['FSC-A'],hlf_g2['PE-Texas Red-A'],size_bins_g2,minimum_n=3)
sizes_g1 = np.array( [ (x+y)/2 for x,y in zip( size_bins_g1[0:-1],size_bins_g1[1:] ) ] )
sizes_g2 = np.array( [ (x+y)/2 for x,y in zip( size_bins_g2[0:-1],size_bins_g2[1:] ) ] )



plt.plot(standardize(sizes_g1),standardize(rb_g1 - bg_rb_g1)), plt.xlim([0,2]); plt.ylim([0,2])
plt.plot(standardize(sizes_g2),standardize(rb_g2 - bg_rb_g2)), plt.xlim([0,2]); plt.ylim([0,2])

plt.plot(standardize(sizes_g1),standardize(actin_g1 - bg_actin_g1)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.plot(standardize(sizes_g2),standardize(actin_g2 - bg_actin_g2)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.xlabel('FSC')
plt.ylabel('Protein total')

plt.legend(['Rb - G1','Rb - G2','Actin - G1','Actin - G2'])






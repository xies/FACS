plt.#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 21 21:36:36 2020

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

# T98G
filename = '/Users/xies/Box/Others/SZ_FACS_RB/SZ-091320 Async OPP_Group_T98G OPP647 2nd488.fcs'
tg_bg = FCMeasurement(ID='Test Sample', datafile=filename)

size_bins = np.linspace(tg_bg['FSC-A'].min(),tg_bg['FSC-A'].max(),num = 26)
bg_rb = get_bin_means(tg_bg['FSC-A'],tg_bg['BL1-A'],size_bins)

filename = '/Users/xies/Box/Others/SZ_FACS_RB/SZ-091320 Async OPP_Group_T98G OPP647 RB488.fcs'
tg = FCMeasurement(ID='Test Sample', datafile=filename)

dapi = tg['VL1-A']
fsc = tg['FSC-A']
plt.scatter(dapi,fsc,alpha=0.002)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((tg['VL1-A'],tg['FSC-A'])).T )
tg_g1 = tg[g1_gateI]

plt.scatter(dapi,fsc,alpha=0.002)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((tg['VL1-A'],tg['FSC-A'])).T )
tg_g2 = tg[g2_gateI]

rb_g1,_ = get_bin_means(tg_g1['FSC-A'],tg_g1['BL1-A'],size_bins)
rb_g2,_ = get_bin_means(tg_g2['FSC-A'],tg_g2['BL1-A'],size_bins)
sizes = np.array( [ (x+y)/2 for x,y in zip( size_bins[0:-1],size_bins[1:] ) ] )


plt.plot(sizes,rb_g1 - bg_rb_g1), plt.xlim([0,2]); plt.ylim([0,2])
plt.plot(sizes,rb_g2 - bg_rb_g2), plt.xlim([0,2]); plt.ylim([0,2])



plt.plot(standardize(sizes_g1),standardize(actin_g1 - bg_actin_g1)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.plot(standardize(sizes_g2),standardize(actin_g2 - bg_actin_g2)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.xlabel('FSC')
plt.ylabel('Protein total')

plt.legend(['Rb - G1','Rb - G2','Actin - G1','Actin - G2'])



#%%


# HMEC
filename = '/Users/xies/Box/Others/SZ_FACS_RB/SZ-091320 Async OPP_Group_HMEC OPP647 2nd488.fcs'
hmec_bg = FCMeasurement(ID='Test Sample', datafile=filename)

dapi = hmec_bg['VL1-A']
fsc = hmec_bg['FSC-A']
plt.scatter(dapi,fsc,alpha=0.002)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((hmec_bg['VL1-A'],hmec_bg['FSC-A'])).T )
bg_g1 = hmec_bg[g1_gateI]
size_bins_g1 = np.linspace(bg_g1['FSC-A'].min(),bg_g1['FSC-A'].max(),num = 26)

plt.scatter(dapi,fsc,alpha=0.002)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((hmec_bg['VL1-A'],hmec_bg['FSC-A'])).T )
bg_g2 = hmec_bg[g2_gateI]
size_bins_g2 = np.linspace(bg_g2['FSC-A'].min(),bg_g2['FSC-A'].max(),num = 26)

bg_rb_g1,_ = get_bin_means(bg_g1['FSC-A'],bg_g1['BL1-A'],size_bins_g1)
bg_rb_g2,_ = get_bin_means(bg_g2['FSC-A'],bg_g2['BL1-A'],size_bins_g2)
bg_actin_g1,_ = get_bin_means(bg_g1['FSC-A'],bg_g1['RL1-A'],size_bins_g1)
bg_actin_g2,_ = get_bin_means(bg_g2['FSC-A'],bg_g2['RL1-A'],size_bins_g2)



filename = '/Users/xies/Box/Others/SZ_FACS_RB/SZ-091320 Async OPP_Group_HMEC OPP647 RB488.fcs'
hmec = FCMeasurement(ID='Test Sample', datafile=filename)


dapi = hmec['VL1-A']
fsc = hmec['FSC-A']
plt.scatter(dapi,fsc,alpha=0.002)
g1_gate = roipoly(color='r')
g1_gate_p = mplPath.Path( np.vstack((g1_gate.x,g1_gate.y)).T )
g1_gateI = g1_gate_p.contains_points( np.vstack((hmec['VL1-A'],hmec['FSC-A'])).T )
hmec_g1 = hmec[g1_gateI]

plt.scatter(dapi,fsc,alpha=0.002)
g2_gate = roipoly(color='r')
g2_gate_p = mplPath.Path( np.vstack((g2_gate.x,g2_gate.y)).T )
g2_gateI = g2_gate_p.contains_points( np.vstack((hmec['VL1-A'],hmec['FSC-A'])).T )
hmec_g2 = hmec[g2_gateI]

rb_g1,_ = get_bin_means(hmec_g1['FSC-A'],hmec_g1['BL1-A'],size_bins_g1)
actin_g1,_ = get_bin_means(hmec_g1['FSC-A'],hmec_g1['RL1-A'],size_bins_g1)
rb_g2,_ = get_bin_means(hmec_g2['FSC-A'],hmec_g2['BL1-A'],size_bins_g2)
actin_g2,_ = get_bin_means(hmec_g2['FSC-A'],hmec_g2['RL1-A'],size_bins_g2)
sizes_g1 = np.array( [ (x+y)/2 for x,y in zip( size_bins_g1[0:-1],size_bins_g1[1:] ) ] )
sizes_g2 = np.array( [ (x+y)/2 for x,y in zip( size_bins_g2[0:-1],size_bins_g2[1:] ) ] )




plt.plot(sizes_g1,rb_g1 - bg_rb_g1), plt.xlim([0,2]); plt.ylim([0,2])
plt.plot(standardize(sizes_g2),rb_g2 - bg_rb_g2), plt.xlim([0,2]); plt.ylim([0,2])



plt.plot(standardize(sizes_g1),standardize(actin_g1 - bg_actin_g1)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.plot(standardize(sizes_g2),standardize(actin_g2 - bg_actin_g2)), plt.xlim([0,2.5]); plt.ylim([0,2.5])
plt.xlabel('FSC')
plt.ylabel('Protein total')

plt.legend(['Rb - G1','Rb - G2','Actin - G1','Actin - G2'])



#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 14:23:51 2018

@author: mimi
"""

import seaborn as sb
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

sb.set()
plt.style.use('ggplot')

filename = '/Users/mimi/Box Sync/mIOs/FACS data/11-20-18 Lgr5 500nM Palbo 6hr 24hr/cell_cycle.xlsx'
df = pd.read_excel(filename)
df['Condition Time'] = df['Time'].map(str) + ' ' + df['Condition'].map(str)

iscs = df[ df['Cell Type'] == 'Ki67+ Lgr5+']
tas = df[ df['Cell Type'] == 'Ki67+ Lgr5-']


iscs.plot(kind='bar',stacked=True,x='Condition Time')
plt.title('ISCs: Ki67+ Lgr5+')


tas.plot(kind='bar',stacked=True,x='Condition Time')
plt.title('TAs: Ki67+ Lgr5-')


plt.close('all')

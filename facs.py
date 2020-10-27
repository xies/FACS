#! /usr/bin/env python
"""
@author: devoncb@stanford.edu, xies@stanford.edu
"""

import numpy as np
import pandas as pd
from scipy import stats
from FlowCytometryTools import FCMeasurement
from FlowCytometryTools import PolyGate
from FlowCytometryTools import ThresholdGate
import glob
import pickle

def load_facs(file_name):
    """
    Load .fcs data and converts into a pandas DataFrame
    
    Parameters
    ----------
        file_name - path string to the .fjo file
        
    Returns
    -------
        df- pandas DataFrame
    
    """
    sample = FCMeasurement(ID='Test Sample', datafile=file_name)
    parameters = list(sample.channel_names)
    for par in range(len(parameters)):
    	parameters[par] = parameters[par].encode('ascii','ignore')
	raw_data = sample.data.values
    # sample_id = sample.meta['SampleID'].encode('ascii','ignore')
    
    return raw_data, np.array(parameters)

"""
FUNCTION: Gates cells based on specified threshold in requested parameters;
general function designed to take any number of requested parameters and provide
universal gating
"""

def gate_cells(data, parameters, gate_pars, gate_thr, cut_max = False, \
	cut_pars = None):
	gate_data = data
	
	for gates in range(len(gate_pars)):
		par = NP.where(NP.array(parameters) == NP.array(gate_pars[gates]))[0][0]
		ind = NP.where(gate_data[:,par] > gate_thr[gates])[0]
		gate_data = gate_data[ind,:]
		
	if cut_max == True:
		for param in cut_pars:
			par = NP.where(NP.array(parameters) == \
				NP.array(param))[0][0]
			ind = NP.where(gate_data[:,par] < NP.max(gate_data[:,par]))[0]
			gate_data = gate_data[ind,:]
			
	return gate_data

"""
FUNCTION: Loads and filters all data in folder.
Finds all files in folder, filters junk cells, separates into clusters and 
concatenates into dictionary with callable gene name.
"""

def filter_all_data():
	all_data = {}
	strain_list = []
	file_list = glob.glob('*.fcs')
	
	raw_data, name, param = load_file(file_list[0])
	
	print param
	print ""
	gate_pars = input("Enter gating parameters as a list: ")
	gate_thr = input("Enter gating thresholds as a list: ")
	max_treat = input("Enter whether max should be cut (True or False): ")
	if max_treat == True:
		cut_pars = input("Enter parameters to remove max: ")
	else:
		cut_pars = None
	
	for file in file_list:
		raw_data, name, parameters = load_file(file)
		inc_data = raw_input("Should {0} be included?: ".format(name))
		if inc_data == 'y':
			filtered_data = gate_cells(raw_data, parameters, gate_pars, \
				gate_thr, cut_max = max_treat, cut_pars = cut_pars)
			all_data[name] = \
				{'strain':name,'parameters':parameters,'data':filtered_data}

	return all_data

"""
FUNCTION: Perform background subtraction based on fits to data, will require
input of background strain - prompted from wrapper
"""

def bg_subtract(all_data, bg_strain, x_param, y_param):
	pos_bg_x = NP.where(all_data[bg_strain]['parameters']==x_param)[0][0]
	pos_bg_y = NP.where(all_data[bg_strain]['parameters']==y_param)[0][0]
	bg_x = all_data[bg_strain]['data'][:,pos_bg_x]	
	bg_y= all_data[bg_strain]['data'][:,pos_bg_y]
	bg_reg = stats.linregress(bg_x,bg_y)
	bg_func = NP.poly1d(NP.array([bg_reg[0],bg_reg[1]]))

	for strain in all_data.keys():
		pos_x = NP.where(all_data[strain]['parameters']==x_param)[0][0]
		pos_y = NP.where(all_data[strain]['parameters']==y_param)[0][0]
		bg_x = all_data[strain]['data'][:,pos_x]	
		bg_y= all_data[strain]['data'][:,pos_y]
		bg_sub_y = bg_y - bg_func(bg_x)
		bg_sub_y = bg_sub_y[:,None]
		all_data[strain]['data'] = \
			NP.hstack((all_data[strain]['data'],bg_sub_y))
		all_data[strain]['parameters'] = \
			NP.append(all_data[strain]['parameters'],y_param + '_bg_sub')
	
	return all_data

"""
FUNCTION: Bin means
"""

def bin_means(x, y, n_bin=100, plot_range = None):
	if plot_range == None:
		bins = NP.linspace(min(x),max(x),n_bin)
	else:
		bins = NP.linspace(plot_range[0],plot_range[1],n_bin)
	plot_data = NP.vstack((x,y)).transpose()
	binned_data = []
	for bin in range(len(bins)-1):
		temp_bin = plot_data[NP.where(plot_data[:,0] >= bins[bin])[0],:]
		temp_bin2 = temp_bin[NP.where(temp_bin[:,0] < bins[bin + 1]),:]
		if len(temp_bin2[0][:,0]) > 0:
			mean = NP.mean(temp_bin2[0][:,1])
			sd = NP.std(temp_bin2[0][:,1])
			binned_data.append([NP.mean([bins[bin], bins[bin + 1]]),mean, sd])
	return NP.array(binned_data)

def bin_all(all_data, x_param, y_param):
	for sample in all_data.keys():
		temp_data = all_data[sample]['data']
		temp_param = all_data[sample]['parameters']
		pos_x = NP.where(all_data[sample]['parameters']==x_param)[0][0]
		pos_y = NP.where(all_data[sample]['parameters']==y_param)[0][0]
		x_data = temp_data[:,pos_x]
		y_data = temp_data[:,pos_y]
		binned_data = bin_means(x_data,y_data)
		all_data[sample]['bin_data'] = binned_data
	
	return all_data

"""
FUNCTION: Plot data plus binned means
"""

def plot_scaling_data(all_data, x_param, y_param, file_type = 'pdf'):
	for sample in all_data.keys():
		temp_data = all_data[sample]['data']
		temp_bin = all_data[sample]['bin_data']
		temp_param = all_data[sample]['parameters']
		pos_x = NP.where(all_data[sample]['parameters']==x_param)[0][0]
		pos_y = NP.where(all_data[sample]['parameters']==y_param)[0][0]
		x_data = temp_data[:,pos_x]
		y_data = temp_data[:,pos_y]
		PL.plot(x_data, y_data, '.', markersize = 12, \
			alpha = .05, label = "{0}".format(sample))
		PL.plot(temp_bin[:,0],temp_bin[:,1],'-k', linewidth = 2)
		PL.legend(loc = "upper left")
		PL.xlabel(x_param + ' (Size)')
		PL.ylabel("Protein Abundance (AU)")
		PL.ylim(0,NP.max(y_data))
		PL.xlim(0,NP.max(x_data))
		savepath = '{0}.{1}'.format(sample,file_type)
		PL.savefig(savepath)
		PL.close()

"""
FUNCTION: Designed to permit manual graph making without typing all lines
"""

def plot_scaling_manual(all_data, file_type = 'tif'):
	print "All available strains: " 
	print all_data.keys()
	sample = raw_input("Enter strain name: ")
	print "Parameters: "
	print all_data[sample]['parameters']
	x_param = raw_input('Enter size parameter: ')
	y_param = raw_input('Enter protein abundance parameter: ')
	x_max = NP.array(input("Enter x axis maximum: "))
	y_max = NP.array(input("Enter y axis maximum: "))
	temp_data = all_data[sample]['data']
	temp_bin = all_data[sample]['bin_data']
	temp_param = all_data[sample]['parameters']
	pos_x = NP.where(all_data[sample]['parameters']==x_param)[0][0]
	pos_y = NP.where(all_data[sample]['parameters']==y_param)[0][0]
	x_data = temp_data[:,pos_x]
	y_data = temp_data[:,pos_y]
	PL.plot(x_data, y_data, '.', markersize = 12, \
		alpha = .05, label = "{0}".format(sample))
	PL.plot(temp_bin[:,0],temp_bin[:,1],'-k', linewidth = 2)
	PL.legend(loc = "upper left")
	PL.xlabel(x_param + ' (Size)')
	PL.ylabel("Protein Abundance (AU)")
	PL.ylim(0,y_max)
	PL.xlim(0,x_max)
	savepath = raw_input("Enter plot file name: ")
	PL.savefig(savepath + '.' + file_type)
	PL.close()

"""
FUNCTION: Wrapper for all procesing
"""

def analyze_data():
	print "Analyzing current folder..."
	print "Gating data..."
	all_data = filter_all_data()
	all_strains = all_data.keys()
	print all_data[all_strains[0]]['parameters']
	x_param = raw_input('Enter size parameter: ')
	y_param = raw_input('Enter protein abundance parameter: ')
	print "Current folder data: {0}".format(all_strains)
	bg_strain = raw_input('Enter background strain: ')
	print "Performing background subtraction..."
	all_data = bg_subtract(all_data, bg_strain, x_param, y_param)
	print "Binning data..."
	data = bin_all(all_data, x_param, y_param + '_bg_sub')
	print "Plotting all data..."
	plot_scaling_data(data, x_param, y_param + '_bg_sub', file_type = "tif")
	print "Saving all data..."
	file_name = raw_input("Enter file name: ")
	pickle.dump(data, open(file_name + ".p", "w"))
	print "Process complete."
	
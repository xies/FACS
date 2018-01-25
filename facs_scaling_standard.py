#! /usr/bin/env python

"""
Script to analyze a folder of facs data run on Kopito facscalibur for scaling analysis.
"""

import fcm as fcm
print "Imported fcm..."
import numpy as NP
print "Imported NumPy..."
import pylab as PL
print "Imported PyLab..."
import scipy as SP
print "Imported SciPy..."
import glob as glob
print "Glob imported..."
import pickle
print "Pickle imported..."

"""
FUNCTION: Loads facs data and converts to numpy array
Only useful for GFP only flow data.  Must be updated if multi-colored analysis to be
performed.  Also, currently uses local folder or full path as given.
"""

def load_file(File_name):
	Path = File_name
	current_data = fcm.loadFCS(Path)
	sample_id = current_data.notes['text']['sample id']
	FSC = current_data[:,0]
	SSC = current_data[:,1]
	GFP = current_data[:,2]
	data = NP.vstack((FSC, SSC, GFP))
	data = NP.transpose(data)
	
#	print "Loaded {0}...".format(File_name)

	return data, sample_id

"""
FUNCTION: Filter out cell debris and saturated values. written by Marie Loustalot.
Intakes data from 
 and removes values below a pre-defined hard threshold.  Also,
removes values =1023 because this is saturated.
"""

def gate_cells(cells,fsc_thr,ssc_thr,cut_max = True, cut_GFP_max = True):
	
	ind = NP.where(cells[:,0]>fsc_thr)
	cells2 = cells[ind[0],:]
	ind2 = NP.where(cells2[:,1] > ssc_thr)
	cells3 = cells2[ind2[0],:]

	if cut_max == True:
		ind3 = NP.where(cells3[:,1] < NP.max(cells3[:,1]))
		cells4 = cells3[ind3[0],:]
	else:
		cells4 = cells3

	if cut_GFP_max == True:
		ind4 = NP.where(cells4[:,2] < NP.max(cells4[:,2]))
		cells5 = cells4[ind4[0],:]
	else:
		cells5 = cells4

	return cells5

"""
FUNCTION: Loads and filters all data in folder.
Finds all files in folder, filters junk cells, separates into clusters and concatenates
into dictionary with callable gene name.
"""

def filter_all_data(fsc_thr, ssc_thr, max_treat = True):
	all_data = {}
	strain_list = []
	file_list = glob.glob('Data*')

	for file in file_list:
		raw_data, name = load_file(file)
		strain_list.append(name)

	strain_list = NP.unique(NP.array(strain_list))

	for file in file_list:
		raw_data, name = load_file(file)
		filtered_data = gate_cells(raw_data, fsc_thr, ssc_thr, cut_max = max_treat)
		all_data[name] = filtered_data

	return all_data, strain_list

"""
FUNCTION: Bin data by size.
Takes all data and bins across range of values.  Constant size bins enforced to permit
background subtraction bin-wise.
"""

def plot_size_bins(x, y, n_bin=100, plot_range = [0, 1023]):
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

def bin_all(data):
	binned_data = {}
	for strain in data.keys():
		x = data[strain][:,1]
		y = data[strain][:,2]
		temp = plot_size_bins(x,y)
		binned_data[strain] = temp
	
	return binned_data

"""
FUNCTION: Background subtract data
Uses binned means to do background subtraction - selected over fitting as it reduces the
parameters of the final plot.  Should be more robust to weird data.  Fixed bin sizes in 
plot_size_bins() make this feasible.
"""

def bg_sub_wrap(raw_data, binned_data, bg_strain, strain_list):
	if type(bg_strain) == str:
		bg_strain = [bg_strain] * len(strain_list)
	if type(bg_strain) == list:
		bg_strain = bg_strain
	
	new_raw_data = {}
	new_binned_data = {}
	for sample in range(len(strain_list)):
		new_temp_binned, new_temp_raw = bg_sub(raw_data, binned_data, bg_strain, sample, \
			strain_list)
		new_binned_data[strain_list[sample]] = new_temp_binned
		new_raw_data[strain_list[sample]] = new_temp_raw
	
	return new_binned_data, new_raw_data

def bg_sub(raw_data, binned_data, bg_strain, sample, strain_list):
	strain = strain_list[sample]
	print strain
	bg = bg_strain[sample]
	print bg
	bg_bin_fit = bg_fit = NP.polyfit(raw_data[bg][:,1],raw_data[bg][:,2],1)
	bg_function = NP.poly1d(bg_fit)
	bg_sub_val = binned_data[strain][:,1] - bg_function(binned_data[strain][:,0])
	all_vals = NP.vstack([NP.transpose(binned_data[strain]), bg_sub_val]).transpose()
	
	raw_sub_val = NP.zeros(len(raw_data[strain][:,1]))
	bg_bin_fit = bg_fit = NP.polyfit(raw_data[bg][:,1],raw_data[bg][:,2],1)
	bg_function = NP.poly1d(bg_fit)
	sample_data = raw_data[strain]
	corrected_GFP = sample_data[:,2] - bg_function(sample_data[:,1])
	corrected_data = \
		NP.vstack((sample_data.transpose(),corrected_GFP)).transpose()
	all_vals_raw = corrected_data
	
	return all_vals, all_vals_raw

"""
FUNCTION: Plot background subtracted graphs.
Takes all strains in current data from binned data and plots them individually.  Saves
plots as PDF.
"""

def plot_save_data(proc_data, binned_data, file_type = 'pdf'):
	for sample in proc_data.keys():
		temp_data = proc_data[sample]
		temp_bin = binned_data[sample]
		PL.plot(proc_data[sample][:,1], proc_data[sample][:,3], '.', markersize = 12, \
			alpha = .05, label = "{0}".format(sample))
		PL.plot(temp_bin[:,0],temp_bin[:,3],'-k', linewidth = 2)
		PL.legend(loc = "upper left")
		PL.xlabel("SSC (Volume)")
		PL.ylabel("GFP Intensity (AU)")
		PL.ylim(0,1050)
		PL.xlim(0,1050)
		savepath = '{0}.{1}'.format(sample,file_type)
		PL.savefig(savepath)
		PL.close()

"""
FUNCTION: Core wrapper, call from interactive python, called by _main_ function.
"""

def analyze_data(fsc_thr, ssc_thr, bg_strain):
	print "Analyzing current folder..."
	print "Gating data..."
	data, strains = filter_all_data(fsc_thr, ssc_thr)
	print "Current folder data: {0}".format(strains)
	print "Binning data..."
	bin_data = bin_all(data)
	print "Performing background subtraction..."
	new_bin, new_raw = bg_sub_wrap(data, bin_data, bg_strain, strains)
	print "Plotting all data..."
	plot_save_data(new_raw, new_bin, file_type = "tif")
	print "Saving all data..."
	pickle.dump(new_raw, open("raw_data.p", "w"))
	pickle.dump(new_bin, open("binned_data.p", "w"))
	print "Process complete."

"""
FUNCTION: Double plotting function.
Plots haploid and diploid strains with trend line for comparison.  Runs independently of
analyze_data() but depends on having the final analyzed data in the current folder.
"""

def plot_dbl_graph(gene_name, file_type = 'pdf', plot_alpha = .01, lo_bound_hap = 0, \
		hi_bound_hap = 1050, lo_bound_dip = 0, hi_bound_dip = 1050):
	raw_data = pickle.load(open("raw_data.p","rb"))
	bin_data = pickle.load(open("binned_data.p","rb"))
	ploidy = ['Haploid', 'Diploid']
	plot_col = ['b', 'r']
	
	for n in range(2):
		sample = gene_name + ' ' + ploidy[n]
		temp_data = raw_data[sample]
		temp_bin = bin_data[sample]
		if ploidy[n] == 'Haploid':
			cut_ind_lo = NP.where(temp_bin[:,0] > lo_bound_hap)[0]
			cut_ind_hi = NP.where(temp_bin[:,0] < hi_bound_hap)[0]
		if ploidy[n] == 'Diploid':
			cut_ind_lo = NP.where(temp_bin[:,0] > lo_bound_dip)[0]
			cut_ind_hi = NP.where(temp_bin[:,0] < hi_bound_dip)[0]
		cut_ind = NP.intersect1d(cut_ind_lo, cut_ind_hi)
		temp_bin = temp_bin[cut_ind,:]
		PL.plot(temp_data[:,1], temp_data[:,3], '.' + plot_col[n], markersize = 12, \
			alpha = plot_alpha, label = "{0}".format(sample))
		PL.plot(temp_bin[:,0],temp_bin[:,3],'-k', linewidth = 2)
	
	PL.legend(loc = "upper left")
	PL.xlabel("SSC (Volume)")
	PL.ylabel("GFP Intensity (AU)")
	PL.ylim(-100,1050)
	PL.xlim(-100,1050)
	savepath = '{0}.{1}'.format(gene_name,file_type)
	PL.savefig(savepath)
	PL.close()

def plot_dbl_err(gene_name, file_type = 'pdf', plot_alpha = .01, lo_bound_hap = 0, \
		hi_bound_hap = 1050, lo_bound_dip = 0, hi_bound_dip = 1050):
	raw_data = pickle.load(open("raw_data.p","rb"))
	bin_data = pickle.load(open("binned_data.p","rb"))
	ploidy = ['Haploid', 'Diploid']
	plot_col = ['b', 'r']
	
	for n in range(2):
		sample = gene_name + ' ' + ploidy[n]
		temp_bin = bin_data[sample]
		if ploidy[n] == 'Haploid':
			cut_ind_lo = NP.where(temp_bin[:,0] > lo_bound_hap)[0]
			cut_ind_hi = NP.where(temp_bin[:,0] < hi_bound_hap)[0]
		if ploidy[n] == 'Diploid':
			cut_ind_lo = NP.where(temp_bin[:,0] > lo_bound_dip)[0]
			cut_ind_hi = NP.where(temp_bin[:,0] < hi_bound_dip)[0]
		cut_ind = NP.intersect1d(cut_ind_lo, cut_ind_hi)
		temp_bin = temp_bin[cut_ind,:]
		PL.plot(temp_bin[:,0],temp_bin[:,3],'-' + plot_col[n], linewidth = 2, \
			label = "{0}".format(sample))
		PL.plot(temp_bin[:,0],temp_bin[:,3] - temp_bin[:,2],'-' + plot_col[n], \
			linewidth = 1)
		PL.plot(temp_bin[:,0],temp_bin[:,3] + temp_bin[:,2],'-' + plot_col[n], \
			linewidth = 1)
	
	PL.legend(loc = "upper left")
	PL.xlabel("SSC (Volume)")
	PL.ylabel("GFP Intensity (AU)")
	PL.ylim(-100,1050)
	PL.xlim(-100,1050)
	savepath = '{0} bins.{1}'.format(gene_name,file_type)
	PL.savefig(savepath)
	PL.close()		
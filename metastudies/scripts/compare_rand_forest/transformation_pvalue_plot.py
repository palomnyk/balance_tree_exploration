#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for comparing random forest output to pvalues
# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os, sys
import time
from scipy import stats
import researchpy as rp
from statistics import mean
from matplotlib import markers
import math as math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from sklearn.metrics import accuracy_score, roc_auc_score
import argparse
import random

# --------------------------------------------------------------------------
print("Reading commmandline input with optparse.", flush = True)
# --------------------------------------------------------------------------

parser = argparse.ArgumentParser(description="This script runs a random forest test on various datasets.")
# parser.add_option("-f", "--file", dest="filename",
#                   help="write report to FILE", metavar="FILE")
parser.add_argument("-d", "--homedir",
                  default=os.path.expanduser(os.path.join("~", "git", "balance_tree_exploration")),
                  help="path to git balance treee exploration git repository", dest="homedir", metavar="homedir")
options, unknown = parser.parse_known_args()

# --------------------------------------------------------------------------
print("Establishing directory layout.", flush = True)
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(options.homedir)
projects = ["Vanderbilt", "Vangay", "Zeller", "Noguera-Julian"]
projects = ["Vanderbilt", "Noguera-Julian"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "pval_acc_vs_acc_python_by_transformation.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
comp_ds = ['alr_DADA2', 'clr_DADA2', 'DaDa2', 'Filtered_IQtree', \
	'Filtered_IQtree_blw.sqrt_enorm', 'Filtered_Silva_DADA2', \
	'Filtered_Silva_DADA2_blw.sqrt_enorm', 'Filtered_UPGMA_DADA2', \
	'Filtered_UPGMA_DADA2_blw.sqrt_enorm', 'lognorm_DADA2', 'Silva_DADA2', \
	'Silva_DADA2_blw.sqrt_enorm']

pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
#set font sizes
plt.rc('font', size=15) 
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('axes', labelsize=20) 
plt.rc('axes', titlesize=50)
for ds1 in comp_ds:
	print(ds1)
	pvalues = []
	ds1_name = []
	ds2_name = []
	ave_diff = [] #ds1_ave - ds2_ave
	transform = []
	proj = []
	for ds2 in comp_ds:
		if (ds1 != ds2):
			print(ds2)
			train_percent = 0.75
			print(f"{ds1} {ds2}")
			for project in projects:
				print(f"Adding project {project}")
				main_output_label = f"sklearn_random_forest_manual_{train_percent}train"
				op_dir = os.path.join(home_dir, project, "output")
				result_fpath = os.path.join(op_dir, "tables", f"{main_output_label}.csv")
				print(result_fpath)
				my_table = pd.read_csv(result_fpath, sep=',', header=0)
				#table 1
				ds1_table = my_table.loc[my_table["dataset"] == ds1,]
				splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
				# print(ds1_table[splits])
				ds1_means = ds1_table[splits].agg(mean, axis = 1)
				ds2_table = my_table.loc[my_table["dataset"] == ds2,]
				ds2_means = ds2_table[splits].agg(mean, axis = 1)
				my_tpval = stats.ttest_ind(ds1_means, ds2_means).pvalue
				pvalues.append(my_tpval)
				ds1_name.append(ds1)
				ds2_name.append(ds2)
				ave_diff.append(sum(ds1_means) - sum(ds2_means)) #ds1_ave - ds2_ave
				proj.append(project)
				transform.append(ds2)
					# stats.ttest_ind(df['bp_after'][df['sex'] == 'Male'],
          #       df['bp_after'][df['sex'] == 'Female'])
	# print(ave_diff)
	print(pvalues )
	# print(ds1_name)
	# print(ds2_name)
	# print(proj)
	print("build graphic")
	fig = plt.figure(figsize=(11,11))
	fig.suptitle(f"Metastudy {train_percent}training {ds1} vs others by accuracy, Python only")
	plt.subplots_adjust(bottom=0.8)
	ax = fig.add_subplot(1,1,1)
	my_markers = ["o", "s", "P", "v", "x", "1", "*", "+", "-", "D", "|"]
	my_transforms = list(set(transform))
	for i in range(len(proj)):
		my_proj = proj[i]
		# my_label
		my_marker = my_markers[projects.index(my_proj)]
		ax.scatter(ave_diff[i], math.log10(pvalues[i]), s=70, label=my_transforms[my_transforms.index(transform[i])])
	# plt.annotate(label, (x_lst[i], y_lst[i]))
	plt.axhline(y = math.log10(0.05), color = 'r', label="p=0.05")
	# ax.plot([0], [0,1], color = "r", label = "expected")
	ax.set_xlabel(f"mean difference in accuracy between {ds1} and others")
	ax.set_ylabel(f"log10 pvalue")
	ax.legend(title="Legend", loc="bottom right", framealpha=1)
	fig.tight_layout()
	print("Saving figure to pdf", flush = True)
	pdf.savefig( fig )
# print("Making seperate legend.")
# fig = plt.figure(figsize=(11,11))
# fig.suptitle(f"Legend {train_percent}training Python RF accuracy vs accuracy")
# plt.subplots_adjust(bottom=0.8)
# ax = fig.add_subplot(1,1,1)
# ax.plot([0,1], [0,1], color = "r", label = "expected")
# ax.plot(ds1_lst, a*ds1_lst+b, color = "green", label = "accuracy polyfit")
# fig.tight_layout()
# my_projects = list(set(list(map(lambda x: x[1], list(ds2_score.values())))))#pulling second element from each dict.value
# my_markers = ["o", "s", "P", "v", "x"]
# for i, label in enumerate(list(ds2_score.keys())):
# 	my_proj = ds2_score[label][1]
# 	print(f"{my_proj} {ds1_lst[i]} {ds2_lst[i]}, {list(ds2_score.keys())[i]}")
# 	my_marker = my_markers[my_projects.index(my_proj)]
# 	ax.scatter(0, 0, s=70, label=list(ds2_score.keys())[i], marker=my_marker)
# 	ax.legend(title="Legend",  loc="center", framealpha=1, mode = "expand", markerscale=2)
# print("Saving figure to pdf", flush = True)
# pdf.savefig( fig, bbox_inches='tight' )

print("Saving pdf", flush = True)
pdf.close()


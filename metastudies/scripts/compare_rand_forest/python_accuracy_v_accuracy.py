#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for comparing random forest output to pvalues
# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os, sys
import time
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

output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "accuracy_vs_accuracy_python_by_transformation.pdf")
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
font1 = {'family':'serif','color':'blue','size':20}
font2 = {'family':'serif','color':'darkred','size':15}
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
	for ds2 in comp_ds:
		if (ds1 != ds2):
			print(ds2)
			train_percent = 0.75
			ds1_score = {}
			ds2_score = {}
			features = []
			proj = []
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
				# print(ds1_table)
				splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
				# print(ds1_table[splits])
				means = ds1_table[splits].agg(mean, axis = 1)
				assert not means.empty, f"{ds1} is not in the table from {project}"
				for feat, ave in zip(list(ds1_table["metadata"].values) ,list(means)):
					ds1_score[f"{project}_{feat}"] = [ave, project]
				#table 2
				ds2_table = my_table.loc[my_table["dataset"] == ds2,]
				splits = ds2_table.columns[ds2_table.columns.str.startswith('split')].tolist()
				means = ds2_table[splits].agg(mean, axis = 1)
				assert not means.empty, f"{ds2} is not in the table from {project}"
				# print(my_table)
				for feat, ave in zip(list(ds2_table["metadata"].values) ,list(means)):
					ds2_score[f"{project}_{feat}"] = [ave, project]
			print("build graphic")
			# same_keys = set(ds2_score.keys()).intersection(set(ds1_score.keys()))
			# ds1_score = {key:ds1_score[key] for key in same_keys}
			# ds2_score = {key:ds2_score[key] for key in same_keys}
			ds1_lst = np.array(list(map(lambda x: x[0], list(ds1_score.values()))))
			ds2_lst = np.array(list(map(lambda x: x[0], list(ds2_score.values()))))
			a, b = np.polyfit(ds1_lst, ds2_lst, 1)
			# print(a, b)
			fig = plt.figure(figsize=(11,11))
			fig.suptitle(f"Metastudy {train_percent}training {ds1} vs {ds2} by accuracy, Python only")
			plt.subplots_adjust(bottom=0.8)
			ax = fig.add_subplot(1,1,1)
			my_projects = list(set(list(map(lambda x: x[1], list(ds2_score.values())))))#pulling second element from each dict.value
			my_markers = ["o", "s", "P", "v", "x"]
			for i, label in enumerate(list(ds2_score.keys())):
				if ds1_lst[i] > 0 and ds2_lst[i] > 0:
					my_proj = ds2_score[label][1]
					print(f"{my_proj} {ds1_lst[i]} {ds2_lst[i]}, {list(ds2_score.keys())[i]}")
					my_marker = my_markers[my_projects.index(my_proj)]
					ax.scatter(ds1_lst[i], ds2_lst[i], s=70, label=list(ds2_score.keys())[i], marker=my_marker)
			# plt.annotate(label, (x_lst[i], y_lst[i]))
			ax.plot([0,1], [0,1], color = "r", label = "expected")
			ax.plot(ds1_lst, a*ds1_lst+b, color = "green", label = "accuracy polyfit")
			ax.set_xlabel(f"Accuracy {ds1}")
			ax.set_ylabel(f"Accuracy {ds2}")
			# ax.legend(title="Legend", loc="upper left", framealpha=1)
			fig.tight_layout()
			print("Saving figure to pdf", flush = True)
			pdf.savefig( fig )
print("Making seperate legend.")
fig = plt.figure(figsize=(11,11))
fig.suptitle(f"Legend {train_percent}training Python RF accuracy vs accuracy")
plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
ax.plot([0,1], [0,1], color = "r", label = "expected")
ax.plot(ds1_lst, a*ds1_lst+b, color = "green", label = "accuracy polyfit")
fig.tight_layout()
my_projects = list(set(list(map(lambda x: x[1], list(ds2_score.values())))))#pulling second element from each dict.value
my_markers = ["o", "s", "P", "v", "x"]
for i, label in enumerate(list(ds2_score.keys())):
	my_proj = ds2_score[label][1]
	print(f"{my_proj} {ds1_lst[i]} {ds2_lst[i]}, {list(ds2_score.keys())[i]}")
	my_marker = my_markers[my_projects.index(my_proj)]
	ax.scatter(0, 0, s=70, label=list(ds2_score.keys())[i], marker=my_marker)
	ax.legend(title="Legend",  loc="center", framealpha=1, mode = "expand", markerscale=2)
print("Saving figure to pdf", flush = True)
pdf.savefig( fig, bbox_inches='tight' )

print("Saving pdf", flush = True)
pdf.close()


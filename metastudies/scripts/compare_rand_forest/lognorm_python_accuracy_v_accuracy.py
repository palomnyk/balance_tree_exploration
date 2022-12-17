#!/usr/bin/env python
# Author: Aaron Yerke, aaronyerke@gmail.com
# This is a script for comparing random forest output to pvalues

print(f"""Running {__file__}.
This is a script for comparing random forest output with pvalues.
""")
# --------------------------------------------------------------------------
print("Loading external libraries.",flush = True)
# --------------------------------------------------------------------------
import os
from statistics import mean
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from scipy.stats import linregress
import argparse

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
projects = ["Jones", "Vangay", "Zeller", "Noguera-Julian"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "lognorm_accuracy_vs_accuracy_python_by_transformation.pdf")
r_sq_boxplot_pdf_fpath = os.path.join(output_dir, "lognorm_rsq_accuracy_v_accuracy_py_by_transf.pdf")
train_percent = 0.75
# --------------------------------------------------------------------------
print("Establishing other constants.", flush = True)
# --------------------------------------------------------------------------
font1 = {'family':'serif','color':'blue','size':20}
font2 = {'family':'serif','color':'darkred','size':15}
# comp_ds = ['alr_DADA2', 'clr_DADA2', 'raw_DADA2', 'Filtered_IQtree_blw.sqrt_enorm', \
# 	"lognorm_Silva_DADA2", 'Filtered_Silva_DADA2', 'Filtered_Silva_DADA2_blw.sqrt_enorm', \
# 	'Filtered_UPGMA_DADA2_blw.sqrt_enorm', 'Silva_DADA2', 'Silva_DADA2_blw.sqrt_enorm']
comp_ds = ['alr_DADA2', 'clr_DADA2', 'raw_DADA2', 'Filtered_IQtree_blw.sqrt_enorm', \
	'Filtered_Silva_DADA2', 'Filtered_Silva_DADA2_blw.sqrt_enorm', \
	'Filtered_UPGMA_DADA2_blw.sqrt_enorm', 'Silva_DADA2', 'Silva_DADA2_blw.sqrt_enorm']


pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
#set font sizes
plt.rc('font', size=15) 
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 
plt.rc('axes', labelsize=20) 
plt.rc('axes', titlesize=25)

r_sq = []

fig = plt.figure(figsize=(15,30))
num_rows = 3
num_cols = 3
ax_count = 1
ds1 = "lognorm_DADA2"
for d2 in range(len(comp_ds)):
	ds2 = comp_ds[d2]
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
		splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
		means = ds1_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"{ds1} is not in the table from {project}"
		for feat, ave in zip(list(ds1_table["metadata"].values) ,list(means)):
			if ave >= 0:
				ds1_score[f"{project}_{feat}"] = [ave, project]
			else:
				ds1_score[f"{project}_{feat}"] = [0, project]
		#table 2
		ds2_table = my_table.loc[my_table["dataset"] == ds2,]
		splits = ds2_table.columns[ds2_table.columns.str.startswith('split')].tolist()
		means = ds2_table[splits].agg(mean, axis = 1)
		assert not means.empty, f"{ds2} is not in the table from {project}"
		# print(my_table)
		for feat, ave in zip(list(ds2_table["metadata"].values) ,list(means)):
			if ave >= 0:
				ds2_score[f"{project}_{feat}"] = [ave, project]
			else:
				ds2_score[f"{project}_{feat}"] = [0, project]
	print(f"preprocessing number of ds1 metadata cats: {len(ds1_score.keys())}")
	print(f"preprocessing number of ds2 metadata cats: {len(ds2_score.keys())}")
	print(f"Throwing out these meta cats:")
	my_dif = list( set(ds1_score.keys()) - set(ds2_score.keys()))
	print(my_dif)
	print("build graphic")
	same_keys = set(ds2_score.keys()).intersection(set(ds1_score.keys()))
	ds1_score = {key:ds1_score[key] for key in same_keys}
	ds2_score = {key:ds2_score[key] for key in same_keys}
	# ds1_score = sorted(ds1_score)
	# ds2_score = sorted(ds2_score)
	ds1_lst = np.array(list(map(lambda x: x[0], list(ds1_score.values()))))
	ds2_lst = np.array(list(map(lambda x: x[0], list(ds2_score.values()))))
	slope, intercept, r_value, p_value, std_err = linregress(ds1_lst, ds2_lst)
	r_sq.append(r_value**2)
	a, b = np.polyfit(ds1_lst, ds2_lst, 1)
	# print(a, b)
	# fig.suptitle(f"Metastudy {train_percent}training {ds1} vs {ds2} by accuracy, Python only")
	plt.subplots_adjust(bottom=0.8)
	ax = fig.add_subplot(num_rows,num_cols, ax_count)
	my_projects = list(set(list(map(lambda x: x[1], list(ds2_score.values())))))#pulling second element from each dict.value
	my_markers = ["o", "s", "P", "v", "x"]
	assert( list(ds2_score.keys()) == set(list(ds2_score.keys())), "Keys are not a set.")
	print(f"number of metata cats {len(ds2_score.keys())}")
	for i, labl in enumerate(list(sorted(ds2_score.keys()))):
		my_proj = ds2_score[labl][1]
		# print(f"{my_proj} {ds1_lst[i]} {ds2_lst[i]}, {list(sorted(ds2_score.keys()))[i]}")
		my_marker = my_markers[my_projects.index(my_proj)]
		ax.scatter(ds1_lst[i], ds2_lst[i], s=100, label=list(sorted(ds2_score.keys()))[i], marker=my_marker)
	# plt.annotate(label, (x_lst[i], y_lst[i]))
	ax.plot([0,1], [0,1], color = "r", label = "expected")
	ax.plot(ds1_lst, a*ds1_lst+b, color = "green", label = "accuracy polyfit")
	ax.text(0.1, 0.8, f"r squared: {round(r_value**2, 4)}", fontsize=20)
	ax.set_xlabel(f"{ds1}")
	ax.set_ylabel(f"{ds2}")
	ax.set_title("Scores")
	ax_count += 1
	# if ax_count > 10:
	# 	fig.tight_layout(pad = 1, h_pad=1, w_pad=2)
	# 	plt.subplots_adjust(left=0.1)
	# 	print("Saving figure to pdf", flush = True)
	# 	pdf.savefig( fig )
	# 	print("Making seperate legend.")
	# 	fig = plt.figure(figsize=(11,26))
	# 	ax = fig.add_subplot(1,1,1)
	# 	ax.legend(title="",  loc="center", framealpha=1, mode = "expand", markerscale=2)
	# 	print("Saving figure to pdf", flush = True)
	# 	pdf.savefig( fig )
	# ax.legend(title="Legend", loc="upper left", framealpha=1)
fig.tight_layout(pad = 1, h_pad=1, w_pad=2)
plt.subplots_adjust(left=0.1)
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )
print("Making seperate legend.")
fig = plt.figure(figsize=(11,26))
ax = fig.add_subplot(1,1,1)
# plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
# fig.tight_layout()
my_markers = ["o", "s", "P", "v", "x"]
for i, label in enumerate(list(sorted(ds2_score.keys()))):
	my_proj = ds2_score[label][1]
	print(f"{my_proj} {ds1_lst[i]} {ds2_lst[i]}, {list(sorted(ds2_score.keys()))[i]}")
	my_marker = my_markers[my_projects.index(my_proj)]
	ax.scatter(ds1_lst[i], ds2_lst[i], s=100, label=list(sorted(ds2_score.keys()))[i], marker=my_marker)
	ax.legend(title="",  loc="center", framealpha=1, mode = "expand", markerscale=2)
print("Saving figure to pdf", flush = True)
pdf.savefig( fig )

print("Saving pdf", flush = True)
pdf.close()

print(f"Accuracy vs accuracy R squared plots")
pdf = matplotlib.backends.backend_pdf.PdfPages(r_sq_boxplot_pdf_fpath)
fig = plt.figure(figsize=(11,11))
fig.suptitle(f"Accuracy vs accuracy R squared")
plt.subplots_adjust(bottom=0.8)
ax = fig.add_subplot(1,1,1)
bp = ax.boxplot(r_sq)
fig.tight_layout()
pdf.savefig( fig, bbox_inches='tight')

print(f"Saving pdf to {plot_pdf_fpath}", flush = True)
pdf.close()

print(f"{__file__} complete!")

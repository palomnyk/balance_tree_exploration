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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import matplotlib.colors as mcolors
from pandas.api.types import is_string_dtype
from requests import head
from sklearn import model_selection
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
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
# projects = ["Vanderbilt", "Vangay", "Zeller", "Noguera-Julian"]
projects = ["Vanderbilt", "Vangay"]
output_dir = os.path.join(home_dir, "metastudies", "output")
assert os.path.exists(output_dir)
plot_pdf_fpath = os.path.join(output_dir, "python_by_transformation.pdf")

comp_ds = ['alr_DADA2', 'clr_DADA2', 'DaDa2', 'Filtered_IQtree', \
	'Filtered_IQtree_blw.sqrt_anorm', 'Filtered_Silva_DADA2', \
	'Filtered_Silva_DADA2_blw.sqrt_anorm', 'Filtered_UPGMA_DADA2', \
	'Filtered_UPGMA_DADA2_blw.sqrt_anorm', 'lognorm_DADA2', 'Silva_DADA2', \
	'Silva_DADA2_blw.sqrt_anorm']

pdf = matplotlib.backends.backend_pdf.PdfPages(plot_pdf_fpath)
for ds1 in comp_ds:
	print(ds1)
	for ds2 in comp_ds:
		if (ds1 != ds2):
			print(ds2)
			train_percent = 0.75
			pvals = {}
			ds1_score = {}
			ds2_score = {}
			features = []
			proj = []
			print(f"{ds1} {ds2}")
			for project in projects:
				print(f"Adding project {project}")
				main_output_label = f"sklearn_random_forest_manual_{train_percent}train_{project}_data"
				op_dir = os.path.join(home_dir, project, "output")
				result_fpath = os.path.join(op_dir, "tables", f"{main_output_label}.csv")
				print(result_fpath)
				my_table = pd.read_csv(result_fpath, sep=',', header=0)
				#table 1
				ds1_table = my_table.loc[my_table["dataset"] == ds1,]
				splits = ds1_table.columns[ds1_table.columns.str.startswith('split')].tolist()
				means = ds1_table[splits].agg(mean, axis = 1)
				assert not means.empty(), f"{ds1} is not in the table from {project}"
				for feat, ave in zip(list(ds1_table["metadata"].values) ,list(means)):
					ds1_score[feat] = ave
				#table 2
				ds2_table = my_table.loc[my_table["dataset"] == ds2,]
				splits = ds2_table.columns[ds2_table.columns.str.startswith('split')].tolist()
				means = ds2_table[splits].agg(mean, axis = 1)
				assert not means.empty(), f"{ds2} is not in the table from {project}"
				pval_fpath = os.path.join(op_dir, "tables", f"{project}_pValuesUnivariate_taxaVmetadata.csv")
				my_table = pd.read_csv(pval_fpath, sep=',', header=0)
				# print(my_table)
				for feat, ave in zip(list(ds2_table["metadata"].values) ,list(means)):
					ds2_score[feat] = ave
					pvals[feat] = min(my_table.loc[(my_table["meta_name"]==feat) & (my_table["taxa_lev"] == "Genus"), "pval"].values)
			#build graphic
			fig = plt.figure(figsize=(11,11))
			fig.suptitle(f"Metastudy {train_percent}training {ds1} vs {ds2} by pvalue")
			plt.subplots_adjust(bottom=0.8)
			ax = fig.add_subplot(1,1,1)
			# print("pvals")
			# print(list(pvals.values()))
			# print("ds1")
			# print(list(ds1_score.values()))
			print("pvals")
			print(pvals)
			print("ds1")
			print(ds1_score)
			print("ds2")
			print(ds2_score)
			print("proj")
			print(proj)
			ax.scatter(pvals.values(), ds1_score.values())
			# colors = list([sublist[-1] for sublist in tables])
			# for patch, color in zip(bp['boxes'], colors):
			# 	patch.set_facecolor(color)
			# ax.axhline(np.nanmean(plot_data), c="r", linestyle="dashed")
			# ax.set_xticklabels(labels = plot_data.columns, rotation=90)
			# ax.tick_params(axis='x', which='major', labelsize=10)
			#for boxplot
			fig.tight_layout()
			print("Saving pdf", flush = True)
			pdf.savefig( fig )
			# sys.exit()
			# print(my_table.columns)
			# print(sorted(set(my_table["dataset"].values),key=lambda s: s.lower()))

print("Saving pdf", flush = True)
pdf.close()


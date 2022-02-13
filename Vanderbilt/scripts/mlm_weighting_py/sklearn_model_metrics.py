#!/usr/bin/env python3

"""
Author: Aaron Yerke (aaronyerke@gmail.com)
For determining if PhILR weighting schemes improve any ML techniques.
This was useful: https://machinelearningmastery.com/compare-machine-learning-algorithms-python-scikit-learn/
Compare Algorithms
"""

# --------------------------------------------------------------------------
print("Loading external libraries.")
# --------------------------------------------------------------------------
from cgi import print_directory
from importlib.resources import path
import os, sys
from posixpath import split
from tkinter.ttk import Style
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

print("Establishing directory layout.")
home_dir = os.path.expanduser(os.path.join("~", "git", "balance_tree_exploration"))#expanduser allows tilda (~) to work
project = "Vanderbilt"
output_dir = os.path.join(home_dir, project, "output")
silva_philr_dir = os.path.join(output_dir, "tables", "silva_philr_weights")

print("Establishing other constants.")
main_output_label = "sklearn_ml_test"
philr_part_weights = ["uniform","gm.counts","anorm","anorm.x.gm.counts","enorm","enorm.x.gm.counts"]
philr_ilr_weights = ["uniform","blw","blw.sqrt","mean.descendants"]

print("Importing data to working env.")
meta_df = pd.read_csv(os.path.join(home_dir, project, "patient_metadata.tsv"), sep='\t', header=0, index_col=0)

# prepare configuration for cross validation test harness
seed = 7
# prepare models
models = []
models.append(('LR', LogisticRegression()))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('KNN', KNeighborsClassifier()))
models.append(('DTREE', DecisionTreeClassifier()))
models.append(('RF', RandomForestClassifier()))
models.append(('NB', GaussianNB()))
models.append(('SVM', SVC()))

col_names = ["ilr_weight", "part_weight", "model", "split1", "split2", "split3", "split4", "split5", "split6", "split7", "split8", "split9", "split10"]
result_fpath = os.path.join(output_dir, "tables", "sklearn_model_metrics.csv")
# os.chdir(os.path.join(output_dir, "tables"))
with open(result_fpath, "w+") as f:
	f.write(",".join(col_names))
	f.write("\n")
	for pw in philr_part_weights:
		for iw in philr_ilr_weights:
			table_fn = f"ref_tree_cln_{iw}_{pw}.csv"
			my_df = pd.read_csv(os.path.join(silva_philr_dir, table_fn), sep=',', header=0, index_col=0)
			print("Dropping any extra values from metadata. Did it work?")
			meta_df = meta_df.loc[list(my_df.index.values)]
			print(list(my_df.index.values) == list(meta_df.index.values))
			spetz_var = meta_df["type"]
			print("evaluate each model in turn.")
			scoring = 'accuracy'
			for name, model in models:
				kfold = model_selection.KFold(n_splits=10, random_state=seed, shuffle=True)
				cv_results = model_selection.cross_val_score(model, my_df, spetz_var, cv=kfold, scoring=scoring)
				# result_str = np.array2string(cv_results, separator=",",suffix="/n")
				result_str = ",".join(map(str, cv_results.tolist()))
				msg = f"{iw},{pw},{name},{result_str}\n"
				f.write(msg)
			# boxplot algorithm comparison
print("Finished recording accuracy.")

result_df = pd.read_csv(result_fpath, sep=',', header=0)
print(result_df.head())
# print(" ".join(result_df.index))

# pdf = matplotlib.backends.backend_pdf.PdfPages(os.path.join(output_dir, "graphics", main_output_label))

algos = list(set(result_df.loc[:,"model"]))
print(algos)
num_cols = 2
num_rows = abs(-len(algos)//num_cols)
print(result_df.head())
overall_mean = np.nanmean(result_df.iloc[:,3:])
fig = plt.figure(figsize=(11,11))

# plt.subplots_adjust(bottom=0.8)
for al in range(len(algos)):
	algo = algos[al]
	print(algo)
	fig_df = pd.DataFrame(result_df[result_df["model"] == algo])
	plot_data = fig_df.iloc[:,3:].transpose()
	print(f"df shape: {fig_df.shape[0]} {fig_df.shape[1]}")
	# boxplot algorithm comparison
	fig.suptitle(f"{project} PhILR weighting machine learning accuracy")
	ax = fig.add_subplot(num_rows,num_cols, al +1)
	ax.boxplot(plot_data)
	ax.title.set_text(f"{algo} by weighting scheme")
	ax.axhline(np.nanmean(plot_data), c="r", linestyle="dashed")
	ax.axhline(overall_mean, c="g", linestyle = ("-."))
	ax.locator_params(axis='y', tight=True, nbins=4)
	# print(fig_df.iloc[:,3:].head())
	# print(fig_df.loc[:,"ilr_weight"].values)
	new_labs = [f"{x}\n{y}" for x,y in zip(fig_df.loc[:,"ilr_weight"].values, fig_df.loc[:,"part_weight"].values)]
	# ax.set_xticklabels(fig_df.loc[:,"ilr_weight"].tolist(), rotation=90)
	ax.set_xticklabels(new_labs, rotation=90)
	ax.tick_params(axis='x', which='major', labelsize=6)

fig.tight_layout()

fig.savefig(os.path.join(output_dir, "graphics", f"bp_{main_output_label}.png"))
# pdf.close()

print("Python script completed.")



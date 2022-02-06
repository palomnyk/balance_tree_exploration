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
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.linear_model import LogisticRegression
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

print("Establishing directory layout.")
home_dir = os.path.join("~", "git", "balance_tree_exploration")
project = "Vanderbilt"
output_dir = os.path.join(home_dir, project, 'output')

print("Establishing other constants.")
main_output_label = "testing_pycaret"

print("Importing data to working env.")
my_df = pd.read_csv(os.path.join(output_dir, "tables", "philr_denovo_tree_UPGMA_1000_sd_filtered.csv"),\
   sep=',', header=0, index_col=0)
meta_df = pd.read_csv(os.path.join(home_dir, project, "patient_metadata.tsv"),\
   sep='\t', header=0, index_col=0)

print("Dropping any extra values from metadata. Did it work?")
meta_df = meta_df.loc[list(my_df.index.values)]
list(my_df.index.values) == list(meta_df.index.values)

spetz_var = meta_df["type"]

# prepare configuration for cross validation test harness
seed = 7
# prepare models
models = []
models.append(('LR', LogisticRegression()))
models.append(('LDA', LinearDiscriminantAnalysis()))
models.append(('KNN', KNeighborsClassifier()))
models.append(('CART', DecisionTreeClassifier()))
models.append(('NB', GaussianNB()))
models.append(('SVM', SVC()))
# evaluate each model in turn
results = []
names = []
scoring = 'accuracy'
for name, model in models:
	kfold = model_selection.KFold(n_splits=10, random_state=seed)
	cv_results = model_selection.cross_val_score(model, my_df, spetz_var, cv=kfold, scoring=scoring)
	results.append(cv_results)
	names.append(name)
	msg = f"{name}: {cv_results.mean()} {cv_results.std()}"
	print(msg)
# boxplot algorithm comparison
fig = plt.figure()
fig.suptitle('Algorithm Comparison')
ax = fig.add_subplot(111)
plt.boxplot(results)
ax.set_xticklabels(names)
plt.savefig(fname=os.path.join(output_dir, "graphics", f"ml_sklearn_{scoring}.png"),\
   format = "png")

print("Python script completed.")

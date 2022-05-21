
# --------------------------------------------------------------------------
print("Loading external libraries.")
# --------------------------------------------------------------------------
import os, sys
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
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.naive_bayes import GaussianNB
from sklearn.metrics import accuracy_score
import argparse
import random

# --------------------------------------------------------------------------
print("Establishing directory layout.")
# --------------------------------------------------------------------------
home_dir = os.path.expanduser(os.path.join("~", "git", "balance_tree_exploration"))
project = "Vanderbilt"
output_dir = os.path.join(home_dir, project, "output")

# --------------------------------------------------------------------------
print("Importing data to working env.")
# --------------------------------------------------------------------------
meta_df = pd.read_csv(os.path.join(home_dir, project, "patient_metadata.tsv"), \
	sep="\t", header=0, index_col="Run")
hashseq_df = pd.read_csv(os.path.join(output_dir, "hashseq", "hashseq.csv"), sep=",", header=0, index_col=0)
asv_table = pd.read_csv(os.path.join(output_dir, "tables", "ForwardReads_DADA2.txt"), sep="\t", header=0, index_col=0)
clr_table = pd.read_csv(os.path.join(output_dir, "tables", "clr_asv.csv"), sep=",", header=0, index_col=0)
alr_table = pd.read_csv(os.path.join(output_dir, "tables", "alr_asv.csv"), sep=",", header=0, index_col=0)
ln_table = pd.read_csv(os.path.join(output_dir, "tables", "lognorm_asv.csv"), sep=",", header=0, index_col=0)
ln_hs_tab = pd.read_csv(os.path.join(output_dir,"tables", "lognorm_hashseq.csv"))
HashSeq_clr = pd.read_csv(os.path.join(output_dir,"tables", "clr_hashseq.csv"))
HashSeq_alr = pd.read_csv(os.path.join(output_dir,"tables", "alr_hashseq.csv"))

# meta_df = meta_df.loc[list(asv_table.index.values)]#drops rows from metadata that aren't in asv_table
# if all(meta_df.index == hashseq_df.index):
# 	print("dataframes are the same.")

# --------------------------------------------------------------------------
print("Establishing other constants")
# --------------------------------------------------------------------------
metad_cols = range(len(meta_df.columns))
seed = 7
scoring = "accuracy"
main_output_label = "sklearn_random_forest_manual"
result_fpath = os.path.join(output_dir, "tables", f"{main_output_label}_{project}_data.csv")
col_names = ["dataset", "metadata", "split1", "split2", "split3", "split4", "split5", "split6", "split7", "split8", "split9", "split10"]
pdf_fpath = os.path.join(output_dir, "graphics", f"bp_{main_output_label}_{project}.pdf")
num_rf_iterations = 10

# ----------------------------------------------------------------------------
print("Setting up tables to feed the random forest model.")
# # --------------------------------------------------------------------------
tables = []
tables.append(("DaDa2", asv_table))
tables.append(("HashSeq", hashseq_df))
tables.append(("lognorm_DADA2", ln_table))
tables.append(("lognorm_HashSeq", ln_hs_tab))
tables.append(("alr_DADA2", alr_table))
tables.append(("alr_HashSeq", HashSeq_alr))
tables.append(("clr_DADA2", clr_table))
tables.append(("clr_HashSeq", HashSeq_clr))

print(len(tables))
philr_part_weights = ["anorm","enorm"]
philr_ilr_weights = ["blw.sqrt","mean.descendants"]
silva_philr_dir = os.path.join(output_dir, "tables", "silva_philr_weights")
if not os.path.exists(silva_philr_dir):
  print(f"{silva_philr_dir} does not exist. Use silva_philr_weights.R to create it.")
  sys.exit()
for pw in philr_part_weights:
	for iw in philr_ilr_weights:
		table_fn = f"ref_tree_cln_{iw}_{pw}.csv"
		my_df = pd.read_csv(os.path.join(silva_philr_dir, table_fn), sep=',', header=0, index_col=0)
		my_label = f"ref_tree_philr_{iw}_{pw}"
		tables.append((my_label, my_df))

# meta_df.head
# meta_df = meta_df.sample(frac=1)
# meta_df.head

# --------------------------------------------------------------------------
print(f"Running random forest model to find {scoring}.")
# # --------------------------------------------------------------------------
with open(result_fpath, "w+") as fl:
	fl.write(",".join(col_names))
	fl.write("\n")
	for meta_c in metad_cols:
		m_c = list(meta_df.columns)[meta_c]
		# meta_df = meta_df.loc[list(my_table.index.values)
		for name, my_table in tables:
			my_accuracy = [0] * num_rf_iterations
			random.seed(10)
			for i in range(num_rf_iterations):
				rand_int = random.randint(0, 1000)
				spetz_var = meta_df.loc[list(my_table.index.values),m_c]#metadata var to test
				pred_train, pred_test, resp_train, resp_test = model_selection.train_test_split(my_table, spetz_var, test_size=0.75, random_state=rand_int, shuffle=True) 
				if is_string_dtype(spetz_var) == True and spetz_var.isnull().sum() < 5:
					clf = RandomForestClassifier(max_depth=2, random_state=0)
					clf.fit(pred_train, resp_train)
					resp_pred = clf.predict(pred_test)
					my_accuracy[i] = accuracy_score(resp_test, resp_pred)
			final_acc = ",".join(map(str, my_accuracy))
			print(final_acc)
			msg = f"{name},{m_c},{final_acc}\n"
			print(msg)
			fl.write(msg)
		print(f"{name}, {m_c} mean is: {mean(my_accuracy)}")
print("Finished recording accuracy.")

# --------------------------------------------------------------------------
print(f"Building boxplot PDF.")
# --------------------------------------------------------------------------
#Setup for building boxplots
result_df = pd.read_csv(result_fpath, sep=',', header=0)
metadata_cats = list(set(result_df["metadata"]))
num_cols = 2
num_rows = abs(-len(tables)//num_cols)

pdf = matplotlib.backends.backend_pdf.PdfPages(pdf_fpath)
for meta_c in metadata_cats:
	fig = plt.figure(figsize=(11,11))
	fig.suptitle(f"{project} random forest manual split 4-fold {scoring} {meta_c}")
	plt.subplots_adjust(bottom=0.8)
	meta_result_df = pd.DataFrame(result_df[result_df["metadata"] == meta_c])
	# flat_num_only = pd.DataFrame(meta_result_df.iloc[:,5:]).to_numpy().flatten()
	plot_data = meta_result_df.iloc[:,2:].transpose()
	f_mean = np.nanmean(plot_data)
	ax = fig.add_subplot(1,1,1)
	ax.boxplot(plot_data)
	ax.axhline(np.nanmean(plot_data), c="r", linestyle="dashed")
	ax.axhline(f_mean, c="g", linestyle = ("-."))
	ax.set_xticklabels(meta_result_df["dataset"].tolist(), rotation=90)
	ax.tick_params(axis='x', which='major', labelsize=15)
	#for boxplot
	fig.tight_layout()
	pdf.savefig( fig )

print("Saving pdf")
pdf.close()

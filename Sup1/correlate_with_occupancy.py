import numpy as np
import os
from matplotlib import pyplot as plt
from scipy import stats as st
import sys

cdt_path = sys.argv[1]
tsv_path = sys.argv[2]
name = sys.argv[3]

values = {}

with open(tsv_path) as tsv:
	for line in tsv:
		line = line.strip().split()
		values[line[0]] = float(line[1])
	tsv.close()

for file_name in os.listdir(cdt_path):
	occupancy = {}
	with open(os.path.join(cdt_path, file_name)) as cdt:
		for line in cdt:
			line = line.strip().split()
			if line[0] != "Uniqe":	#skip header
				gene = line[0]
				if "%" in gene:	#remove trailing characters
					gene = gene.split("%")[0]
				occupancy[gene] = sum([float(line[i]) for i in range(1,len(line))])
		cdt.close()

	xs = []
	ys = []
	for gene in occupancy.keys():
		if gene in values.keys():
			if occupancy[gene] != 0 and values[gene] != 0 and np.isfinite(values[gene]):
				xs.append(np.log2(occupancy[gene]))
				ys.append(np.log2(values[gene]))

	r, p = st.pearsonr(xs, ys)

	plt.scatter(xs, ys)
	plt.xlabel("Occupancy (tag count, log2)")
	plt.ylabel("{}\n(log2)".format(name))
	plt.title("n={}\nR={}".format(len(xs), r))
	plt.savefig("{}_{}".format(name, file_name[0:len(file_name)-5]))
	plt.close()

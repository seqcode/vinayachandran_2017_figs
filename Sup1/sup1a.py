import sys
import os
import numpy as np
from sklearn.decomposition import PCA
from matplotlib import pyplot as plt

def calculate_sums(infile, sums_by_gene, index, num_files):
	with open(infile) as cdt:
		for i, line in enumerate(cdt):
			if i > 0:	#skip header
				line = line.strip().split()
				gene = line[0]
				if gene not in sums_by_gene.keys():
					sums_by_gene[gene] = np.zeros(num_files)	#initialize empty row
				sums_by_gene[gene][index] = sum(map(float, line[1:len(line)]))
		cdt.close()


in_dir = sys.argv[1]		
file_names = []
for file_name in os.listdir(in_dir):
	if file_name.endswith(".cdt"):
		file_names.append(file_name)

n = len(file_names)	#number of experiments
sums_by_gene = {}
for i, file_name in enumerate(file_names):
	calculate_sums(in_dir + "/" + file_name, sums_by_gene, i, n)
	print "finished processing {}".format(file_name)

genes = sums_by_gene.keys()
m = len(genes)
all_sums = np.zeros((n,m))

for i, gene in enumerate(genes):
	sums = sums_by_gene[gene]	
	for j, count_sum in enumerate(sums):
		all_sums[j,i] = count_sum

pca = PCA(n_components=2)
principal_components = pca.fit_transform(all_sums)
plt.scatter(principal_components[:,0], principal_components[:,1])
plt.show()

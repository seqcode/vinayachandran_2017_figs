import numpy as np
import os
from matplotlib import pyplot as plt
from scipy import stats as st

transcription_rates = {}

with open("../shared_files/holstege.tsv") as holstege:
	for line in holstege:
		line = line.strip().split()
		if line[1] != "NAN":
			transcription_rates[line[0]] = float(line[1])
	holstege.close()

for file_name in os.listdir("b_CDT"):
	occupancy = {}
	with open("b_CDT/" + file_name) as cdt:
		for line in cdt:
			line = line.strip().split()
			if line[0] != "Uniqe":	#skip header
				occupancy[line[0]] = sum([float(line[i]) for i in range(1,len(line))])
		cdt.close()

	xs = []
	ys = []
	for gene in occupancy.keys():
		if gene in transcription_rates.keys():
			if occupancy[gene] != 0 and transcription_rates[gene] != 0:
				xs.append(np.log2(occupancy[gene]))
				ys.append(np.log2(transcription_rates[gene]))

	r, p = st.pearsonr(xs, ys)

	plt.scatter(xs, ys)
	plt.xlabel("Occupancy (tag count, log2)")
	plt.ylabel("Transcription rate (mRNA/hr)\n(log2)")
	plt.title("n={}\nR={}".format(len(xs), r))
	plt.savefig(file_name[0:len(file_name)-4])
	plt.show()			

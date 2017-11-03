import numpy as np
from matplotlib import pyplot as plt

def get_tags(tab_file):
	#initialize with zeroes
	tags = {}
	with open("../shared_files/sacCer3.chrom.sizes") as chrom_sizes:
		for line in chrom_sizes:
			data = line.strip().split()
			chrom = data[0]
			for i in range(1, int(data[1])):
				key = "{}:{}".format(chrom, i)
				tags[key] = 0
		chrom_sizes.close()

	#fill with tags at each genomic coordinate 
	with open(tab_file) as tab:
		for line in tab:
			data = line.strip().split()
			chrom = data[0]
			index = data[1]
			key = "{}:{}".format(chrom, index)
			tags[key] = float(data[4]) 
		tab.close()

	return tags

window_size = 50

mhs_tags = get_tags("tab_files_c/Normalized_tab_files/50428sacCer3_.tab")
hs_tags = get_tags("tab_files_c/Normalized_tab_files/50429sacCer3_.tab")

counts = []

for orientation in ("upstream", "downstream"):
	#get genomic coordinates
	chroms = []
	coords = []
	with open("../shared_files/Hsf1-union-Xu-TSS-divergent-{}.gff".format(orientation)) as gff:
		for line in gff:
			data = line.split()
			chroms.append(data[0])
			coords.append(int(data[3]))
		gff.close()

	#tag counts at each TSS
	mhs_tss_counts = np.zeros_like(chroms, dtype=float)	
	hs_tss_counts = np.zeros_like(chroms, dtype=float)

	for i, (chrom, coord) in enumerate(zip(chroms, coords)):
		start = coord - window_size
		end = coord + window_size
	
		#sum tags in window
		mhs_tss_counts[i] = sum([mhs_tags["{}:{}".format(chrom, j)] for j in range(start,end)])	 
		hs_tss_counts[i] = sum([hs_tags["{}:{}".format(chrom, j)] for j in range(start,end)])
	
	#fold change
	log_ratios = np.zeros_like(mhs_tss_counts)
	for i, (mhs_tss_count, hs_tss_count) in enumerate(zip(mhs_tss_counts, hs_tss_counts)):
		if hs_tss_count != 0 and mhs_tss_count != 0:
			log_ratios[i] = np.log2(hs_tss_count/mhs_tss_count)

	upregulated = 0
	downregulated = 0
	for log_ratio in log_ratios:
		if log_ratio > 0:
			upregulated += 1
		elif log_ratio < 0:
			downregulated += 1

	counts.append(downregulated)
	counts.append(upregulated)

plt.bar(range(len(counts)), counts)
plt.xticks(range(len(counts)), ("Down", "Up", "Down", "Up"))
plt.xlabel("HS-induced changes in TFIIB occupancy")
plt.ylabel("# of genes")
plt.savefig("Fig4C.png")
plt.show()

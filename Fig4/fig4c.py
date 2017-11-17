from matplotlib import pyplot as plt
import numpy as np

counts = []

for orientation in ("upstream", "downstream"):
	#tag counts at each TSS
	mhs_tss_counts = np.loadtxt("tab_files_c/Normalized_tab_files/50428sacCer3_{}.txt".format(orientation))	
	hs_tss_counts = np.loadtxt("tab_files_c/Normalized_tab_files/50429sacCer3_{}.txt".format(orientation))

	#fold change
	upregulated = 0
	downregulated = 0
	for mhs_tss_count, hs_tss_count in zip(mhs_tss_counts, hs_tss_counts):
		if hs_tss_count != 0 and mhs_tss_count != 0:
			ratio = hs_tss_count/mhs_tss_count
			if ratio > 1:
				upregulated += 1
			elif ratio < 1:
				downregulated += 1

	counts.append(downregulated)
	counts.append(upregulated)

plt.bar(range(len(counts)), counts)
plt.xticks(range(len(counts)), ("Down", "Up", "Down", "Up"))
plt.xlabel("HS-induced changes in TFIIB occupancy")
plt.ylabel("# of genes")
plt.savefig("Fig4C.png")

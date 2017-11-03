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
	mhs_tss_counts = []	
	hs_tss_counts = []

	for chrom, coord in zip(chroms, coords):
		start = coord - window_size
		end = coord + window_size
	
		#sum tags in window
		mhs_tss_counts.append(sum([mhs_tags["{}:{}".format(chrom, j)] for j in range(start,end)]))	 
		hs_tss_counts.append(sum([hs_tags["{}:{}".format(chrom, j)] for j in range(start,end)]))
	
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
plt.show()

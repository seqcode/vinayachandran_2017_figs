import numpy as np

all_motifs = {}

with open("closest_tss_within_500bp.gff") as gff:
	for line in gff:
		line = line.split()
		geneid = line[17]
		if geneid in all_motifs.keys():
			all_motifs[geneid].append((line[0], line[3], line[4], line[5], line[6]))	#add motif information to list of motifs for this gene
		else:
			all_motifs[geneid] = [(line[0], line[3], line[4], line[5], line[6])]
	gff.close()

with open("assignment_list.tsv", "w") as out:
	for geneid in all_motifs.keys():
		motifs = np.array(all_motifs[geneid])
		occupancies = [motif[3] for motif in motifs]	#occupancy for each motif at this gene
		sort = np.argsort(occupancies)	
		sort = np.flip(sort, 0)	#sort by decreasing occupancy
		sorted_motifs = motifs[sort]
		motif1 = sorted_motifs[0]
		if len(sorted_motifs) > 1:
			motif2 = sorted_motifs[1]
			out.write("\t".join((geneid, motif1[0], motif1[1], motif1[2], motif1[4], motif2[0], motif2[1], motif2[2], motif2[4])))
		else:
			out.write("\t".join((geneid, motif1[0], motif1[1], motif1[2], motif1[4], "None", "None", "None", "None")))
		out.write("\n")
	out.close() 

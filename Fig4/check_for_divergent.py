import numpy as np

with open("upstream_opposite_strand_within_500bp.gff") as in_file:
	with open("divergent_list.tsv", "w") as out_file:
		for line in in_file:
			line = line.strip().split()
			geneid1 = line[8]
			geneid2 = line[17]
			if line[10] == "Xu_2009_ORFs":	#ignore non-coding features
				strand1 = geneid1[6]
				strand2 = geneid2[6]
				if strand1 != strand2:	#must be on opposite strands
					out_file.write(geneid1 + "\t" + geneid2 + "\n")
				else:
					out_file.write(geneid1 + "\tNone\n")
			else:
				out_file.write(geneid1 + "\tNone\n")
		out_file.close()
	in_file.close()

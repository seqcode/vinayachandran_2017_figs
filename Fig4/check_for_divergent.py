import numpy as np

with open("opposite_strand.gff") as in_file:
	with open("divergent_list.tsv", "w") as out_file:
		for line in in_file:
			line = line.strip().split()
			geneid1 = line[8]
			geneid2 = line[17]
			num1 = int(geneid1[3:6])
			strand1 = geneid1[6]
			num2 = int(geneid2[3:6])
			strand2 = geneid2[6]
			if strand1 != strand2 and np.abs(num1 - num2) == 1:
				out_file.write(geneid1 + "\t" + geneid2 + "\n")
			else:
				out_file.write(geneid1 + "\tNone\n")
		out_file.close()
	in_file.close()

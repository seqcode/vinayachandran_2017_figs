import numpy as np
import sys

def fill_dict(path):
	dct = {}
	with open(path) as in_file:
		for line in in_file:
			line = line.strip().split()
			if len(line) > 1:
				dct[line[0]] = line[1]
			else:
				dct[line[0]] = None
		in_file.close()
	return dct

hsf1_flank = 100
spt3_flank = 150
tfiib_flank = 150

header = ("chr", "status", "TSS1 geneID", "TSS1 class", "TSS1 coordinate", "TSS1 strand", "TSS2 geneID", "TSS2 class", "TSS2 coordinate", "TSS2 strand", "MHS TFIIB TSS1 occupancy", "HS3 TFIIB TSS1 occupancy", "TFIIB TSS1 occupancy start coordinate", "TFIIB TSS1 occupancy end coordinate", "MHS TFIIB TSS2 occupancy", "HS3 TFIIB TSS2 occupancy", "TFIIB TSS2 occupancy start coordinate", "TFIIB TSS2 occupancy end coordinate", "motif1 start coordinate", "motif1 end coordinate", "motif1 strand", "motif1 sequence", "motif2 start coordinate", "motif2 end coordinate", "motif2 strand", "motif2 sequence", "MHS motif1 Hsf1 occupancy", "HS3 motif1 Hsf1 occupancy", "motif1 Hsf1 occupancy start coordinate", "motif1 Hsf1 occupancy end coordinate", "MHS motif2 Hsf1 occupancy", "HS3 motif2 Hsf1 occupancy", "motif2 Hsf1 occupancy start coordinate", "motif2 Hsf1 occupancy end coordinate", "MHS motif1 Spt3 occupancy", "HS3 motif1 Spt3 occupancy", "motif1 Spt3 occupancy start coordinate", "motif1 Spt3 occupancy end coordinate", "MHS motif2 Spt3 occupancy", "HS3 motif2 Spt3 occupancy", "motif2 Spt3 occupancy start coordinate", "motif2 Spt3 occupancy end coordinate")

assignment_dict = {}
with open("assignment_list.tsv") as assignments:
	for line in assignments:
		line = line.strip().split()
		geneid = line[0]
		assignment_dict[geneid] = line
	assignments.close()

xu_dict = {}
with open("../shared_files/Xu_2009_ORF_TSS_TES_V64.gff") as xu:
	for line in xu:
		line = line.strip().split()
		geneid = line[8]
		xu_dict[geneid] = line
	xu.close()

seq_dict = {}
with open("motifs.fa") as fasta:
	for line in fasta:
		line = line.strip()
		if line[0] == ">":
			line = line.split(">")[1].split(":")
			chrom = line[0]
			coords = line[1].split("-")
			start = int(coords[0]) + 1	#fasta is off by 1
			end = coords[1]
			key = "{}:{}-{}".format(chrom, str(start), end)
		else:
			seq_dict[key] = line
	fasta.close()

gene_class_dict = fill_dict("gene_classes.tsv")
divergent_dict = fill_dict("divergent_list.tsv")
tfiib_mhs_occupancy_dict = fill_dict("TFIIB_MHS_occupancy.tsv")
tfiib_hs3_occupancy_dict = fill_dict("TFIIB_HS3_occupancy.tsv")
hsf1_mhs_occupancy_dict = fill_dict("Hsf1_MHS_occupancy.tsv")
hsf1_hs3_occupancy_dict = fill_dict("Hsf1_HS3_occupancy.tsv")
spt3_mhs_occupancy_dict = fill_dict("Spt3_MHS_occupancy.tsv")
spt3_hs3_occupancy_dict = fill_dict("Spt3_HS3_occupancy.tsv")

geneids = assignment_dict.keys()

table = np.zeros((len(geneids), len(header)), dtype=object)

for i, geneid in enumerate(geneids):
	#chrom; fill from Xu
	chrom = xu_dict[geneid][0]
	table[i,0] = chrom
	
	#status; fill from divergent list
	if divergent_dict[geneid] == "None":
		table[i,1] = "Hsf1 Motif_Hsf1-bound_Non-divergent"
	else:
		if assignment_dict[geneid][5] == "None":
			table[i,1] = "One Hsf1 Motif_Hsf1-bound_Divergent Gene Pair"
		else:
			table[i,1] = "Two Hsf1 Motif_Hsf1-bound_Divergent Gene Pair" 

	#TSS1 geneID
	if "%" in geneid:	#remove weird trailing characters
		table[i,2] = geneid.split("%")[0]
	else:
		table[i,2] = geneid

	#TSS1 gene class
	if geneid in gene_class_dict.keys():
		table[i,3] = gene_class_dict[geneid]

	#TSS1 coordinate; fill from Xu
	tss1_coord = int(xu_dict[geneid][3])
	table[i,4] = tss1_coord

	#TSS1 strand; fill from Xu
	table[i,5] = xu_dict[geneid][6] 

	#TSS2 geneID; fill from divergent list
	divergentid = divergent_dict[geneid]
	if "%" in divergentid:	#remove weird trailing characters
		table[i,6] = divergentid.split("%")[0]
	else:
		table[i,6] = divergentid

	#TSS2 gene class
	if divergentid != "None" and divergentid in gene_class_dict.keys():
		table[i,7] = gene_class_dict[divergentid]

	#TSS2 coordinate; fill from Xu
	if divergentid != "None":
		tss2_coord = int(xu_dict[divergentid][3])
		table[i,8] = tss2_coord

	#TSS2 strand; fill from Xu
	if divergentid != "None":
		table[i,9] = xu_dict[divergentid][6]

	#TSS1 TFIIB_MHS_occupancy
	table[i,10] = tfiib_mhs_occupancy_dict[geneid]

	#TSS1 TFIIB_HS3_occupancy
	table[i,11] = tfiib_hs3_occupancy_dict[geneid]

	#TSS1 TFIIB occupancy start; fill from TSS1 coordinate
	table[i,12] = tss1_coord - tfiib_flank

	#TSS1 TFIIB occupancy end; fill from TSS1 coordinate
	table[i,13] = tss1_coord + tfiib_flank

	#TSS2 TFIIB_MHS_occupancy
	if divergentid != "None":
		table[i,14] = tfiib_mhs_occupancy_dict[divergentid]

	#TSS2 TFIIB_HS3_occupancy
	if divergentid != "None":
		table[i,15] = tfiib_hs3_occupancy_dict[divergentid]

	#TSS2 TFIIB occupancy start; fill from TSS2 coordinate
	if divergentid != "None":
		table[i,16] = tss2_coord - tfiib_flank

	#TSS2 TFIIB occupancy end; fill from TSS2 coordinate
	if divergentid != "None":
		table[i,17] = tss2_coord + tfiib_flank

	#motif 1 start coordinate; fill from assignment list
	motif1_start = assignment_dict[geneid][2]
	table[i,18] = motif1_start

	#motif 1 end coordinate; fill from assignment list
	motif1_end = assignment_dict[geneid][3]
	table[i,19] = motif1_end

	#motif 1 strand; fill from assignment list
	table[i,20] = assignment_dict[geneid][4]

	#motif 1 sequence; fill from fasta
	table[i,21] = seq_dict["{}:{}-{}".format(chrom, motif1_start, motif1_end)]

	#motif 2 start coordinate; fill from assignment list
	motif2_start = assignment_dict[geneid][6]
	table[i,22] = motif2_start

	#motif 2 end coordinate; fill from assignment list
	motif2_end = assignment_dict[geneid][7]
	table[i,23] = motif2_end

	#motif 2 strand; fill from assignment list
	table[i,24] = assignment_dict[geneid][8]

	#motif 2 sequence; fill from fasta
	if motif2_start != "None":
		table[i,25] = seq_dict["{}:{}-{}".format(chrom, motif2_start, motif2_end)]

	#motif1 Hsf1_MHS_occupancy
	table[i,26] = hsf1_mhs_occupancy_dict["{}:{}-{}".format(chrom, motif1_start, motif1_end)]

	#motif1 Hsf1_HS3_occupancy
	table[i,27] = hsf1_hs3_occupancy_dict["{}:{}-{}".format(chrom, motif1_start, motif1_end)]

	#motif1 Hsf1 occupancy start coordinate; fill from motif1 start coordinate
	table[i,28] = int(motif1_start) - hsf1_flank	

	#motif1 Hsf1 occupancy end coordinate; fill from motif1 end coordinate
	table[i,29] = int(motif1_end) + hsf1_flank	

	#motif2 Hsf1_MHS_occupancy
	if motif2_start != "None":
		table[i,30] = hsf1_mhs_occupancy_dict["{}:{}-{}".format(chrom, motif2_start, motif2_end)]

	#motif2 Hsf1_HS3_occupancy
	if motif2_start != "None":
		table[i,31] = hsf1_hs3_occupancy_dict["{}:{}-{}".format(chrom, motif2_start, motif2_end)]

	#motif2 Hsf1 occupancy start coordinate; fill from motif2 start coordinate
	if motif2_start != "None":
		table[i,32] = int(motif2_start) - hsf1_flank	

	#motif2 Hsf1 occupancy end coordinate; fill from motif2 end coordinate
	if motif2_start != "None":
		table[i,33] = int(motif2_end) + hsf1_flank	

	#motif1 Spt3_MHS_occupancy
	table[i,34] = spt3_mhs_occupancy_dict["{}:{}-{}".format(chrom, motif1_start, motif1_end)]

	#motif1 Spt3_HS3_occupancy
	table[i,35] = spt3_hs3_occupancy_dict["{}:{}-{}".format(chrom, motif1_start, motif1_end)]

	#motif1 Spt3 occupancy start coordinate; fill from motif1 start coordinate
	table[i,36] = int(motif1_start) - spt3_flank

	#motif1 Spt3 occupancy end coordinate; fill from motif1 end coordinate
	table[i,37] = int(motif1_end) + spt3_flank

	#motif2 Spt3_MHS_occupancy
	if motif2_start != "None":
		table[i,38] = spt3_mhs_occupancy_dict["{}:{}-{}".format(chrom, motif2_start, motif2_end)]

	#motif2 Spt3_HS3_occupancy
	if motif2_start != "None":
		table[i,39] = spt3_hs3_occupancy_dict["{}:{}-{}".format(chrom, motif2_start, motif2_end)]

	#motif2 Spt3 occupancy start coordinate; fill from motif2 start coordinate
	if motif2_start != "None":
		table[i,40] = int(motif2_start) - spt3_flank	

	#motif2 Spt3 occupancy end coordinate; fill from motif2 end coordinate
	if motif2_start != "None":
		table[i,41] = int(motif2_end) + spt3_flank

with open("Table_S1.tsv", "w") as out:
	out.write("\t".join(header))
	out.write("\n")
	for row in table:
		out.write("\t".join([str(val) for val in row]))
		out.write("\n")
	out.close()

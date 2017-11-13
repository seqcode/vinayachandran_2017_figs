import sys

gff_path = sys.argv[1]
genes_list_path = sys.argv[2]

gene_info = {}
with open(gff_path) as gff:
	for line in gff:
		gene = line.strip().split()[8]
		gene_info[gene] = line
	gff.close()

with open(genes_list_path[0:len(genes_list_path)-3] + "gff", "w") as out:
	with open(genes_list_path) as genes_list:
		for line in genes_list:
			gene = line.strip()
			out.write(gene_info[gene])
		genes_list.close()
	out.close()

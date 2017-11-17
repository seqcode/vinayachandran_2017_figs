import argparse

def extract_occupancy(tab_file, gff_file, genome_file, out_file, upstream, downstream):
	#initialize tag counts with zeroes
	tags = {}
	with open(genome_file) as chrom_sizes:
		for line in chrom_sizes:
			line = line.strip().split()
			chrom = line[0]
			for i in range(1, int(line[1])):
				key = "{}:{}".format(chrom, i)
				tags[key] = 0
		chrom_sizes.close()

	#fill with tags at each genomic coordinate 
	with open(tab_file) as tab:
		for line in tab:
			if line[0] != "#":
				line = line.strip().split()
				if line[0] != "chrom":
					chrom = line[0]
					index = line[1]
					key = "{}:{}".format(chrom, index)
					tags[key] = float(line[4]) 
		tab.close()

	#get tag counts at each genomic feature
	with open(out_file, "w") as out:
		with open(gff_file) as gff:
			for line in gff:
				line = line.split()
				chrom = line[0]
				coord = int(line[3])
				start = coord - upstream
				end = coord + downstream
				count = sum([tags["{}:{}".format(chrom, j)] for j in range(start,end)])
				out.write(str(count) + "\n")
			gff.close()
		out.close()

def main():
	parser = argparse.ArgumentParser(description="Get tag occupancy in window surrounding each genomic feature.")
	parser.add_argument("tab_file", help="path to tab file with tag counts")
	parser.add_argument("gff_file", help="path to GFF file with genome features")
	parser.add_argument("genome_file", help="path to file with chromosome sizes")
	parser.add_argument("out_file", help="path to file to write to")
	parser.add_argument("upstream", type=int, help="number of bp upstream of feature")
	parser.add_argument("downstream", type=int, help="number of bp downstream of feature")
	
	args = parser.parse_args()

	extract_occupancy(args.tab_file, args.gff_file, args.genome_file, args.out_file, args.upstream, args.downstream)

if __name__ == "__main__":
	main()

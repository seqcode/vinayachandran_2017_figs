import sys

with open(sys.argv[1]) as gff:
	with open(sys.argv[2], "w") as out:
		for line in gff:
			if not line.startswith("#"):	#skip header
				line_array = line.strip().split()
				stddev = float(line_array[8].split(";")[2].split("=")[1])
				if stddev > 0:
					out.write(line)
		out.close()
	gff.close()

# 2016/10/09 Bongso Park
# Original script from Rohit Reja
# Center for Eukaryotic Gene Regulation
# The Penn State University
# Edited by Lila Rieber - 9/18/17

from itertools import tee, izip
import os
import argparse
import sys

# Process file script will process the tab, or idx files from the working folder.
def process_file(idxDir, genFile, stranded, outdir, tmpdir):
	
	full_data = {}
	filenames = []
	# Reading all the tab files and computing the bottom 25% quantile.
	for fname in os.listdir(idxDir):
		# TAB, IDX, or SCIDX file formats are required.
		if not fname.endswith(".tab") or fname.endswith(".idx") or fname.endswith(".scidx"):
			continue
		
		
		idxData = {}
		filenames.append(fname)
		print "[INFO]: Reading the index file = "+fname
		with open(os.path.join(idxDir,fname)) as IN:
			for line in IN:
				if line.startswith("#") or line.startswith("chrom"):
					continue
				cols = line.rstrip().split("\t")
				if stranded:
					idxData[cols[0]+":"+cols[1]+"_F"] = float(cols[2])
					idxData[cols[0]+":"+cols[1]+"_R"] = float(cols[3])
				else: 
					idxData[cols[0]+":"+cols[1]] = float(cols[2]) + float(cols[3])	
		IN.close()

		print "[INFO]: Creating the genomic intervals for = "+fname
		if stranded: 
			outtmp_forward = open(os.path.join(tmpdir,fname+"_forward"),"w")
			outtmp_reverse = open(os.path.join(tmpdir,fname+"_reverse"),"w")
			tmpfile_forward = os.path.join(tmpdir,os.path.splitext(fname)[0]+"_forward.tmp")
			tmpfile_reverse = os.path.join(tmpdir,os.path.splitext(fname)[0]+"_reverse.tmp")
		else:
			outtmp = open(os.path.join(tmpdir,fname),"w")
			tmpfile = os.path.join(tmpdir,os.path.splitext(fname)[0]+"_.tmp")	
		with open(genFile) as ING:
			nooflines = 0
			for line in ING:
				chrom, END = line.rstrip().split("\t")
				for j in list(xrange(1, int(END)+1)):
					key =  chrom+":"+str(j)
					if stranded:
						key_F =  chrom+":"+str(j)+"_F"
						key_R =  chrom+":"+str(j)+"_R"
						if idxData.has_key(key_F) and idxData.has_key(key_R):
							outtmp_forward.write(key+"\t"+str(idxData[chrom+":"+str(j)+"_F"])+"\n")
							outtmp_reverse.write(key+"\t"+str(idxData[chrom+":"+str(j)+"_R"])+"\n")
						else:
							outtmp_forward.write(key+"\t0\n")
							outtmp_reverse.write(key+"\t0\n")
					else:
						if key in idxData:
							outtmp.write(key+"\t"+str(idxData[chrom+":"+str(j)])+"\n")	
						else:
							outtmp.write(key+"\t0\n")
					nooflines += 1
				 
		ING.close()

		if stranded: 
			outtmp_forward.close()
			outtmp_reverse.close()
		else:
			outtmp.close()
	
		print "Sorting the file"
		if stranded: 
			os.system("sort -nrk 2 "+os.path.join(tmpdir,fname+"_forward")+" >"+tmpfile_forward)
			os.system("sort -nrk 2 "+os.path.join(tmpdir,fname+"_reverse")+" >"+tmpfile_reverse)
		else:
			os.system("sort -nrk 2 "+os.path.join(tmpdir,fname)+" >"+tmpfile)	
		print "Completed Sorting !"

	totals = [0] * nooflines
	if stranded: 
		strand = ["forward","reverse"]
	else:
		strand = [""]
	#for fname in os.listdir(tmpdir):
	#	if not fname.endswith(".tmp"):
	#		continue
	#	IN = open(os.path.join(tmpdir,fname))
	#	list1 = []
	#	for line in IN:
	#		cols = line.rstrip().split("\t")
	#		list1.append(float(cols[1]))
	#	IN.close()
	#	totals = map(sum, izip(list1, totals))
	
	for ele in strand:
		for fname in os.listdir(tmpdir):
			if not fname.endswith(ele+".tmp"):
				continue
			IN = open(os.path.join(tmpdir,fname))
			list1 = []
			for line in IN:
				cols = line.rstrip().split("\t")
				list1.append(float(cols[1]))
			IN.close()
			totals = map(sum, izip(list1, totals))
	
		for fname in os.listdir(tmpdir):
			if not fname.endswith("_"+ele+".tmp"):
				continue
			IN = open(os.path.join(tmpdir,fname))
			final_name = os.path.join(outdir,os.path.splitext(fname)[0]+".tab")
			out = open(os.path.join(outdir,os.path.splitext(fname)[0]+".tmp"),"w")
			rank = 0
			for line in IN:
				cols = line.rstrip().split("\t")
				if float(cols[1]) > 0:
					chrom = cols[0].split(":")[0]
					start = cols[0].split(":")[1]
					#Quantile normalized values
					val = float(totals[rank])/len(filenames)
					if ele == "forward":
						out.write(chrom+"\t"+start+"\t"+str(val)+"\t0\t"+str(val)+"\n")
					else:
						out.write(chrom+"\t"+start+"\t0\t"+str(val)+"\t"+str(val)+"\n")
					rank = rank + 1
				else:
					break
			out.close()
			IN.close()
			os.system("sort -k1,1 -k2,2n "+os.path.join(outdir,os.path.splitext(fname)[0]+".tmp")+" >"+final_name)
			os.system("rm "+os.path.join(outdir,fname))
	os.system("rm -r " + tmpdir)


def windows(iterable, size):
	iters = tee(iterable, size)
	for i in xrange(1, size):
		for each in iters[i:]:
			next(each, None)
	return izip(*iters)


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('input_directory',
					  help='Path to directory containing *.tab files')
	parser.add_argument('gen_file',
					  help='The file containing all chromosome number and length.')
	parser.add_argument('--stranded', action="store_true", 
					   help="use strand information")
	 
	args = parser.parse_args()
		
	outdir = os.path.join(args.input_directory,"Normalized_tab_files")
	if not os.path.exists(outdir): os.makedirs(outdir)
	
	tmpdir = os.path.join(args.input_directory,"tmpdir")
	if not os.path.exists(tmpdir): os.makedirs(tmpdir)
	
	process_file(args.input_directory, args.gen_file, args.stranded, outdir, tmpdir)
	
	
if __name__ == "__main__":
	run() 

# Composite plot script by Rohit Reja
# Pugh Lab, Center for Eukaryotic Gene Regulation
# Edited by Lila Rieber - 11/6/17

# Import libraries
import os
import argparse
from pylab import *
import numpy as np
import matplotlib as plt

list1 = {}
## color schema Dark-RED-ROY-G-BIV Black, 
colors = ["#FF3333","#FF9933","#FFFF00","#4C9900","#0000FF","#4B0082","#9F00FF","#000000","#FF3333","#FF9933","#FFFF00","#4C9900","#0000FF","#4B0082","#9F00FF","#000000","#FF3333","#FF9933","#FFFF00","#4C9900",
"#0000FF","#4B0082","#9F00FF","#000000","#FF3333","#FF9933","#FFFF00","#4C9900","#0000FF","#4B0082","#9F00FF","#000000","#FF3333","#FF9933","#FFFF00","#4C9900","#0000FF","#4B0082","#9F00FF","#000000",
"#FF3333","#FF9933","#FFFF00","#4C9900","#0000FF","#4B0082","#9F00FF","#000000"]

def process_file(infile):
	print "processing "+infile
	X = []
	noL = 0
	with open(infile) as in_sense:
		for line in in_sense:
			if line.startswith("Uniqe") or line.startswith("gene"):
				tmp = line.rstrip().split("\t")[2:]
				X = [float(x) for x in tmp]
				xmin = min(X)
				xmax = max(X)
				Y = [0]*len(X)
				continue
		
			noL += 1
			tmplist = line.rstrip().split("\t")[2:]
			newList = [float(x) for x in tmplist]
			Y = map(add,Y,newList)
		in_sense.close()

	return X, Y, xmin, xmax, noL	

def plot_graph(X, Y1, Y2, xmin, xmax, window_size, y, shaded, normalize, ax, label, count, noL):
	if shaded:
		# matplotlib v2.0
		plt.rcParams["font.family"] = "sans-serif"
		plt.rcParams["font.serif"] = "Ubuntu"
		plt.rcParams["font.monospace"] = "Ubuntu Mono"
		plt.rcParams["font.size"] = 10
		plt.rcParams["axes.labelsize"] = 10
		plt.rcParams["axes.labelweight"] = "bold"
		plt.rcParams["xtick.labelsize"] = 8
		plt.rcParams["ytick.labelsize"] = 8
		plt.rcParams["legend.fontsize"] = 10
		plt.rcParams["figure.titlesize"] = 12

	X = movingaverage(X, window_size)
	Y1 = movingaverage(Y1, window_size)
	if Y2 is not None:
		Y2 = movingaverage(Y2, window_size)
	if normalize:
		Y1 = [float(x)/max(Y1) for x in Y1]
		if Y2 is not None:
			Y2 = [float(x)/max(Y2) for x in Y2]

	if shaded:
		ax.plot(X, Y1, color=colors[count],label=label,lw=3.0, zorder=2)
		ax.fill_between(X,Y1,0,facecolor=colors[count],edgecolor=colors[count])	
	else: 
		ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)

	if Y2 is not None:
		if shaded:	 
			ax.plot(X, -np.array(Y2), color=colors[count],label=label,lw=3.0, zorder=2)
			ax.fill_between(X,-np.array(Y2),0,facecolor=colors[count],edgecolor=colors[count])	
		else: 
			ax.plot(X, -np.array(Y2), color=colors[count],label=label,lw=3.0)

	if y is not None:
		if Y2 is None:
			ax.set_ylim(0, y)
		else:
			ax.set_ylim(-y, y)
	
# Moving Average smoothing
def movingaverage(interval, window_size):
	window= np.ones(int(window_size))/float(window_size)
	return np.convolve(interval, window, "valid")

def run():
	parser = argparse.ArgumentParser()
	parser.add_argument("input_directory",
						help="Path to directory containing *.tab files")
	parser.add_argument("-w", type=int, default=5,
						help="Window size of moving average, default=5")
	parser.add_argument("-y", action="store", type=int,
					  help="y-axis bound")
	parser.add_argument("--shaded", action="store_true", 
				help="Creates shaded plot")
	parser.add_argument("--normalize", action="store_true", 
				help="Divides by max to normalize")

	args = parser.parse_args()

	antisense_files = []
	sense_files = []
	unstranded_files = []
	
	output_folder = os.path.join(args.input_directory,"_composite/") 
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	outfile = os.path.join(output_folder,"composite_plot_all_factors.svg")
	
	if not os.path.exists(args.input_directory):
		parser.error("Path {} does not exist.".format(args.input_directory))
	if os.path.isdir(args.input_directory):
		for fname in os.listdir(args.input_directory):
			if fname.endswith(".cdt"):
				fpath = os.path.join(args.input_directory, fname)
				if fpath.endswith("_forward.cdt"):
					sense_files.append(fpath)
				elif fpath.endswith("_reverse.cdt"): 
					antisense_files.append(fpath)
				else:
					unstranded_files.append(fpath)
		
		sense_files.sort()
		antisense_files.sort()
			
		# Declaring plotting parameters
		f, ax = subplots(1,1,sharex="all")
		if args.shaded:
			 f.subplots_adjust(hspace=0)
			
		for i, unstranded_file in enumerate(unstranded_files):
			X, Y, xmin, xmax, noL = process_file(unstranded_file)
			label = os.path.basename(unstranded_file).split("_")[0]
			plot_graph(X, Y, None, xmin, xmax, args.w, args.y, args.shaded, args.normalize, ax, label, i, noL)

		for i, (sense_file, antisense_file) in enumerate(zip(sense_files, antisense_files)):
			X, Y1, xmin1, xmax1, noL = process_file(sense_file)
			X, Y2, xmin2, xmax2, noL = process_file(antisense_file)
			xmin = min((xmin1, xmin2))
			xmax = max((xmax1, xmax2))
			label = os.path.basename(sense_file).split("_")[0]
			plot_graph(X, Y1, Y2, xmin, xmax, args.w, args.y, args.shaded, args.normalize, ax, label, i, noL)

		ax.legend()
		savefig(outfile)
	   
# Execute the main function -> run() 
if __name__ == "__main__":
	run() 

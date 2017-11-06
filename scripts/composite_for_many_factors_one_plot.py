import sys, os, re, difflib
from operator import add
from optparse import OptionParser, IndentedHelpFormatter
from pylab import *
import numpy as np
from matplotlib.font_manager import FontProperties

### color schema for Figure 2A
colors = ["#FBDB0C","#FF00FF","#0000FF","#336600","#FF0000","#000000","#00FF00","#660066","#FBDB0C","#FF00FF","#0000FF","#336600","#FF0000","#000000","#00FF00","#660066","#FBDB0C","#FF00FF","#0000FF","#336600",
"#FF0000","#000000","#00FF00","#660066","#FBDB0C","#FF00FF","#0000FF","#336600","#FF0000","#000000","#00FF00","#660066","#FBDB0C","#FF00FF","#0000FF","#336600","#FF0000","#000000","#00FF00","#660066"] 

list1 = []
def  process_files(sense,anti,options,output_folder,ax,count,color_count):
    print "processing "+sense+" and "+anti
    
    label = os.path.basename(sense).split("_")[0]
    X =[]
    # Process Sense file first
    in_sense = open(sense)
    for line in in_sense:
        if line.startswith("Uniqe"):
            tmp = line.rstrip().split("\t")[2:]
            X = [int(x) for x in tmp]
            xmin = min(X)
            xmax = max(X)
            Y1 = [0]*len(X)
            continue
        if line.startswith("EWEIGHT"):
            continue
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y1 = map(add,Y1,newList)
    
    # Process anti-sense file
    in_anti = open(anti,"rt")
    for line in in_anti:
        if line.startswith("Uniqe"):
            Y2 = [0]*len(X)
            continue
        if line.startswith("EWEIGHT"):
            continue
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y2 = map(add,Y2,newList)
    
    # Uncomment here to normalize by number of genes.
    #Y1 = [float(x)/max(Y1) for x in Y1]
    #Y2 = [float(x)/max(Y2) for x in Y2]

    #### To print out the X and Y coordiantes.
    out1 = open(os.path.join(output_folder,label+"_sense.txt"),"w")
    out2 = open(os.path.join(output_folder,label+"_anti.txt"),"w")
    
    for x,y in zip(X,Y1):
        if x >= -200 and x <= 200:
            out1.write(str(x)+"\t"+str(y)+"\n")
            
    for x,y in zip(X,Y2):
        if x >= -200 and x <= 200:
            out2.write(str(x)+"\t"+str(y)+"\n")
            
    
    
    plot_graph(X,Y1,Y2,xmin,xmax,options,ax,label,count,color_count)

def process_onestrand_files(infile,options,output_folder,ax,count,color_count):
    print "processing "+infile
    X = []
    label = os.path.basename(infile).split("_")[0]
    # Process the only CDT file 
    in_sense = open(infile,"rt")
    line_count = 0
    for line in in_sense:
        if line.startswith("Uniqe") or line.startswith("ID") or line.startswith("gene"):
            tmp = line.rstrip().split("\t")[2:]
            X = [int(x) for x in tmp]
            xmin = min(X)
            xmax = max(X)
            Y = [0]*len(X)
            continue
        if line.startswith("EWEIGHT") or line.startswith("chr1:1-100"):
            continue
        
        tmplist = line.rstrip().split("\t")[2:]
        newList = [float(x) for x in tmplist]
        Y = map(add,Y,newList)
        line_count = line_count + 1
    
    ####### Create excel files for each factor #######
    out3 = open(os.path.join(output_folder,label+".txt"),"w")
    for x,y in zip(X,Y):
        if x >= -500 and x <= 500:
            out3.write(str(x)+"\t"+str(y)+"\n")
    #######
    plot_graph(X,Y,0,xmin,xmax,options,ax,label,count,color_count)
    

def plot_graph(X,Y1,Y2,xmin,xmax,options,ax,label,count,color_count):
    
    # Converting the arrays to numpy array
    X = np.array(X)
    Y1 = np.array(Y1)
    if not Y2 == 0:
        
        Y2 = np.array(Y2)
        ## Smoothing using moving average, REQUIRED WHEN PLOTTING COMPOSITE ON +Y AND  -Y
        X = movingaverage(X,options.window)
        Y1 = movingaverage(Y1,options.window)
        Y2 = movingaverage(Y2,options.window)
        #Y1 = [float(x)/max(Y1) for x in Y1]
        #Y2 = [float(x)/max(Y2) for x in Y2]
        Y2 = [-1*x for x in Y2]
       
        ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)
        ax.plot(X, Y2, color=colors[count],lw=3.0)

        for tick in ax.yaxis.get_major_ticks():
            tick.label.set_fontsize(12)
            
        for tick in ax.xaxis.get_major_ticks():
            tick.label.set_fontsize(12)
    else:
        Y1 = movingaverage(Y1,options.window)
        #Y1 = [float(x)/max(Y1) for x in Y1]
        X = movingaverage(X,options.window)
        
        ## To print the area under the curve
        #print sum(Y1)
        #Y1 = [float(x)/max(Y1) for x in Y1]
        
        
        #index = Y1.index(max(Y1)) - 500
        #print "Max value at "+str(index)
        #X = smoothListGaussian(X,options.window)
        #Y1 = smoothListGaussian(Y1,options.window)
        
        if label == "Rap1" or label == "Normal":
            ax.plot(X, Y1, color="#AAAAAA",label=label,lw=3.0)
            ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
            
        else:
            ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)

    if options.ylim is not None:
       ax.set_ylim(-options.ylim, options.ylim)
    
def movingaverage(interval, window_size):
    # suggestion from: http://argandgahandapandpa.wordpress.com/2011/02/24/python-numpy-moving-average-for-data/
    #extended_interval = np.hstack([[interval[0]] * (window_size- 1), interval])
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')


def smoothListGaussian(list,degree):  

     window=degree*2-1  
     weight=np.array([1.0]*window)  
     weightGauss=[]  
     for i in range(window):  
         i=i-degree+1  
         frac=i/float(window)  
         gauss=1/(np.exp((4*(frac))**2))  
         weightGauss.append(gauss)  
     weight=np.array(weightGauss)*weight  
     smoothed=[0.0]*(len(list)-window)  
     for i in range(len(smoothed)):  
         smoothed[i]=sum(np.array(list[i:i+window])*weight)/sum(weight)
     new_smooth = [int(i) for i in smoothed]
     return new_smooth  





usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python composite_for_many_factors_one_plot.py /usr/local/folder_containing_CDT_files
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-w', action='store', type='int', dest='window', default = 5,
                      help='Window size of moving average., Default=5')
    parser.add_option('-y', action='store', type='int', dest='ylim', default = None,
                      help='y-axis bounds, Default=None (i.e. will be set from data)')

    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    antisense_files = []
    sense_files = []
    
    output_folder = os.path.join(args[0],'_composite/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    outfile = os.path.join(output_folder,"composite_plot_all_factors.svg")
    
    
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    if os.path.isdir(args[0]):
        for fname in os.listdir(args[0]):
            print fname
            if fname.endswith("cdt") or fname.endswith("txt"):
                # intelligently join paths without worrying about '/'
                fpath = os.path.join(args[0], fname)
                print fpath
                if re.match(r'(.*)reverse(.*)',fname):
                    antisense_files.append(fpath)
                elif re.match(r'(.*)forward(.*)',fname):
                    sense_files.append(fpath)
                    
        # No of subplots to plot
        nof = len(sense_files)
         # Declaring plotting parameters
        #f,ax = subplots(1,1,1)
	f,ax = subplots(1,1)	#TEST
        count = -1
        color_count = 0
        
        
        if len(antisense_files) == 0:
            for f in sense_files:
                count += 1
                process_onestrand_files(f,options,output_folder,ax,count,color_count)
                color_count = color_count + 2
            lfp = FontProperties()
            lfp.set_size(12)
            ax.legend(prop=lfp)
            savefig(outfile)
        else:   
            for f1 in antisense_files:
                count += 1
                Oratio = 0
                file1_anti = ""
                file2_sense = ""
                for f2 in sense_files:
                    # Calculating the distance ration between two file names and grouping together based upon the highest ratio
                    # Requirement: Antisense file should have anti word in it and both sense and antisense files should have similar names
                    Nratio = difflib.SequenceMatcher(None,f1,f2).ratio()
                    if(Nratio >= Oratio):
                        file1_anti = f1
                        file2_sense = f2
                        Oratio = Nratio
                    
                process_files(file2_sense,file1_anti,options,output_folder,ax,count,color_count)
                color_count = color_count + 2
            lfp = FontProperties()
            lfp.set_size(12)
            ax.legend(prop=lfp)
            savefig(outfile)
    
    
if __name__ == "__main__":
    run() 

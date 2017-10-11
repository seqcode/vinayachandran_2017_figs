# Composite plot script by Rohit Reja
# Pugh Lab, Center for Eukaryotic Gene Regulation
# Edited by Lila Rieber - 9/15/17

# Import libraries
import sys, os
import argparse
from pylab import *
import numpy as np
import matplotlib as plt
from scipy import stats

list1 = {}
## color schema Dark-RED-ROY-G-BIV Black, 
colors = ['#FF3333','#FF9933','#FFFF00','#4C9900','#0000FF','#4B0082','#9F00FF','#000000','#FF3333','#FF9933','#FFFF00','#4C9900','#0000FF','#4B0082','#9F00FF','#000000','#FF3333','#FF9933','#FFFF00','#4C9900',
'#0000FF','#4B0082','#9F00FF','#000000','#FF3333','#FF9933','#FFFF00','#4C9900','#0000FF','#4B0082','#9F00FF','#000000','#FF3333','#FF9933','#FFFF00','#4C9900','#0000FF','#4B0082','#9F00FF','#000000',
'#FF3333','#FF9933','#FFFF00','#4C9900','#0000FF','#4B0082','#9F00FF','#000000']

def process_onestrand_files(infile, window_size, y, shaded, normalize, output_folder, ax, count):
    print 'processing '+infile
    X = []
    label = os.path.basename(infile).split('_')[0]
    # Process the only CDT file 
    noL = 0
    with open(infile) as in_sense:
        for line in in_sense:
            if line.startswith('Uniqe') or line.startswith('ID') or line.startswith('gene'):
                tmp = line.rstrip().split('\t')[2:]
                X = [int(x) for x in tmp]
                xmin = min(X)
                xmax = max(X)
                Y = [0]*len(X)
                continue
        
            noL += 1
            tmplist = line.rstrip().split('\t')[2:]
            newList = [float(x) for x in tmplist]
            Y = map(add,Y,newList)
        in_sense.close()
        
    list1[label] = smoothListGaussian(Y, window_size)
    ##### REMOVE THE HASH WHEN YOU WANT TO DIVIDE BY NO OF GENES.
    #Y = [float(x)/noL for x in Y]
    plot_graph(X, Y, 0, xmin, xmax, window_size, y, shaded, normalize, ax, label, count, noL)
    

def plot_graph(X, Y1, Y2, xmin, xmax, window_size, y, shaded, normalize, ax, label, count, noL):
    if shaded:
        # matplotlib v2.0
        plt.rcParams['font.family'] = 'sans-serif'
        plt.rcParams['font.serif'] = 'Ubuntu'
        plt.rcParams['font.monospace'] = 'Ubuntu Mono'
        plt.rcParams['font.size'] = 10
        plt.rcParams['axes.labelsize'] = 10
        plt.rcParams['axes.labelweight'] = 'bold'
        plt.rcParams['xtick.labelsize'] = 8
        plt.rcParams['ytick.labelsize'] = 8
        plt.rcParams['legend.fontsize'] = 10
        plt.rcParams['figure.titlesize'] = 12

    X = movingaverage(X, window_size)
    Y1 = movingaverage(Y1, window_size)
    if normalize:
        Y1 = [float(x)/max(Y1) for x in Y1]

    if shaded:     
        if label == 'Hsf-MHS':	
            ax.plot(X, Y1, color='#AAAAAA',label=label,lw=3.0, zorder=1)	
            ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
        else:
            ax.plot(X, Y1, color=colors[count],label=label,lw=3.0, zorder=2)	
    else: 
        if label == 'Nap1-wt-140-180-Rep2' or label == 'Heatshock': 
            ax.plot(X, Y1, color='#AAAAAA',label=label,lw=3.0)
            ax.fill_between(X,Y1,0,facecolor='#AAAAAA',edgecolor='#AAAAAA')
            ax.set_ylim(0,1)
        else:
            ax.plot(X, Y1, color=colors[count],label=label,lw=3.0)
    if y is not None:
       ax.set_ylim(0, y)
    
# Moving Average smoothing
def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'valid')

# Gaussian smoothing
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

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_directory',
                        help='Path to directory containing *.tab files')
    parser.add_argument('-w', type=int, default=5,
                        help='Window size of moving average, default=5')
    parser.add_argument('-y', action='store', type=int, default = None,
                      help='y-axis bound, Default=None (i.e. will be set from data)')
    parser.add_argument('--shaded', action='store_true', 
		        help='Creates shaded plot')
    parser.add_argument('--normalize', action='store_true', 
		        help='Divides by max to normalize')

    args = parser.parse_args()

    antisense_files = []
    sense_files = []
    
    output_folder = os.path.join(args.input_directory,'_composite/') 
    if not os.path.exists(output_folder): os.makedirs(output_folder)
    outfile = os.path.join(output_folder,'composite_plot_all_factors.svg')
    
    if not os.path.exists(args.input_directory):
        parser.error('Path {} does not exist.'.format(args.input_directory))
    if os.path.isdir(args.input_directory):
        for fname in os.listdir(args.input_directory):
            if fname.endswith('.txt') or fname.endswith('.cdt'):
                # intelligently join paths without worrying about '/'
                fpath = os.path.join(args.input_directory, fname)
                sense_files.append(fpath)
                    
        # The number of the subplots
        nof = len(sense_files)
        
        # Declaring plotting parameters
        f,ax = subplots(1,1,sharex='all')
        if args.shaded:
             f.subplots_adjust(hspace=0)
        count = -1        

    for f in sense_files:
        count += 1
        process_onestrand_files(f, args.w, args.y, args.shaded, args.normalize, output_folder, ax, count)
    for k1,v1 in list1.items():
        for k2,v2 in list1.items():
            if k1 == k2:
                continue
            print 'pvalue of KS-test between {} and {} = {}'.format(k1, k2, stats.ks_2samp(v1,v2))
   
    if args.shaded:
        csfont = {'fontname':'Arial MS'}
        ax.legend(loc=1,prop={'size':12})
    else:
        ax.legend()
    savefig(outfile)
       
# Execute the main function -> run() 
if __name__ == '__main__':
    run() 

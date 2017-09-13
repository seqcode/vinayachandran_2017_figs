# 2016/10/09 Bongso Park
# Original script from Rohit Reja
# Center for Eukaryotic Gene Regulation
# The Penn State University
from itertools import tee, izip, islice
from operator import itemgetter
import itertools
import sys, os, math, operator
from optparse import OptionParser , IndentedHelpFormatter
from collections import OrderedDict

# Process file script will process the tab, or idx files from the working folder.
def process_file(idxDir,options,outdir,tmpdir):
    
    full_data = {}
    filenames = []
    ## Reading all the tab files and computing the bottom 25% quantile.
    for fname in os.listdir(idxDir):
        # TAB, IDX, or SCIDX file formats are required.
        if not fname.endswith(".tab") or fname.endswith(".idx") or fname.endswith(".scidx"):
            continue
        
        idxData = {}
        filenames.append(fname)
        print "[INFO]: Reading the index file = "+fname
        IN = open(os.path.join(idxDir,fname),"rt")
        for line in IN:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            idxData[cols[0]+":"+cols[1]] = float(cols[2]) + float(cols[3])
        IN.close()
        
        print "[INFO]: Creating the genomic intervals for = "+fname
        outtmp = open(os.path.join(tmpdir,fname),"w")
        tmpfile = os.path.join(tmpdir,os.path.splitext(fname)[0]+".tmp")
        ING = open(options.genFile,"rt")
        nooflines = 0
        for line in ING:
            chrom, END = line.rstrip().split("\t")
            for j in list(xrange(1, int(END)+1)):
                key =  chrom+":"+str(j)
                if key in idxData:
                    #intervals[key] = idxData[chrom+":"+str(j)]
                    outtmp.write(key+"\t"+str(idxData[chrom+":"+str(j)])+"\n")
                else:
                    #intervals[key] = 0
                    outtmp.write(key+"\t0\n")
                nooflines = nooflines + 1
                 
        ING.close()
        outtmp.close()
        print "Sorting the file"
        os.system("sort -nrk 2 "+os.path.join(tmpdir,fname)+" >"+tmpfile)
        print "Completed Sorting !"
        
    totals = [0] * nooflines
    for fname in os.listdir(tmpdir):
        if not fname.endswith(".tmp"):
            continue
        IN = open(os.path.join(tmpdir,fname),"rt")
        list1 = []
        for line in IN:
            cols = line.rstrip().split("\t")
            list1.append(float(cols[1]))
        IN.close()
        totals = map(sum, itertools.izip(list1, totals))
    
    for fname in os.listdir(tmpdir):
        if not fname.endswith(".tmp"):
            continue
        IN = open(os.path.join(tmpdir,fname),"rt")
        final_name = os.path.join(outdir,os.path.splitext(fname)[0]+".tab")
        out = open(os.path.join(outdir,os.path.splitext(fname)[0]+".tmp"),"w")
        rank = 0
        for line in IN:
            cols = line.rstrip().split("\t")
            if float(cols[1]) > 0:
                chrom = cols[0].split(":")[0]
                start = cols[0].split(":")[1]
                val = float(totals[rank])/len(filenames)
                out.write(chrom+"\t"+start+"\t0\t"+str(val)+"\t"+str(val)+"\n")
                rank = rank + 1
            else:
                break
        out.close()
        IN.close()
        os.system("sort -k1,1 -k2,2n "+os.path.join(outdir,os.path.splitext(fname)[0]+".tmp")+" >"+final_name)
    os.system("rm -fr "+tmpdir)
    os.system("rm "+os.path.join(outdir,os.path.splitext(fname)[0]+".tmp"))
        
            
        

    

def windows(iterable, size):
    iters = tee(iterable, size)
    for i in xrange(1, size):
        for each in iters[i:]:
            next(each, None)
    return izip(*iters)




usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python quantile_normalization.py [OPTIONS] /path/to/idx/directory/

# The total reads in the Input should be less than the total reads in
the treatment.
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-g', action='store', type='string', dest='genFile',
                      help='The file containing all chromosome no and length.')
    #parser.add_option('-s', action='store', type='int', dest='window',default=5,
    #                  help='Size of the windows to take across the genome., Default=5')
    #parser.add_option('-q', action='store', type='int', dest='quant',default=25,
    #                  help='The quantile to consider for normalization. default=25')
    #parser.add_option('-f', action='store', type='int', dest='flag',default = 1,
    #                  help='for HS/MHS sample = 0/ For Treatment/control sample = 1.Default = 1')
    
    (options, args) = parser.parse_args()
    
        
    outdir = os.path.join(args[0],"Normalized_tab_files")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    tmpdir = os.path.join(args[0],"tmpdir")
    if not os.path.exists(tmpdir): os.makedirs(tmpdir)
    
    process_file(args[0],options,outdir,tmpdir)
    
    
if __name__ == "__main__":
    run() 

# 2016/10/09 Bongso Park
# Original script from Rohit Reja
# Center for Eukaryotic Gene Regulation
# The Penn State University
# Strandness implementation in progress not completed yet
#
# 9/13/2017 - edited by Lila Rieber

from itertools import tee, izip
import os
from optparse import OptionParser, IndentedHelpFormatter

# Process file script will process the tab, or idx files from the working folder.
def process_file(idxDir,options,outdir,tmpdir):
    
    full_data = {}
    filenames = []
    ## Reading all the tab files and computing the bottom 25% quantile.
    for fname in os.listdir(idxDir):
        # TAB, IDX, or SCIDX file formats are required.
        if not fname.endswith(".tab") or fname.endswith(".idx") or fname.endswith(".scidx"):
            continue
        
        # In the updated version, we are using strand separated version.
        idxData = {}
        filenames.append(fname)
        print "[INFO]: Reading the index file = "+fname
        IN = open(os.path.join(idxDir,fname),"rt")
        for line in IN:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            idxData[cols[0]+":"+cols[1]+"_F"] = float(cols[2])
            idxData[cols[0]+":"+cols[1]+"_R"] = float(cols[3])
        IN.close()
        
        print "[INFO]: Creating the genomic intervals for = "+fname
        outtmp_forward = open(os.path.join(tmpdir,fname+"_forward"),"w")
        outtmp_reverse = open(os.path.join(tmpdir,fname+"_reverse"),"w")
        tmpfile_forward = os.path.join(tmpdir,os.path.splitext(fname)[0]+"_forward.tmp")
        tmpfile_reverse = os.path.join(tmpdir,os.path.splitext(fname)[0]+"_reverse.tmp")
        ING = open(options.genFile,"rt")
        nooflines = 0
        for line in ING:
            chrom, END = line.rstrip().split("\t")
            for j in list(xrange(1, int(END)+1)):
                key =  chrom+":"+str(j)
                key_F =  chrom+":"+str(j)+"_F"
                key_R =  chrom+":"+str(j)+"_R"
                if idxData.has_key(key_F) and idxData.has_key(key_R):
                    outtmp_forward.write(key+"\t"+str(idxData[chrom+":"+str(j)+"_F"])+"\n")
                    outtmp_reverse.write(key+"\t"+str(idxData[chrom+":"+str(j)+"_R"])+"\n")
                else:
                    outtmp_forward.write(key+"\t0\n")
                    outtmp_reverse.write(key+"\t0\n")
                nooflines = nooflines + 1
                 
        ING.close()
        outtmp_forward.close()
        outtmp_reverse.close()
        print "Sorting the file"
        os.system("sort -nrk 2 "+os.path.join(tmpdir,fname+"_forward")+" >"+tmpfile_forward)
        print "Completed Sorting - Forward tags!"        
        os.system("sort -nrk 2 "+os.path.join(tmpdir,fname+"_reverse")+" >"+tmpfile_reverse)
        print "Completed Sorting - Reverse tags!"        
    
    totals = [0] * nooflines
    strand = ["forward","reverse"]

    for ele in strand:
        for fname in os.listdir(tmpdir):
            if not fname.endswith(ele+".tmp"):
                continue
            IN = open(os.path.join(tmpdir,fname),"rt")
            list1 = []
            for line in IN:
                cols = line.rstrip().split("\t")
                list1.append(float(cols[1]))
            IN.close()
            totals = map(sum, izip(list1, totals))
    
        for fname in os.listdir(tmpdir):
            if not fname.endswith("_"+ele+".tmp"):
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
    #leave the temporary files to examine
    os.system("rm "+os.path.join(outdir,os.path.splitext(fname)[0]+"_"+ele+".tmp"))
        

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
    
    (options, args) = parser.parse_args()
    
        
    outdir = os.path.join(args[0],"Normalized_tab_files")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
    tmpdir = os.path.join(args[0],"tmpdir")
    if not os.path.exists(tmpdir): os.makedirs(tmpdir)
    
    process_file(args[0],options,outdir,tmpdir)
    
    
if __name__ == "__main__":
    run() 

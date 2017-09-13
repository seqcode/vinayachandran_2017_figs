from optparse import OptionParser , IndentedHelpFormatter
import sys, os, re, operator
from collections import OrderedDict


def  process_file(cdtDir,options,outdir):
    input = open(options.sortbyFile,"rt")
    order = OrderedDict()
    for line in input:
        if line.startswith("Uniqe") or line.startswith("EWEIGHT") or line.startswith("gene") or line.startswith("#"):
            continue
        cols = line.rstrip().split("\t")
        if options.criteria == 0:
            key = cols[0]
            #key = line.rstrip()
            order[key] = 1
        if options.criteria == 2:
            key = cols[8]
            order[key] = 1
        elif options.criteria == 1:
            key = cols[0]+":"+cols[3]+"-"+cols[4]
            order[key] = 1
    input.close()
    
    
    
    for fname in os.listdir(cdtDir):
        data = {}
        if not (fname.endswith('.cdt') or fname.endswith('.gff') or fname.endswith('.gtf')):
            continue
        outfile = os.path.join(outdir,os.path.splitext(fname)[0]+".cdt")
        out = open(outfile,"w")
        in1 = open(os.path.join(cdtDir,fname))
        
        for line in in1:
            if line.startswith("Uniqe") or line.startswith("EWEIGHT") or line.startswith("Gene ID"):
                out.write(line)
                continue
            tmp = line.rstrip().split("\t")
            length = len(tmp)
            data[tmp[0]] = line
            
            
        for k,v in order.items():
            try:
                out.write(data[k])
            except KeyError:
                #out.write(k+"\t0\n")
                print "the following key does not exist "+k
                
        
        if options.extra > 0:
            # Create a line with fake coordiantes and all 0 tag counts
            line = "chr1:1-100"
            for i in range(1,length-1):
                line = line+"\t0"
            # Run the loop from 0 to "no of lines to be inserted" and print it at the end of the file.
            for j in range(0,options.extra):
                out.write(line+"\n")
            
        out.close()
        in1.close()
        print "completed file "+fname     
    


usage = '''
input_paths may be:
- a directory to run on all files in them

example usages:
python sort_CDT_by_given_file.py  <path_to_CDT_folder>
'''.lstrip()


 
# We must override the help formatter to force it to obey our newlines in our custom description
class CustomHelpFormatter(IndentedHelpFormatter):
    def format_description(self, description):
        return description


def run():
    parser = OptionParser(usage='%prog [options] input_paths', description=usage, formatter=CustomHelpFormatter())
    parser.add_option('-i', action='store', type='string', dest='sortbyFile',
                      help='Input file which will decide the sorting order')
    parser.add_option('-o', action='store', type='int', dest='criteria',default=0,
                      help='Sorting criteria, 0=> sort by first column of any file, 1=> sort by given chr:start-end order in the file, 2 => last column of gff')
    parser.add_option('-n', action='store', type='int', dest='extra',default=0,
                      help='Extra N lines to add to the end of CDT')
    
    (options, args) = parser.parse_args()
    
    if not args:
        parser.print_help()
        sys.exit(1)
        
    outdir = os.path.join(os.path.dirname(args[0]),"sorted")
    if not os.path.exists(outdir): os.makedirs(outdir)
    
     
    if not os.path.exists(args[0]):
        parser.error('Path %s does not exist.' % args[0])
    
    process_file(args[0],options,outdir)
    
    
    
if __name__ == "__main__":
    run() 

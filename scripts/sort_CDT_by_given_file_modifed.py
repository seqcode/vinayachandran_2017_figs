from optparse import OptionParser , IndentedHelpFormatter
import sys, os

def  process_file(cdtDir,options,outdir):
    input = open(options.sortbyFile,"rt")
    order = []
    for line in input:
        cols = line.rstrip().split("\t")
        order.append(cols[0])
       
    input.close()
    
    
        
    for fname in os.listdir(cdtDir):
        data = {}
        if not (fname.endswith('.cdt') or fname.endswith('.txt')):
            continue
        outfile = os.path.join(outdir,os.path.splitext(fname)[0]+".cdt")
        out = open(outfile,"w")
        in1 = open(os.path.join(cdtDir,fname))
        
        for line in in1:
            if line.startswith("Uniqe") or line.startswith("EWEIGHT") or line.startswith("gene"):
                out.write(line)
                continue
            tmp = line.rstrip().split("\t")
            length = len(tmp)
            data[tmp[0]] = line
            
            
        for k in order:
            try:
                out.write(data[k])
            except KeyError:
                continue

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
    parser.add_option('-n', action='store', type='int', dest='extra',default=50,
                      help='Number of lines to add in between, default = 50')
    
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

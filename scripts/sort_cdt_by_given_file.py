import os
from collections import OrderedDict
import argparse
import sys

def  process_file(args):
    outdir = os.path.join(args.input_directory, 'sorted')
    if not os.path.exists(outdir): os.makedirs(outdir)
    
     
    if not os.path.exists(args.input_directory):
        parser.error('Path {} does not exist.'.format(args.input_directory))

    input = open(args.sort_file)
    order = OrderedDict()
    for line in input:
        if line.startswith('Uniqe') or line.startswith('EWEIGHT') or line.startswith('gene') or line.startswith('#'):
            continue
        cols = line.rstrip().split('\t')
        if args.criteria == 0:
            key = cols[0]
            order[key] = 1
        if args.criteria == 2:
            key = cols[8]
            order[key] = 1
        elif args.criteria == 1:
            key = cols[0]+':'+cols[3]+'-'+cols[4]
            order[key] = 1
    input.close()
    
    for fname in os.listdir(args.input_directory):
        data = {}
        if not (fname.endswith('.cdt') or fname.endswith('.gff') or fname.endswith('.gtf')):
            continue
        outfile = os.path.join(outdir,os.path.splitext(fname)[0]+'.cdt')
        out = open(outfile,'w')
        in1 = open(os.path.join(args.input_directory,fname))
        
        for line in in1:
            if line.startswith('Uniqe') or line.startswith('EWEIGHT') or line.startswith('Gene ID'):
                out.write(line)
                continue
            tmp = line.rstrip().split('\t')
            length = len(tmp)
            data[tmp[0]] = line
        
        for k,v in order.items():
            try:
                out.write(data[k])
            except KeyError:
                #print 'the following key does not exist '+k
                continue
        
        if args.extra > 0:
            # Create a line with fake coordinates and all 0 tag counts
            line = 'chr1:1-100'
            for i in range(1,length-1):
                line = line+'\t0'
            # Run the loop from 0 to 'no of lines to be inserted' and print it at the end of the file.
            for j in range(0,args.extra):
                out.write(line+'\n')
            
        out.close()
        in1.close()
        print 'completed file '+fname     
    

def run():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_directory',
                      help='Path to directory containing *.cdt files')
    parser.add_argument('sort_file',
                      help='Path to file describing sort order')
    parser.add_argument('-n', type=int, dest='extra',default=50,
                     help='Number of lines to add in between, default=50')
    parser.add_argument('-o', type=int, dest='criteria',default=0,
                     help='Sorting criteria, 0=> sort by first column of any file, 1=> sort by given chr:start-end order in the file, 2 => last column of gff')
     
    args = parser.parse_args()

    process_file(args)
    
    
if __name__ == '__main__':
    run() 

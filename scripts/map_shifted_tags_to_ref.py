import os
from collections import defaultdict
import argparse #Lila

def process_file(args):
    outdir = os.path.join(os.path.dirname(args.input_directory),"_CDT")
    if not os.path.exists(outdir): os.makedirs(outdir)
     
    if not os.path.exists(args.input_directory):
       parser.error('Path {} does not exist.'.format(args.input_directory))    

    idxDir=args.input_directory	#Lila

    # Process reference file first:
    in0 = open(args.ref,"rt")
    tmpref = os.path.join(outdir,"tmpref.gff")
    outref = open(tmpref,"w")
    for line in in0:
        if line.startswith("#") or line.startswith("chrom"):
            continue
        cols = line.rstrip().split("\t")
        if cols[6] == "+" or cols[6] == ".":
            start = int(cols[3]) - args.up
            end = int(cols[3]) + args.down
            if start <= 0:
                outref.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t1\t"+str(end)+"\t"+cols[5]+"\t"+cols[6]+"\t"+cols[7]+"\t"+cols[8]+"\n")
            else:
                outref.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t"+str(start)+"\t"+str(end)+"\t"+cols[5]+"\t"+cols[6]+"\t"+cols[7]+"\t"+cols[8]+"\n")
        elif cols[6] == "-":
            start = int(cols[3]) - args.down
            end = int(cols[3]) + args.up
            if start <= 0:
                outref.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t1\t"+str(end)+"\t"+cols[5]+"\t"+cols[6]+"\t"+cols[7]+"\t"+cols[8]+"\n")
            else:
                outref.write(cols[0]+"\t"+cols[1]+"\t"+cols[3]+"\t"+str(start)+"\t"+str(end)+"\t"+cols[5]+"\t"+cols[6]+"\t"+cols[7]+"\t"+cols[8]+"\n")
    outref.close()
    in0.close()
    
    for fname in os.listdir(idxDir):
        
        if not (fname.endswith(".idx") or fname.endswith(".tab")):
            continue
        outcdt = os.path.join(outdir,os.path.splitext(fname)[0]+".cdt")
        tmptab = os.path.join(outdir,"tmptab.gff")
        out = open(tmptab,"w")
        in1 = open(os.path.join(idxDir,fname),"rt")
        print "INFO: Processing = "+fname
        for line in in1:
            if line.startswith("#") or line.startswith("chrom"):
                continue
            cols = line.rstrip().split("\t")
            total = float(cols[2]) + float(cols[3])
            if int(cols[1]) == 0:
                continue
            out.write(cols[0]+"\t.\t.\t"+cols[1]+"\t"+cols[1]+"\t"+str(total)+"\t.\t.\t.\n")
        
        in1.close()
        out.close()
        intersect = os.path.join(outdir,"intersect.txt")
        os.system(args.location+"intersectBed"+" -wao -a "+tmpref+" -b "+tmptab+" >"+intersect)
        
        # Remove tmptab as its not required anymore.
        os.system("rm "+tmptab)
        
        process_intersect_file(intersect,args,outcdt)
        
        # Remove intersect as its not required anymore.
        os.system("rm "+intersect)
        
        
    os.system("rm "+tmpref)

def process_intersect_file(infile,args,outcdt):
    cdt_dict = defaultdict(list)
    in2 = open(infile,"rt")
    out = open(outcdt,"w")
    print_header(out,args)
    for line in in2:
        cols = line.rstrip().split("\t")
        dist = int(cols[12]) - int(cols[2])
        if cols[12] == -1:
            cdt_dict[cols[8]].append("NA")
            continue
        if cols[6] == "+":
            cdt_dict[cols[8]].append(str(dist)+":"+cols[14])
        elif cols[6] == "-":
            # performing flipping here!
            new_dist = (-1)*dist
            cdt_dict[cols[8]].append(str(new_dist)+":"+cols[14])
    in2.close()
    
    for k,v in cdt_dict.items():
        tmpdict = {}
        line = k+"\t1"
        totaltag = 0
        if v == "NA":
            for i in range((-1)*args.up,args.down+1):
                line = line+"\t0"
                totaltag = totaltag + 0
                
        for val in v:
            pos = int(val.split(":")[0])
            tag = (val.split(":"))[1]
            tmpdict[pos] = tag
        for i in range((-1)*args.up,args.down+1):
            if i in tmpdict:
                line = line+"\t"+tmpdict[i]
                totaltag = totaltag + float(tmpdict[i])
                
            else:
                line = line+"\t0"
        # Comment this line for total tags.
        out.write(line+"\n")
    
    out.close()

def print_header(out,args):
    line = "Uniqe ID\tGWEIGHT"
    for i in range((-1)*args.up,args.down+1):
        line = line+"\t"+str(i)
    out.write(line+"\n")

def run():
     #Lila
     parser = argparse.ArgumentParser()
     parser.add_argument('input_directory',
                      help='Path to directory containing *.tab files')
     parser.add_argument('ref',
                      help='Stranded reference file in gff format.')
     parser.add_argument('-u', type=int, dest='up',default=500,
                     help='Upstream distance, default=500')
     parser.add_argument('-d', type=int, dest='down',default=500,
                      help='Downstream distance, default=500')
     parser.add_argument('-l', dest='location',default='',
                      help='BEDTools location (if not in path)')
     
     args = parser.parse_args()
        
     process_file(args)
    
if __name__ == "__main__":
    run() 

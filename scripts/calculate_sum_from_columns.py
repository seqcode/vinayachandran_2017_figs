#!/usr/bin/env python
# encoding: utf-8
usage = """
Usage:
It calculates the sum of the range of data relative to your coor

***
<required>:
-s: suffix for all the input files
-c: where the coordinate of your reference point (ex. if this dataset is mapped to -500 to 1000 relative to TSS, then -c 500)
-u: upstream limitation
-d: downstream limitation
-o: output file

Example: python calculate_sum_from_columns.py -s .tab -c coor_location -u -150 -d 150 -o outfile 
***

Created by Kuangyu Yen on 2012-04-27.
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
"""

import sys, getopt, os, csv, re, operator
from itertools import ifilter

#
# define functions
#
def getFiles(suffix):
	all_files = os.listdir(os.getcwd())
	return [x for x in all_files if x.endswith(suffix)]

def sortDictValue(dict):
  id_list = sorted(dict.iteritems(), key=operator.itemgetter(0), reverse=False)
  sort_id, corresponding_value = [], []
  
  for row in id_list:
    sort_id.append(row[0])
    corresponding_value.append(row[1])
  
  return sort_id, corresponding_value


def calculate_sum(infile, left, right):
  reader, dict_sum = csv.reader(open(infile, "rU"), delimiter="\t"), {}
  
  for no, row in enumerate(reader):
    if no > 0: dict_sum[row[0]] = sum(map(float, row[left:right]))
    else: pass 
  
  genelist, values = sortDictValue(dict_sum)
  
  return genelist, values


if __name__ == '__main__':
  if len(sys.argv) < 2 or not sys.argv[1].startswith("-"): sys.exit(usage)
  
  # get argument
  optlist, alist = getopt.getopt(sys.argv[1:], "hs:c:u:d:o:")
  for opt in optlist:
    if opt[0] == "-h": sys.exit(usage)
    elif opt[0] == "-s": suffix = opt[1]
    elif opt[0] == "-c": ref_point = int(opt[1])+1
    elif opt[0] == "-u": up_lim = int(opt[1])
    elif opt[0] == "-d": dn_lim = int(opt[1])
    elif opt[0] == "-o": out_file = opt[1]
    
  left, right = ref_point + up_lim, ref_point + dn_lim
  infiles, dict_exp_value, result = getFiles(suffix), {}, []
  
  for infile in infiles:
    genelist, values = calculate_sum(infile, left, right)
    dict_exp_value[infile] = values
    
    print "finish processing %s" %infile
  
  i = 0
  while i < len(genelist):
    tmp = []
    tmp.append(genelist[i])
    for j in range(len(dict_exp_value.keys())):
		tmp.append(dict_exp_value[dict_exp_value.keys()[j]][i])
      
    result.append(tmp)
    i = i + 1

  
  exp_name = dict_exp_value.keys()
  exp_name.insert(0, "gene")
  writer = csv.writer(open(out_file, "w"), delimiter="\t")
  writer.writerow(exp_name)
  for row in result:
    writer.writerow(row)
  
  
  print "done!"
  

    
 

#!/usr/bin/env python
# encoding: utf-8
usage = """
Usage:
It calculates the sum of the range of data relative to your coor

***
<required>:
suffix for all the input files
where the coordinate of your reference point (ex. if this dataset is mapped to -500 to 1000 relative to TSS, then this is 500)
upstream limit
downstream limit
output file

Example: python calculate_sum_from_columns.py .tab 500 -150 150 outfile 
***

Created by Kuangyu Yen on 2012-04-27.
Copyright (c) 2012 __PughLab@PSU__. All rights reserved.
Edited by Lila Rieber on 10-11-2017
"""

import argparse, os, csv, operator

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
  parser = argparse.ArgumentParser()
  parser.add_argument('suffix',
                        help='Suffix of files to process')
  parser.add_argument('ref_point', type=int,
                        help='Reference point')
  parser.add_argument('up_lim', type=int,
                      help='Upper limit')
  parser.add_argument('dn_lim', type=int, 
                      help='Lower limit')
  parser.add_argument('out_file', 
                      help='File to write to')
    
  args = parser.parse_args()   

  left, right = args.ref_point + args.up_lim, args.ref_point + args.dn_lim
  infiles, dict_exp_value, result = getFiles(suffix), {}, []
  
  for infile in infiles:
    genelist, values = calculate_sum(infile, left, right)
    dict_exp_value[infile] = values
    
    print "finish processing {}".format(infile)
  
  i = 0
  while i < len(genelist):
    tmp = []
    tmp.append(genelist[i])
    for j in range(len(dict_exp_value.keys())):
      tmp.append(dict_exp_value[dict_exp_value.keys()[j]][i])
      
    result.append(tmp)
    i += 1

  
  exp_name = dict_exp_value.keys()
  exp_name.insert(0, "gene")
  writer = csv.writer(open(out_file, "w"), delimiter="\t")
  writer.writerow(exp_name)
  for row in result:
    writer.writerow(row)
  
  
  print "done!"

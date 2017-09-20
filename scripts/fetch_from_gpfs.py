import os
import sys

wd = sys.argv[1]
with open("../{}/tab_files_list".format(wd)) as tab_files:
	for line in tab_files:
		tab_file = line.strip()
		id_string = tab_file[0:5]
		try:
			int(id_string)	#test if proper ID
			id_start = int(id_string[0:2])
			if id_start >=1 and id_start <= 51:
				os.system("rsync lur159@lionxg.rcc.psu.edu:/gpfs/cyberstar/pughhpc/archive/2013_Backup_files/2013_Backup_tab_files/Q{}/{}* ../{}/tab_files".format(id_start, tab_file, wd))
			elif id_start <= 64:
				os.system("rsync lur159@lionxg.rcc.psu.edu:/gpfs/cyberstar/pughhpc/archive/2014_Pugh_lab/{}* ../{}/tab_files".format(tab_file, wd))
		except:
			print "Warning: {} cannot be found on GPFS".format(tab_file)
tab_files.close()

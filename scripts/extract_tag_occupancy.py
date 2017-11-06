import sys

cdt_prefix = sys.argv[1]
u = int(sys.argv[2])
d = int(sys.argv[3])

#x = {}
#y = {}
z = {}

cnt = 0
with open("{}_forward.cdt".format(cdt_prefix)) as forward:
	for line in forward:
		data = line.strip().split("\t")
		key = data[0].split(";")
		cnt += 1
		if cnt > 2:
			#x.update({key[0]:data})
	
			cnt2 = 0
			occ = 0.0
			# 1003 is the TSS midpoint
			start = 1003 + u
			end = 1003 + d
			for ele in data:
				cnt2 += 1
				# Polymerase Occupancy will be calculated by specific interval
				# Rpb3 : -100 to 500
				# Ssl2 : -200 to 200
				# Sua7 : -200 to 200
				if cnt2 > 2 and cnt2 >= start and cnt2 <= end:
					occ += float(ele)
			z.update({key[0]:[occ,0.0]})
forward.close()

cnt = 0
with open("{}_reverse.cdt".format(cdt_prefix)) as reverse:
	for line in reverse:
		data = line.strip().split("\t")
		cnt += 1
		if cnt > 2:
			#y.update({data[0]:data})
		
			cnt2 = 0
			occ = 0.0
			for ele in data:
				cnt2 += 1
				if cnt2 > 2:
					occ += float(ele)
			try:
				z[data[0]][1] = occ
			except:
				z.update({data[0]:[0.0,occ]})
reverse.close()

with open("{}_occupancy.txt".format(cdt_prefix), "w") as out:
	for the_key in z:
		if z[the_key][0] != 0.0 and z[the_key][1] != 0.0:
			print "test"
			out.write(",".join((the_key, str(z[the_key][0]), str(z[the_key][1]))))
			out.write("\n")
out.close()

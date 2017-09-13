import sys

fn = sys.argv[1]
u = int(sys.argv[2])
d = int(sys.argv[3])

f1 = open("Pol2/"+fn+".cdt","r")
f2 = open("mRNA/"+fn+".cdt","r")
x = {}
y = {}
z = {}
cnt = 0
for line in f1:
	line = line.strip()
	data = line.split("\t")
	key = data[0].split(";")
	cnt += 1
	if cnt > 2:
		x.update({key[0]:data})
	
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
cnt = 0
for line in f2:
	line = line.strip()
	data = line.split("\t")
	cnt += 1
	if cnt > 2:
		y.update({data[0]:data})
		
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
f2.close()
f1.close()

print "PolII Occupancy:",len(x)
print "mRNA Occupancy:",len(y)
print "Total dataset:",len(z)

f = open(fn+"_occupancy.txt","w")
f.write("Gene name,PolII,mRNA\n")
for the_key in z:
	if z[the_key][0] != 0.0 and z[the_key][1] != 0.0:
		f.write(the_key+","+str(z[the_key][0])+","+str(z[the_key][1])+"\n")
f.close()

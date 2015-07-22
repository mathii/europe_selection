import fileinput

pops=[]
for line in fileinput.input():
	pops.append(line[:-1])

for i in range(len(pops)-1):
	for j in range(i+1, len(pops)):
		print pops[i]+","+pops[j]

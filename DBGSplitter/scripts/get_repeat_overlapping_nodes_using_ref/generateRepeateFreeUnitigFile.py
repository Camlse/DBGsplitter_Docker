import sys

doNotOutputUnitigs=set()
with open(sys.argv[2]) as fin:
	for line in fin:
		line = line.strip()
		doNotOutputUnitigs.add(int(line))

with open(sys.argv[1]) as fin:
	while True:
		header = fin.readline().strip()
		if header=="":
			break
		seq = fin.readline().strip()
		unitig=int(header[1:])
		if unitig not in  doNotOutputUnitigs:
			print header
			print seq

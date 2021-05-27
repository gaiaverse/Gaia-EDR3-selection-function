import numpy as np
import random
import copy 
import os
from subprocess import check_output
global binDict
global M
global lBar
import sys

rootToFiles = "../../SlimmedData/"
M = 20

if len(sys.argv) > 1:
	rootToFiles = str(sys.argv[1])
	M = int(sys.argv[2])

binDict = {};

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])

def size(filename):
	return os.path.getsize(filename)

def getUniqueFiles(rootToFiles):
	list = os.listdir(rootToFiles)

	baseBranch= {}
	for file in list:
		if ".csv" in file:
			root = file.split(".")[0]
			if "_" in root:
				root = file.split("_")[0]
			
			if root in baseBranch.keys():
				baseBranch[root].append(file)
			else:
				baseBranch[root] = [file]
		
	uniqueFiles = [];
	for files in baseBranch.values():
		if len(files) == 1:
			uniqueFiles.append(files[0])
		else:
			for file in files:
				if "_" in file:
					uniqueFiles.append(file)
	
	
	# ~ uniqueFiles = [];
	# ~ for i in range(0,142):
		# ~ uniqueFiles.append( str(i) + ".csv")
	
	return uniqueFiles

def generateBinDict(rootToFiles):
	files = getUniqueFiles(rootToFiles)
	
	for file in files:
		v = size(rootToFiles + file)
		binDict[file] = v

def Load(assignment):
	bins = np.zeros(M)
	if len(assignment) > 0:
		for key, value in assignment.items():
			bins[value] += binDict[key]
	else:
		bins[0] = sum(binDict.values())
	return bins
	
	
def Score(assignment):

	
	
	bins = Load(assignment)
	
	diff = bins/lBar - 1.0

	return np.linalg.norm(diff)
	
	
def sortShiftAllocation():

	keys = list(binDict.keys())
	values = list(binDict.values() )
	
	keys = [x for _, x in sorted(zip(values,keys),reverse = True)]
	values = sorted(values, reverse = True)
	
	
	assign = {}
	coreLoad = np.zeros(M)
	
	sortedID = 0
	coreList = np.arange(0,M)
	baseList = np.arange(0,M)
	corePos = 0
	while sortedID < len(binDict):
		fileName = keys[sortedID]
		core = coreList[0]
		fileSize = values[sortedID]
		
		assign[fileName] = core
		coreLoad[core] = coreLoad[core] + values[sortedID]
		
		sortedID += 1
		# ~ corePos += 1
		# ~ if corePos == M:
			# ~ corePos = 0

		coreList = [x for _,x in sorted(zip(coreLoad,baseList))]	

	
	
	b = Load(assign)
	
	mu = sum(b)/len(b)
	
	minLoad = min(b)
	rootIDX = np.argmin(b)

	if rootIDX != 0:

		for fileName, core in assign.items():
			if core == 0:
				assign[fileName] = rootIDX
			if core == rootIDX:
				assign[fileName] = 0

	b = Load(assign)
	print("\n\tCore load:")
	for i in range(0,len(b)):
		print("\t\tCore%d -> %d  (%.3f)" % (i,b[i],(b[i]-lBar)/lBar))
	
	print("\n\tDifference from perfect load = %.4f\n\n" % (Score(assign)))
	return assign

def nameSubstruct(name):
	stripped = name[0:len(name) - 4]
	
	if "_" in name:
		f = stripped.split("_")
		bin = int(f[0])
	else:
		bin = int(stripped)
	
	# ~ #OVERRIDE BIN ALLOCATION
	
	#bin = 0;
	return [name,bin]

def printAllocation(fileAssign):
	
	coreAssigns = [None] * M
	
	for fileName, core in fileAssign.items():
		if coreAssigns[core] == None:
			coreAssigns[core] = [fileName]
		else:
			coreAssigns[core].append(fileName)
			
	delim = ","
	output = ""
	for j in range(0,M):
		s = str(j)+ delim
		
		if coreAssigns[j] != None:
			i = 0
			for file in coreAssigns[j]:
				d = nameSubstruct(file)		
				
				s = s + d[0] + delim + str(d[1])
				i += 1
				if i < len(coreAssigns[j]):
					s= s + delim
		output = output + s+ "\n"
	f = open("../../ModelInputs/coreAssignments.dat","w")
	f.write(output)
	f.close()
	
generateBinDict(rootToFiles)

	

lBar = sum(binDict.values())/M



assign = sortShiftAllocation()
printAllocation(assign)


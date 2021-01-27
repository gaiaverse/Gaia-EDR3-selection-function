import numpy as np
import random
import copy 
import os
from subprocess import check_output
global binDict
global M
global lBar
# ~ binDict = {'0.csv': 2168,
 # ~ '1.csv': 77934358,
 # ~ '2.csv': 79709135,
 # ~ #'3.csv': 146956668,
 # ~ '3_a.csv': 4898556,
 # ~ '3_b.csv': 4898556,
 # ~ '3_c.csv': 4898556,
 # ~ '4.csv': 41517690,
 # ~ '5.csv': 46894754,
 # ~ '6.csv': 53163352,
 # ~ '7.csv': 60561416,
 # ~ '8.csv': 69169274,
 # ~ '9.csv': 37907256,
 # ~ '10.csv': 40019183,
 # ~ '11.csv': 42256042,
 # ~ '12.csv': 44534527,
 # ~ '13.csv': 46853127,
 # ~ '14.csv': 49201568,
 # ~ '15.csv': 51614651,
 # ~ '16.csv': 54042223,
 # ~ '17.csv': 56471020,
 # ~ '18.csv': 58869471,
 # ~ '19.csv': 61171633,
 # ~ '20.csv': 63456388,
 # ~ '21.csv': 65835340,
 # ~ '22.csv': 68318737,
 # ~ '23.csv': 70720962,
 # ~ '24.csv': 72765555,
 # ~ '25.csv': 74004982,
 # ~ '26.csv': 73484455,
 # ~ '27.csv': 69381536,
 # ~ '28.csv': 58732195,
 # ~ '29.csv': 39446813,
 # ~ '30.csv': 18588763,
 # ~ '31.csv': 6334581,
 # ~ '32.csv': 2264779,
 # ~ '33.csv': 1385543,
 # ~ '34.csv': 2684287}

rootToFiles = "../../MainData/"


binDict = {};

def wc(filename):
    return int(check_output(["wc", "-l", filename]).split()[0])

def getUniqueFiles(rootToFiles):
	list = os.listdir(rootToFiles)

	baseBranch= {}
	for file in list:
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
	
	return uniqueFiles

def generateBinDict(rootToFiles):
	files = getUniqueFiles(rootToFiles)
	
	for file in files:
		v = wc(rootToFiles + file)
		print(file + " has " + str(v))
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
	
def dumbAllocation():
	assignment= {}
	dictSize = len(binDict)
	assignSave = {}
	nTries = 100000
	minScore = Score(assignment)
	for n in range(0,nTries):
		
		
		
		for key in binDict.keys():

			assignment[key] = random.randint(0,M-1)
		
		score = Score(assignment)
		
	
		if score  < minScore:
	
			print("New best at %d%% = %.3f" % (int(float(n)/nTries*100),float(score)))
			minScore = score
			assignSave = copy.deepcopy(assignment)
	
	print("Optimal Random Assignment:")
	print(assignSave)
	
	b = Load(assignSave)
	mu = sum(b)/len(b)
	print("\nCore load:")
	for i in range(0,len(b)):
		print("Core%d -> %d  (%.3f)" % (i,b[i],(b[i]-lBar)/lBar))
	
	print("\nDifference from perfect load = %.4f" % (Score(assignSave)))
	minLoad = min(b)
	rootIDX = np.argmin(b)
	print("\nRoot should be assigned as core %d" % rootIDX)
	
def sortShiftAllocation():
	sortedDict = dict(sorted(binDict.items(), key=lambda item: item[1], reverse=True))
	keys = list(sortedDict.keys())
	values = list(sortedDict.values() )
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

	
	print("Optimal Smart Assignment:")
	print(assign)
	
	b = Load(assign)
	print(b)
	mu = sum(b)/len(b)
	print("\nCore load:")
	for i in range(0,len(b)):
		print("Core%d -> %d  (%.3f)" % (i,b[i],(b[i]-lBar)/lBar))
	
	print("\nDifference from perfect load = %.4f" % (Score(assign)))
	minLoad = min(b)
	rootIDX = np.argmin(b)
	print("\nRoot should be assigned as core %d" % rootIDX)
	if rootIDX != 0:
		print(assign)
		print("\nReassigning root to 0....")
		for fileName, core in assign.items():
			if core == 0:
				assign[fileName] = rootIDX
			if core == rootIDX:
				assign[fileName] = 0
		print(assign)
	return assign

def nameSubstruct(name):
	stripped = name[0:len(name) - 4]
	
	if "_" in name:
		f = stripped.split("_")
		bin = int(f[0])
	else:
		bin = int(stripped)
	
	return [name,bin]

def printAllocation(fileAssign):
	
	coreAssigns = [None] * M
	
	for fileName, core in fileAssign.items():
		if coreAssigns[core] == None:
			coreAssigns[core] = [fileName]
		else:
			coreAssigns[core].append(fileName)
			
	delim = ",\t"
	output = ""
	for j in range(0,M):
		s = str(j)+ delim
		
		if coreAssigns[j] != None:
			for file in coreAssigns[j]:
				d = nameSubstruct(file)
				s = s + d[0] + delim + str(d[1]) + delim
		output = output + s+ "\n"
	print(output)
	f = open("coreAssignments.dat","w")
	f.write(output)
	f.close()
	
generateBinDict(rootToFiles)
print(binDict)

	
M = 20
lBar = sum(binDict.values())/M



assign = sortShiftAllocation()
printAllocation(assign)


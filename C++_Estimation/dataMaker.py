import random
import math


a = 3
b = 2
c = 10
sigmaFrac = 0.2

nFile = 10
nDataPerFile = 1500


xRange = (-10,10)
yRange = (-10,10)
for fileID in range(0,nFile):
	
	
	name = "Data/MockData_" + str(fileID) + ".dat"
	
	f = open(name,'w')
	
	for n in range(0,nDataPerFile):
		
		x = random.uniform(xRange[0],xRange[1])
		y = random.uniform(yRange[0],yRange[1])
		
		zTrue = a * x + b*y + c
		sigma = max(abs(zTrue) * sigmaFrac,0.1)
		
		zMeasure = random.gauss(zTrue,sigma)
		
		line = "%.3f, %.3f, %.3f, %.3f\n" % (x,y,zMeasure,sigma)
		f.write(line)
		
	f.close()



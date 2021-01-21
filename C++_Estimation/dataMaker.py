import random
import math


a = 3
b = 2
c = -10
d = 7500
sigmaFrac = 100

nFile = 20
nDataPerFile = 800000000

box = 100000
xRange = (-box,box)
yRange = (-box,box)
wRange = (-box,box)

for fileID in range(0,nFile):
        
        
        name = "Data/MockData_" + str(fileID) + ".dat"
        print(name)
        f = open(name,'w')
                
        prevFrac = 0
        fracGap = 0.5
        for n in range(0,nDataPerFile):
                
                x = random.uniform(xRange[0],xRange[1])
                y = random.uniform(yRange[0],yRange[1])
                w = random.uniform(wRange[0],wRange[1])
                zTrue = a * x + b*y + c*w + d
                sigma = max(abs(zTrue) * sigmaFrac,0.1)
                
                zMeasure = random.gauss(zTrue,sigma)
                
                line = "%.3f, %.3f, %.3f, %.3f, %.3f\n" % (x,y,w,zMeasure,sigma)
                f.write(line)
                
                progress = round(100*float(n)/nDataPerFile,1)
                if progress >= prevFrac + fracGap:
                    prevFrac = progress
                    print("\t%.1f %%" % progress)
        f.close()

        


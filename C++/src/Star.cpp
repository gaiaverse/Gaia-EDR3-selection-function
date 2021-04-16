#include "Star.h"
Star::Star()
{
	//designed to throw an error when it uses an unitialized star to access memory
	gBin = -1;
}
Star::Star(std::vector<std::string> data, int bin)
{
	
	//the gaia data is stored as a csv with columns as follows:
	//Column 0: Number of Measurements reported by Gaia
	//Column 1: Number of visits predicted 
	//Column 2->?, the predicted times at which the star was visited (variable number of columns) 
	
	nMeasure = stoi(data[1]);
	nVisit = stoi(data[2]);
	
	//predictions not perfect so perform an adjustment to prevent BAD memory access issues
	if (nMeasure > nVisit)
	{
		nMeasure = nVisit;
	}
	
	TimeSeries = std::vector<unsigned int>(data.size() - 2,0);
	for (int i = 3; i < data.size(); ++i)
	{
		TimeSeries[i -3] = stoi(data[i]);
	}
	
	//the bin # is derived from the file name, so has to be inserted manually
	gBin = bin;
}

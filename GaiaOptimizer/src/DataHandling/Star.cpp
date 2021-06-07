#include "Star.h"
Star::Star()
{
	//designed to throw an error when it uses an unitialized star to access memory
	gBin = -9999999999999;
}
Star::Star(const std::vector<std::string> & data, int bin, const std::vector<int> & gapsStart, const std::vector<int> & gapsEnd)
{
	//the gaia data is stored as a csv with columns as follows:
	//Column 0: Number of Measurements reported by Gaia
	//Column 1: Number of visits predicted 
	//Column 2->?, the predicted times at which the star was visited (variable number of columns) 
	
	int gapsPassed = 0;
	int nGaps = gapsStart.size();
	int nextGapStart = gapsStart[0];
	int nextGapEnd = gapsEnd[0];
	int obsInGaps = 0;
	for (int i = 2; i < data.size(); ++i)
	{
		int t = stoi(data[i]);
		
		while (gapsPassed < nGaps && nextGapEnd < t)
		{
			++gapsPassed;
			nextGapStart = gapsStart[gapsPassed];
			nextGapEnd = gapsEnd[gapsPassed];
		}
		if ( (gapsPassed < nGaps) && ( (t>=nextGapStart) && (t <= nextGapEnd)))
		{
			++obsInGaps;
		}
		TimeSeries.push_back(t);
	}
	
	
	nMeasure = stoi(data[0]);
	nVisit = stoi(data[1]);
	int nEff = nVisit - obsInGaps;
	if (nMeasure > nEff)
	{
		nMeasure = nEff;
	}
	
	if (nMeasure < 5)
	{
		std::cout << "Oh shit - Pipeline problems!" << std::endl;
		exit(-999);
	}

	gBin = bin;
}

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
	
	bool allGapsPassed = false;
	int currentGap = 0;
	int nGaps = gapsStart.size();
	int obs = 0;
	
	for (int i = 2; i < data.size(); ++i)
	{
		int t = stoi(data[i]);
		
		bool inGap = false;
		if (!allGapsPassed)
		{
			if (t > gapsEnd[currentGap])
			{
				while (currentGap < nGaps && t > gapsEnd[currentGap])
				{
					++currentGap;
				}
			}
			
			if (currentGap == nGaps)
			{
				allGapsPassed = true;
			}
			else
			{
				if (t >= gapsStart[currentGap] && t<=gapsEnd[currentGap])
				{
					inGap = true;
				}
			}
		}
		
		if (!inGap)
		{
			++obs;
		}
		else
		{
			TimeSeries.push_back(t);
		}
	}
	
	
	int fileK = stoi(data[0]);
	nVisit = TimeSeries.size();
	
	nMeasure = fileK; //std::min(fileK,obs);
	
	
	if (nMeasure < 5)
	{
		std::cout << "Oh shit - Pipeline problems!" << std::endl;
		exit(-999);
	}

	gBin = bin;
}

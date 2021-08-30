#include "Star.h"

Star::Star(const std::vector<std::string> & data, int bin)
{
	//the gaia data is stored as a csv with columns as follows:
	//Column 0: Number of Measurements reported by Gaia
	//Column 1: Number of visits predicted 
	//Column 2->?, the predicted times at which the star was visited (variable number of columns) 
	
	
	for (int i = 2; i < data.size(); ++i)
	{
		int t = stoi(data[i]);
		
		TimeSeries.push_back(t);
			
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

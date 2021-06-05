#include "Star.h"
Star::Star()
{
	//designed to throw an error when it uses an unitialized star to access memory
	gBin = -9999999999999;
}
Star::Star(const std::vector<std::string> & data, int bin)
{
	//the gaia data is stored as a csv with columns as follows:
	//Column 0: Number of Measurements reported by Gaia
	//Column 1: Number of visits predicted 
	//Column 2->?, the predicted times at which the star was visited (variable number of columns) 
	
	nMeasure = stoi(data[0]);
	nVisit = stoi(data[1]);
	if (nMeasure > nVisit)
	{
		nMeasure = nVisit;
	}
	
	TimeSeries = std::vector<unsigned int>(data.size() - 2,0);
	for (int i = 2; i < data.size(); ++i)
	{
		TimeSeries[i -2] = stoi(data[i]);
	}
	if (TimeSeries.size() != nVisit)
	{
		std::cout << "Something weird happened with a stellar time series!" << std::endl;
		exit(-1000);
	}

	gBin = bin;
}

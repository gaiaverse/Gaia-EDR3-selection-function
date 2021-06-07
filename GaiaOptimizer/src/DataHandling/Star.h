#pragma once
#include <vector>
#include <string>
#include <iostream>
// The star class (once fully integrated) will be used to store line-by-line data from the Gaia data files
// There will be ~1.8bn of these objects, so minimising the memory they occupy will be a high priority
class Star
{
	public:
		unsigned short int nMeasure;
		unsigned short int nVisit;
		unsigned int gBin;
		std::vector<unsigned int> TimeSeries;
		
		Star();
		Star(const std::vector<std::string> &data, int bin);
};

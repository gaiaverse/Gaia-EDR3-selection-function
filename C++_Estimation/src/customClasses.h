#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>




class IndexedFile
{
	public:
		int NLines;
		std::string CurrentLine;
		
		IndexedFile(std::string targetFileName);
		
		std::string GetLine(int n);
	private:
		int lineID;
		std::vector<std::size_t> linePositions;
		std::ifstream FileStream;
		void JumpToLine(int n);
};


class Star
{
	public:
		unsigned short int nMeasure;
		unsigned short int nVisit;
		unsigned int gBin;
		std::vector<unsigned int> TimeSeries;
		
		double x;
		double y;
		double z;
		double err;
		
		Star();
		Star(std::vector<std::string> data);
};




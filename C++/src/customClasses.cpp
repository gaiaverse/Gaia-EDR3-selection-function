#include "customClasses.h"


IndexedFile::IndexedFile(std::string targetFileName)
{
	FileStream.open(targetFileName);
	if (!FileStream.is_open() )
	{
		std::cout << "Could not open file " << targetFileName << " - critical error, program terminating" << std::endl;
		exit(1);
	}
	
	linePositions.resize(0);
	const size_t suitablyLargeNumber = std::numeric_limits<std::streamsize>::max();

	do
	{
		linePositions.push_back(FileStream.tellg());
		
		//up to 3x quicker to ignore the line than to save it!
		FileStream.ignore(suitablyLargeNumber,'\n');
		//std::getline(FileStream, CurrentLine);
	} while (FileStream);
	
	NLines = linePositions.size() - 1; // the final entry is the EOF, so not counted as a line
	FileStream.clear();
	//JumpToLine(0);
	lineID = 0;
	JumpToLine(1);
}

void IndexedFile::JumpToLine(int n)
{
	//if doing sequential access, this can speed up access by 100000s factor!
	if (n != lineID)
	{
		FileStream.seekg(linePositions[n]);
	}
	std::getline(FileStream,CurrentLine);
	lineID =n + 1;
}

std::string IndexedFile::GetLine(int n)
{
	if (n < NLines && n >= 0)
	{
		JumpToLine(n);
		
		return CurrentLine;
	}
	else
	{
		std::cout << "Invalid file line provided -- cannot access line " << n << " - throwing error" << std::endl;
		exit(2);
	}
}


Star::Star()
{
	//~ x = 0;
	//~ y = 0;
	//~ w = 0;
	//~ z = 0;
	//~ err = 0;
}
Star::Star(std::vector<std::string> data)
{
	//~ x = std::stod(data[0]);
	//~ y = std::stod(data[1]);
	//~ w = std::stod(data[2]);
	//~ z = std::stod(data[3]);
	//~ err = std::stod(data[4]);
	
	//how is
	//k = # measure
	//n = # visists
	// t_1....t_n (variable size)
	
	 
	// if k > n, then set k = n 
}

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


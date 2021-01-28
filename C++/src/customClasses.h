#pragma once
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <limits>


//custom classes is a nursery for new classes as I build them and they move out into their own properly maintained files
//anything here is probably temporary or badly written....


//NOT CURRENTLY IN USE
//a class for rapidly indexing the line breaks in large files, allowing for non-linear memory access into large files
//once indexed, it is particularly optimized to jump to a line and read sequentially from that point

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





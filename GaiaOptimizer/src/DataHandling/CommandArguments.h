#pragma once
#include <string>
#include "../Main/GlobalVariables.h"
#include "../libs/JSL/JSL.h"
class CommandArgs
{
	public:
		int RandomSeed;
		bool LoadInStartVector;
		std::string StartVectorLocation;
		double GradLim;
		int MaxSteps;
		
		int FreezeSteps;
		std::string DataSource;
		std::string OutputDirectory;
		
		CommandArgs();
		
		void ReadArguments(int argc, char *argv[],int ProcessRank);
		
	private:
	
};

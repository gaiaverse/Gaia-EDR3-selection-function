#pragma once
#include <string>
#include "../GenericFunctions/FileHandler.h"
#include "../Main/GlobalVariables.h"

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

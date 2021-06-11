#pragma once
#include <string>
#include "../Main/GlobalVariables.h"
#include "../libs/JSL/JSL.h"
using JSL::Argument;
class CommandArgs
{
	public:
		//initialisation stuff
		Argument<int> RandomSeed = Argument<int>(0,"random-seed");
		Argument<std::string> StartVectorLocation = Argument<std::string>("__null_location__","restart");
		
		
		//optimiser properties
		Argument<double> GradLim = Argument<double>(0,"gradlim");
		Argument<int> MaxSteps = Argument<int>(1000,"max-steps");
		Argument<int> FreezeSteps = Argument<int>(0,"burnin");
		
		//save locations
		Argument<std::string> DataSource= Argument<std::string>("../../Data/ShuffledData","data");
		Argument<std::string> OutputDirectory = Argument<std::string>("Output","output");
		
		
		std::vector<JSL::ArgumentInterface *> argPointers = {&RandomSeed, &StartVectorLocation, &GradLim, &MaxSteps, &FreezeSteps, &DataSource, &OutputDirectory};
		CommandArgs();
		
		void ReadArguments(int argc, char *argv[],int ProcessRank);
		
	private:
	
};

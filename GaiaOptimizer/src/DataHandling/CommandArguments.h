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
		Argument<int> SaveSteps = Argument<int>(1,"save-steps");
		Argument<bool> SaveAllSteps = Argument<bool>(false,"unique-temp-save");
		Argument<int> Minibatches = Argument<int>(64,"minibatch");
		
		//save locations
		Argument<std::string> DataSource= Argument<std::string>("../../Data/ShuffledData","data");
		Argument<std::string> OutputDirectory = Argument<std::string>("Output","output");
		
		
		//put pointers in here so can easily loop over the (heterogenous) array
		std::vector<JSL::ArgumentInterface *> argPointers = {&RandomSeed, &StartVectorLocation, &GradLim, &MaxSteps, &FreezeSteps, &DataSource, &OutputDirectory,&SaveSteps,&SaveAllSteps,&Minibatches};
		
		CommandArgs(){};
		
		void ReadArguments(int argc, char *argv[],int ProcessRank)
		{
			//check if using config file or command line arguments
			Argument<std::string> ConfigFile = Argument<std::string>("__null_location__","config",argc,argv);
			Argument<char> ConfigDelimiter = Argument<char>(' ',"config-delim",argc,argv);
			
			if ((std::string)ConfigFile == "__null_location__")
			{
				for (int i = 0; i < argPointers.size(); ++i)
				{
					argPointers[i]->ListParse(argc,argv);
				}
			}
			else
			{
				for (int i = 0; i < argPointers.size(); ++i)
				{
					argPointers[i]->Configure(ConfigFile,ConfigDelimiter);
				}
			}
			
			
			//if necessary, create the directories necessary for output saving
			if (ProcessRank == RootID)
			{
				JSL::mkdirReturn dirReport = JSL::mkdirSafely(OutputDirectory);
				JSL::mkdirReturn dirReport2= JSL::mkdirSafely((std::string)OutputDirectory + "/" + (std::string)TempDirName);
				if ((dirReport.Successful || dirReport2.Successful) == false)
				{
					ERROR(1,"Could not locate or create the output directory " + (std::string)OutputDirectory + " or subdirectories therein.");
				}
			}
	
		}
		
	private:
	
};

#include "CommandArguments.h"

CommandArgs::CommandArgs()
{
	RandomSeed = time(NULL);
	LoadInStartVector = false;
	StartVectorLocation = "";
	GradLim = 1e-2;;
	DataSource = "../../TestSets/magnitudes/";
	OutputDirectory = "Output";	
	MaxSteps = 5000;
}


void CommandArgs::ReadArguments(int argc, char * argv[],int ProcessRank)
{
	bool outDirFlag = false;
	bool seedFlag = false;
	bool gradFlag = false;
	bool targetFlag = false;
	bool startFlag = false;
	bool stepFlag = false;
	for (int i = 1; i < argc; ++i)
	{
		
		std::string arg = argv[i];
		
		if (outDirFlag == true)
		{
			OutputDirectory = arg;
			outDirFlag = false;
		
		}
		if (seedFlag == true)
		{
			RandomSeed = std::stoi(arg);
			if (ProcessRank == RootID)
			{
				GlobalLog(2,
					std::cout << "Root reports random seed set to " << RandomSeed << "\n";
				);
			}
			seedFlag = false;
		}
		
		if (gradFlag == true)
		{
			GradLim = std::stod(arg);
			if (ProcessRank == RootID)
			{
				GlobalLog(2,
					std::cout << "Root reports gradient convergence limit set to " << GradLim << "\n";
				);
			}
			gradFlag = false;
		}
		if (targetFlag == true)
		{
			DataSource = arg;
			if (ProcessRank == RootID)
			{
				GlobalLog(2,
					std::cout << "Root reports data source set to " << DataSource << "\n";
				);
			}
			targetFlag = false;
		}
		if (startFlag == true)
		{
			LoadInStartVector = true;
			StartVectorLocation = arg;
			startFlag = false;
		}
		if (stepFlag == true)
		{
			MaxSteps = stoi(arg);
			stepFlag = false;
		}
		
		if (arg == "-f")
		{
			outDirFlag = true;
		}
		if (arg == "-s")
		{
			seedFlag = true;
		}
		if (arg == "-t")
		{
			targetFlag = true;
		}
		if (arg == "-g")
		{
			gradFlag = true;
		}
		if (arg == "-r")
		{
			startFlag = true;
		}
		if (arg == "-l")
		{
			stepFlag = true;
		}
		if (arg == "-h" || arg == "--help")
		{
			forLineVectorInFile("src/GenericFunctions/commandHelpFile.txt",'|',
				for (int i = 0; i < FILE_LINE_VECTOR.size(); ++i)
				{
					std::cout << std::setw(10) << std::left << FILE_LINE_VECTOR[i];
				}
				std::cout << "\n";
			);
			exit(-1);
			
		}
		
	}
	
	mkdirReturn dirReport = mkdirSafely(OutputDirectory);
	mkdirReturn dirReport2= mkdirSafely(OutputDirectory + "/" + TempDirName);
	if ((dirReport.Successful || dirReport2.Successful) == false)
	{
		ERROR(1,"Could not locate or create the output directory " + OutputDirectory + " or subdirectories therein.");
	}
	
}

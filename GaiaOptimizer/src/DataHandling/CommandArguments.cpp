#include "CommandArguments.h"

CommandArgs::CommandArgs()
{

}


void CommandArgs::ReadArguments(int argc,char * argv[],int ProcessRank)
{
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

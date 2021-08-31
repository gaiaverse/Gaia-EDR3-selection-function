#pragma once
#include <string>
#include "../Main/GlobalVariables.h"
#include "../libs/JSL/JSL.h"
using JSL::Argument;
class CommandArgs
{
	public:
		
		//!The value passed to `srand() <https://www.cplusplus.com/reference/cstdlib/srand/>`_ for reproducible randomness
		Argument<int> RandomSeed = Argument<int>(0,"random-seed");
		
		//!The directory to search in for a valid savefile configuration for relaunch
		Argument<std::string> StartVectorLocation = Argument<std::string>("__null_location__","restart");
		
		
		//!The maximum value of :math:`\nabla\mathcal{L}` which will be considered ‘converged’ via #ADABADAM::StopConditions::GradientThreshold
		Argument<double> GradLim = Argument<double>(0,"gradlim");
		
		//!The maximum number of epochs the optimizer may use before exiting, via #ADABADAM::StopConditions::MaxSteps
		Argument<int> MaxSteps = Argument<int>(1000,"max-steps");
		
		//!The number of steps between saving locations via #ADABADAM::OptimizerProperties::StepsPerPositionSave
		Argument<int> SaveSteps = Argument<int>(1,"save-steps");
		
		//!If true, the temporary vectors are saved uniquely via via #ADABADAM::OptimizerProperties::UniquePositionSaves. Recommended to set this to false to prevent huge amounts of data generation.
		Argument<bool> SaveAllSteps = Argument<bool>(false,"unique-temp-save");
		
		//!The maximum number of batches used per epoch in the SGD prescription via via #ADABADAM::OptimizerProperties::MiniBatches
		Argument<int> Minibatches = Argument<int>(64,"minibatch");
		
		//!The factor by which step sizes are reduced when the harness is active, via #ADABADAM::OptimizerProperties::MaxHarnessFactor
		Argument<double> HarnessSlowDown = Argument<double>(10,"harness-slow");
		
		//!The number of full epochs over which the step size recovers from the HarnessSlowDown, via #ADABADAM::OptimizerProperties::HarnessReleaseSteps
		Argument<int> HarnessRelease = Argument<int>(5,"harness-release");
		
		//!The directory for the stellar data lists
		Argument<std::string> DataSource= Argument<std::string>("../../Data/ShuffledData","data");
		
		//!The director for the output data (created if it doesn’t already exist)
		Argument<std::string> OutputDirectory = Argument<std::string>("Output","output");
		
		
		//! Pointers list so can easily loop over the (heterogenous) array for assigments
		std::vector<JSL::ArgumentInterface *> argPointers = {&RandomSeed, &StartVectorLocation, &GradLim, &MaxSteps, &DataSource, &OutputDirectory,&SaveSteps,&SaveAllSteps,&Minibatches,&HarnessRelease,&HarnessSlowDown};
		
		//!Default constructor....doesn't do anything as the arguments self-initialise
		CommandArgs(){};
		
		
		/*!
		 * Initialise the command-line arguments, and check if a configuration file is requested. Note that there is no checking for repeat arguments or multiply defined trigger strings, so multiple assignments are perfectly possible. 
		 * \param argc The system-provided command argument count
		 * \param *argv[] The system-provided list of command arguments
		 * \param ProcessRank The MPI-provided ID of the current process
		 * \returns Initialises the object against the provided parameters
		 */
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

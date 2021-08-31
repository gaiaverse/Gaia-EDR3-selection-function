#pragma once
#include <ctime>
#include <vector>

namespace ADABADAM
{
	
	//!General properties of the Optimiser - along with StopConditions, these are initialised by the user before optimisation starts.
	struct OptimizerProperties
	{
		/*!The number of \verbatim embed:rst:inline :doc:`minibatches <minibatching>` used by the optimiser \endverbatim. This number is mutable and changes as downsteps are invoked */
		int MiniBatches;
		
		//! The multiplicative prefactor all changes to the position have. In a dumb optimiser it would be the magnitude of the step-vector, but the ADAM stuff messes around with that.
		double StepSize;
		
		//!The number which #MiniBatches is divided by when a downstep is called for
		double MinibatchDownStep;
		
		//! The maximum value of #ProgressTracker::Harness, the rate at which the optimizer slows its progress when a change (such as initialisation, or a minibatch downstep) occurs. Bigger values mean a slower movement. Usually set by the user from the command line via CommandArgs::HarnessSlowDown
		double MaxHarnessFactor;
		
		//! The number of complete epochs over which the Harness disables. The harness changes value every minibatch loop, but will always complete after a set number of epochs. Usually set by the user from the command line via CommandArgs::HarnessRelease
		int HarnessReleaseSteps;
	};
	
	//! A series of conditions and variables which assert when and how the optimiser will quit, and if it will consider itself to have `converge' when it does so.
	struct StopConditions
	{
		//! The maximum number of complete epochs that the optimizer can take. If ProgressTracker::CurrentSteps exceeds this number, the optimizer quits and OptimizerStatus::Converged is set to false.
		int MaxSteps;
		
		//! The smallest value of |dx| that the optimiser can take in an epoch before the OptimizerStatus::Converged is set to true. If set to 0, this condition is ignored.
		double PositionChangeThreshold;
		
		//! The smallest value of |dF/dX| that the Optimizer can take in an epoch before the OptimizerStatus::Converged is set to true. If set to 0, this condition is ignored
		double GradientThreshold;
		
		//! The smallest value of |dF| that the Optimizer can take in an epoch before the OptimizerStatus::Converged is set to true. If set to 0, this condition is ignored
		double FunctionChangeThreshold;
		
		//! After the downstep to 1 minibatch occurs, force at least this number of epochs to occur before the system can terminate (implemented to prevent premature termination after the full information is given to the system). 
		int SingleBatchStepThreshold;
		
		//! If true, looks in the #DownStepFile and #TerminationFile at the end of each epoch, and executes those instructions.
		bool UseExternalInstructions;
		
		//!The location of the file which Optimizer::CheckExternalFiles() searches. If the file contains a 1, triggers a downstep.
		std::string DownStepFile;
		
		//!The location of the file which Optimizer::CheckExternalFiles() searches. If the file contains a 1, triggers the optimisation to complete and OptimizerStatus::Converged is set to false.
		std::string TerminationFile;
		
	};
	
	//! Large scale status variables for the optimizer - 
	struct OptimizerStatus
	{
		bool Continues;
		bool Converged;
		bool TooManySteps;
		bool ReachedGradConvergence;
		bool ReachedStepConvergence;
		bool ReachedFunctionConvergence;
		bool CarryingOnRegardless;
		bool ExternalDownStep;
		bool ExternalTermination;
	};
	struct ProgressTracker
	{
		int CurrentSteps;
		double MovingAverage;
		double Harness;
		std::vector<double> SpeedController;
		int StepsPerPositionSave;
		bool UniquePositionSaves;
		
		bool BufferFileOpened;
		std::string SaveLocation;
		
		int Hashes;
		int MaxHashes;
		int SlowdownTriggers;
		int TimeSinceSingleBatch;
	};
	struct MemoryBuffer
	{
		int Size;
		int Position;
		std::chrono::time_point<std::chrono::system_clock> StartTime;
		std::chrono::time_point<std::chrono::system_clock> LastSaveTime;
		int OverrideTime;
		
		std::vector<double> Fs;
		std::vector<double> Gradnorms;
		std::vector<double> DXs;
		std::vector<double> Times;
		std::vector<int> MiniBatches;
		std::vector<int> Epochs;
		std::vector<int> Batches;
		
		int AnalysisSteps;
		int AnalysisSize;
		std::vector<double> Analysis;
	};
}

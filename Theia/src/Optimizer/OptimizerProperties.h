#pragma once
#include <ctime>
#include <vector>

namespace ADABADAM
{
	
	//!General properties of the Optimiser - along with StopConditions, these are initialised by the user before optimisation starts.
	struct OptimizerProperties
	{
		/*!The number of \verbatim embed:rst:inline :ref:`minibatches <minibatching>` used by the optimizer \endverbatim. This number is mutable and changes as downsteps are invoked */
		int MiniBatches;
		
		//! The multiplicative prefactor all changes to the position have. In a dumb optimizer it would be the magnitude of the step-vector, but the ADAM stuff messes around with that.
		double StepSize;
		
		//!The number which #MiniBatches is divided by when a downstep is called for
		double MinibatchDownStep;
		
		//! The maximum value of #ProgressTracker::Harness, the rate at which the optimizer slows its progress when a change (such as initialisation, or a minibatch downstep) occurs. Bigger values mean a slower movement. Usually set by the user from the command line via CommandArgs::HarnessSlowDown
		double MaxHarnessFactor;
		
		//! The number of complete epochs over which the Harness disables. The harness changes value every minibatch loop, but will always complete after a set number of epochs. Usually set by the user from the command line via CommandArgs::HarnessRelease
		int HarnessReleaseSteps;
		
		//!The number of epochs taken by the optimizer before Optimizer::SaveProgress() is called. Usually set by the user from the command line via CommandArgs::SaveSteps
		int StepsPerPositionSave;
		
		//!If true, the temporary, raw vectors are saved uniquely. Recommended to set this to false to prevent huge amounts of data generation. Usually set by the user from the command line via CommandArgs::SaveAllSteps
		bool UniquePositionSaves;

		//! The beta1 parameter of a normal ADAM optimiser -- the first moment weighting factor
		double adamBeta1;
		
		//! The beta2 parameter of a normal ADAM optimiser -- the second moment weighting factor
		double adamBeta2;
		
		//! The dimensionality of the optimizing space. Determined when Optimizer::Minimize() is called. 
		int Dimensions;
	};
	
	//! A series of conditions and variables which assert when and how the optimizer will quit, and if it will consider itself to have `converge' when it does so.
	struct StopConditions
	{
		//! The maximum number of complete epochs that the optimizer can take. If ProgressTracker::CurrentSteps exceeds this number, the optimizer quits and OptimizerStatus::Converged is set to false.
		int MaxSteps;
		
		//! The smallest value of |dx| that the optimizer can take in an epoch before the OptimizerStatus::Converged is set to true. If set to 0, this condition is ignored.
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
	
	//! Large scale status variables for the optimizer - mostly used to communicate if and why the optimser stopped.
	struct OptimizerStatus
	{
		//! If true, the optimizer continues on to the next loop. If false, quits, unless CarryingOnRegardless gets in the way.
		bool Continues;
		
		//! If true, the optimizer is considered 'converged'. If false and the optimizer stopped, then it was forced to despite still moving. 
		bool Converged;
		
		//! A flag which is set which forces the optimizer to ignore a Continues == false flag
		bool CarryingOnRegardless;
		
		//! Flag which is set when StopConditions::MaxSteps was triggered
		bool TooManySteps;
		
		//! Flag which is set when StopConditions::GradientThreshold was triggered
		bool ReachedGradConvergence;
		
		//! Flag which is set when StopConditions::PositionChangeThreshold was triggered
		bool ReachedStepConvergence;
		
		//! Flag which is set when StopConditions::FunctionChangeThreshold was triggered
		bool ReachedFunctionConvergence;
		
		//! Flag which is set when StopConditions::DownStepFile was triggered
		bool ExternalDownStep;
		
		//! Flag which is set when StopConditions::TerminationFile was triggered
		bool ExternalTermination;
	};
	
	//! Enumerates the current progress and position of the optimizer. Contains information about various trackers and info about how the optimizer is moving.
	struct ProgressTracker
	{
		//! The current number of complete epochs 
		int CurrentSteps;
		
		//! The current number of minibatches, cannot be higher than OptimizerProperties::MiniBatches
		int CurrentMinibatches;
		
		//! The current speed multiplier. Usually set to 1, if an event occurs, the harness is set to 1/#OptimizerProperties::MaxHarnessValue and hence slows down the rate of change of the optimizere, and hence protects it from rapid changes induced by changes in the optimizer, rather than changes in the function.
		double Harness;
		
		//! If false, when Optimizer::SaveProgress() is called, overwrites the existing file and inserts appropriate headers.
		bool BufferFileOpened;
		
		//! The location into which the MemoryBuffer object is saved
		std::string SaveLocation;
		
		//! The current number of hashes in the progress bar
		int Hashes;
		
		//! The number of hashes in a full progress bar
		int MaxHashes;
		
		//! Then number of times that the Harness has been invoked due to i.e. the function increasing in value. 
		int SlowdownTriggers;
		
		//! The number of epochs since a downstep to 1 minibatch. Used to determine if StopConditions::SingleBatchStepThreshold applies
		int EpochsSinceSingleBatch;
		
		//! The current learning rate -- very closely related to OptimizerProperties::StepSize, but can be reduced if a step is taken which increases the functional value. 
		double LearningRate;
		
		//! Previous functional value so per-epoch progress can be compared
		double PreviousEpoch;
		//! Previous functional value so per-step progress can be compared
		double PreviousMinibatch;
		
		//! Used to indicate that the zeroth position has been saved to file
		bool InitialSaveComplete;
		
		//!Stops the harness from slowing it down too much
		int EpochsSinceLastHarness;
	};
	
	//! A set of internal memory buffers and controls for when to execute them. Used to hold a buffer for Optimizer::SaveProgress() and also the long-term analysis in Optimizer::CheckConvergence()
	struct MemoryBuffer
	{
		//! The allowed size of the memory buffer -- when full, saves the file and loops back around
		int Size;
		
		//! The current position of the memory buffer. initialised to 0 and increments each calculation.
		int Position;
		
		//! The time at which the optimizer began the current optimization run
		std::chrono::time_point<std::chrono::system_clock> StartTime;
		
		//! The last time at which a Optimizer::SaveProgress() call was made
		std::chrono::time_point<std::chrono::system_clock> LastSaveTime;
		
		//! If LastSaveTime exceeds this parameter, an Optimizer::SaveProgress() call is made, even if the buffer is not full. Usually set to 5 minutes (in seconds)
		int OverrideTime;
		
		//! A vector containing the last-calculated value of the functional value
		std::vector<double> Fs;
		
		//! A vector containing the last-calculated value of the absolute value of the gradient
		std::vector<double> Gradnorms;
		
		//! A vector containing the last-calculated value of the absolute value |dx|
		std::vector<double> DXs;
		
		//! A vector containing the times (measured in seconds from initialisation) that the other values were calculated at
		std::vector<double> Times;
		
		//! A vector containing the minibatch ID for the associated calculations. -1 is used to indicate a whole-epoch average. 
		std::vector<int> MiniBatches;
		
		//! The epoch IDs for the associated calculations
		std::vector<int> Epochs;
		
		//! The number of minibatches being used for each of the associated calculations
		std::vector<int> Batches;
		
		//! The current position of the AnalysisBuffer
		int AnalysisSteps;
		
		//! The size of the AnalysisBuffer
		int AnalysisSize;
		
		//! A buffer of previously calculated values of F. Stored separately from #Fs as the Analysis buffer is used to determine downstep properties so it is beneficial for its size to be small, indepenedent of the buffer.
		std::vector<double> Analysis;
	};
}

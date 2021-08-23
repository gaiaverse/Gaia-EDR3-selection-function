#pragma once
#include <ctime>
#include <vector>

namespace ADABADAM
{
	struct OptimizerProperties
	{
		int MiniBatches;
		int BurnInSteps;
		double StepSize;
		double MinibatchDownStep;
		double MaxHarnessFactor;
		int HarnessReleaseSteps;
	};
	struct StopConditions
	{
		int MaxSteps;
		double PositionChangeThreshold;
		double GradientThreshold;
		double FunctionChangeThreshold;
		int SingleBatchStepThreshold;
		bool UseExternalInstructions;
		std::string DownStepFile;
		std::string TerminationFile;
		
	};
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

#pragma once
#include <vector>
#include <math.h>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <iostream>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"

#include "OptimizerProperties.h";
using Eigen::VectorXd;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
std::vector<int> randomShuffle(int n)
{
	std::vector<int> order;
	for (int i=0; i< n; ++i) 
	{
		order.push_back(i);
	}
	std::random_shuffle ( order.begin(), order.end() );

	return order;
}

template<class T>
class Optimizer
{
	private:
		T & Functor;

		MemoryBuffer Buffer;
		
	public:
		
		OptimiserStatus Status;
		OptimiserProperties Properties;
		StopConditions HaltConditions;
		ProgressTracker Progress;
		
		Optimizer<T>(T& functor) : Functor(functor)
		{
			SetDefaults();	
		}
		
		void SetDefaults()
		{
			HaltConditions.MaxSteps = 10000;
			HaltConditions.PositionChangeThreshold = 1e-5;
			HaltConditions.GradientThreshold = 1e-6;
			HaltConditions.FunctionChangeThreshold = 1e-7;
			
			Properties.Mode = OptimiserModes::ADABADAM;
			Properties.MiniBatches = 1;
			Properties.BurnInSteps = 0;
			Properties.StepSize = 0.1;			
			
			Properties.HarnessReleaseFactor = 0.25;
			Properties.MaxHarnessFactor = 100;
			
			Properties.MinibatchDownStep = 4;
			
			Buffer.Size = 10;
			Buffer.AnalysisSize = 10;
			Buffer.OverrideTime = 300;
			Progress.SaveLocation = "";
			Progress.MaxHashes = 16;
			Progress.StepsPerPositionSave = 1;
		}
		
		void Initialise()
		{
			Status.TooManySteps = false;
			Status.Converged = false;
			Status.ReachedFunctionConvergence = false;
			Status.ReachedGradConvergence = false;
			Status.ReachedStepConvergence = false;
			Progress.Harness = 1.0/Properties.MaxHarnessFactor;
			Progress.BufferFileOpened = false;
			
			Buffer.Position = 0;
			Buffer.StartTime = std::chrono::system_clock::now();
			Buffer.LastSaveTime = std::chrono::system_clock::now();
	
			int n = Buffer.Size;
			
			Buffer.Batches = std::vector<int>(n,0);
			Buffer.MiniBatches = std::vector<int>(n,0);
			Buffer.Fs = std::vector<double>(n,0);
			Buffer.DFs = std::vector<double>(n,0);
			Buffer.Gradnorms = std::vector<double>(n,0);
			Buffer.Times = std::vector<double>(n,0);
			Buffer.Epochs = std::vector<int>(n,0);
			
			Buffer.AnalysisSteps = 0;
			Buffer.Analysis = std::vector<double>(Buffer.AnalysisSize,0);
			
			Progress.Hashes = 0;
			Progress.CurrentSteps = 0;
			Status.Continues = true;
		}
			
		void Minimize(VectorXd & x)
		{

			Initialise();
			
			switch(Properties.Mode)
			{
				case OptimiserModes::ADABADAM:
				{
					ADABADAM(x);
					break;
				}
			}
		}
		
		void ADABADAM(VectorXd &x)
		{
			int EffectiveBatches = Properties.MiniBatches;
			int Dimensions = x.size();
			//initialise ADAM vectors
			VectorXd m = VectorXd::Zero(Dimensions);
			VectorXd v = VectorXd::Zero(Dimensions);
			VectorXd epochGradient = VectorXd::Zero(Dimensions);
			VectorXd oldX = x;
			
			//ADAM Variables
			double beta1 = 0.9;
			double beta2 = 0.99;
			double eps = 1e-8;
			double learningRate = Properties.StepSize;
			
			double previousEpoch = 99999999;
			double previousMinibatch = 9999999;
			
			int t = 1;

			bool burnInStopped = false;
			while (Status.Continues)
			{
				int epochs = Progress.CurrentSteps + 1;
				double epochL = 0;
				epochGradient -= epochGradient;
				
				if (!burnInStopped && Progress.CurrentSteps >= Properties.BurnInSteps)
				{
					burnInStopped = true;
					Functor.Unfreeze();
				}

				
				std::vector<int> batchOrder = randomShuffle(EffectiveBatches);

		
				for (int batches = 0; batches < EffectiveBatches; ++batches)
				{
					std::cout << Progress.Harness << std::endl;
					int currentBatch = batchOrder[batches];
							
					Functor.Calculate(x,currentBatch,EffectiveBatches);
					
					//save initial position
					if (batches == 0 && epochs == 1)
					{
						Functor.SavePosition(false,0,Progress.UniquePositionSaves,x);
						previousEpoch = Functor.Value;
					}
					
					double b1Mod = 1.0/(1.0 - pow(beta1,t));
					double b2Mod = 1.0/(1.0 - pow(beta2,t));			
					
					double gNorm = 0;
					double dxNorm = 0;
					for (int i = 0; i < Dimensions; ++i)
					{
						double g = Functor.Gradient[i];
						gNorm += g*g;
						m[i] = beta1 * m[i] + (1.0 - beta1)*g;
						v[i] = beta2 * v[i] + (1.0 - beta2) * (g*g);
						
						double dx_i = -b1Mod * m[i] /  ((sqrt(v[i]*b2Mod) + eps) ) * learningRate * Progress.Harness;
						dxNorm += dx_i * dx_i;
						x[i] += dx_i;
						epochGradient[i] += g;
					}	
					++t;
					
					epochL += Functor.Value;

					if (EffectiveBatches > 1)
					{
						double df_mini = Functor.Value - previousMinibatch;
						previousMinibatch = Functor.Value;
						double sqrtgNorm = sqrt(gNorm);
						UpdateProgress(batches,EffectiveBatches,Functor.Value,sqrtgNorm,df_mini,sqrt(dxNorm),x);
					}
					Progress.Harness = std::min(1.0,Progress.Harness * (1.0 + Properties.HarnessReleaseFactor));
				}
				
				epochL/=EffectiveBatches;
				epochGradient *= 1.0/EffectiveBatches;
				double df = epochL - previousEpoch;
				previousEpoch = epochL;
				double epochGradNorm = epochGradient.norm();
				double epochDx = (x-oldX).norm();
				oldX = x;
				
				++Progress.CurrentSteps;
				UpdateProgress(-1,EffectiveBatches,epochL,epochGradNorm,df,epochDx,x);
				
				double newBatches = CheckMinibatches(df,EffectiveBatches);
				
				if (newBatches < EffectiveBatches)
				{
					EffectiveBatches = newBatches;
					Progress.Harness = 1.0/Properties.MaxHarnessFactor;
					learningRate = std::min(learningRate*1.1,2*Properties.StepSize);
					std::cout << "\t\t\t\tThe stepsize has been reduced to " << EffectiveBatches << " with a learning rate " << learningRate << std::endl;
				}
				
				
				CheckConvergence(epochGradient,df,epochDx);
				if (Status.Continues == false && EffectiveBatches > 1)
				{
					Status.Continues = true;
					EffectiveBatches = std::max(1,EffectiveBatches/4);
					Progress.Harness = 1.0/Properties.MaxHarnessFactor;
					std::cout << "\t\t\t\tThe stepsize has been reduced to " << EffectiveBatches << " with a learning rate " << learningRate << std::endl;
				}
				
				
			}
			SaveProgress(Buffer.Position);
			Functor.SavePosition(true,0,Progress.UniquePositionSaves,x);
		}
			
		void CheckConvergence(const VectorXd & dg, double df, double dx)
		{
			if (Progress.CurrentSteps > HaltConditions.MaxSteps)
			{
				std::cout << "STeps! " << Progress.CurrentSteps << std::endl;
				Status.TooManySteps = true;
				Status.Continues = false;
			}
			if (HaltConditions.PositionChangeThreshold > 0 && dx < HaltConditions.PositionChangeThreshold)
			{
				std::cout << "Position! " << dx << std::endl;
				Status.ReachedStepConvergence = true;
				Status.Converged = true;
				Status.Continues = false;
			}
			if (HaltConditions.GradientThreshold > 0 && dg.norm() < HaltConditions.GradientThreshold)
			{
				std::cout << "Gradient! " << dg.norm() << std::endl;
				Status.ReachedGradConvergence = true;
				Status.Converged = true;
				Status.Continues = false;
			}
			if (HaltConditions.FunctionChangeThreshold > 0  && abs(df) < HaltConditions.FunctionChangeThreshold)
			{
				std::cout << "Function! " << df << std::endl;
				Status.ReachedFunctionConvergence = true;
				Status.Converged = true;
				Status.Continues = false;
			}
			if (Progress.MovingAverage > 0)
			{
				Status.ReachedFunctionConvergence = true;
				Status.Converged = true;
				Status.Continues = false;
			}
		} 

		int CheckMinibatches(double df,int currentSize)
		{
		
			int analysisPos = Buffer.AnalysisSteps % Buffer.AnalysisSize; 
			
			Buffer.Analysis[analysisPos] = df;
			++Buffer.AnalysisSteps;	

			double newSize = currentSize;
			
			if (Buffer.AnalysisSteps >= Buffer.AnalysisSize)
			{
				if (NeedsBatchReduction())
				{	
					Buffer.AnalysisSteps = 0;
					newSize = currentSize / Properties.MinibatchDownStep;
					if (newSize < 1)
					{
						newSize = 1;
					}
				}
			
			}		
			return newSize;
		}

		bool NeedsBatchReduction()
		{
			bool batchesAreAProblem = false;
			
			int N = Buffer.AnalysisSize;
			
			double sum = 0;
			int signChanges = 0;
	
			for (int i = 0; i <N; ++i)
			{
				double x = Buffer.Analysis[i];

				sum += x;
				if (i > 0 && (sgn(x) != sgn(Buffer.Analysis[i-1])) &&  abs(x) > 0)
				{
					++signChanges;
				} 
			}

			double mean = sum / N;
		
			int problematicSignChanges = std::max(2,N/3);
			if ((mean > 0) || signChanges >= problematicSignChanges)
			{

				batchesAreAProblem = true;
			}


			return batchesAreAProblem;
		}

		void UpdateProgress(int batch, int nBatches,double F, double G, double dF,double dxNorm,const VectorXd & x)
		{
			int i = Buffer.Position;
			Buffer.MiniBatches[i] = batch;
			Buffer.Batches[i] = nBatches;
			Buffer.Fs[i] = F;
			Buffer.Gradnorms[i] = G;
			Buffer.DFs[i] = dF;
			Buffer.Epochs[i] = Progress.CurrentSteps + 1;
			auto time = std::chrono::system_clock::now();
			std::chrono::duration<double> diff = time - Buffer.StartTime;
			
			Buffer.Times[i] = diff.count();
			
			++Buffer.Position;
			
			std::chrono::duration<double> savediff = time - Buffer.LastSaveTime;
			if (Buffer.Position >= Buffer.Size || savediff.count() >= Buffer.OverrideTime)
			{
				SaveProgress(std::min(Buffer.Size,Buffer.Position));
				Buffer.Position = 0;
				Buffer.LastSaveTime = time;
			}
			
			if (batch == -1)
			{
				double frac = 0.9;
				Progress.MovingAverage = frac*Progress.MovingAverage + (1.0- frac)*dF;
				
				
				if (nBatches == 1)
				{
					std::cout << "\t\tEpoch " << Progress.CurrentSteps;
				}
				else
				{
					std::cout << "]";
				}
				std::cout << " complete, at Calculation Evaluation " << Functor.LoopID << "\n";
				std::cout << "\t\t\t(L,Gradnorm,dL,|dx|,nBatch) = (" << std::setprecision(10) << F << ", " <<  std::setprecision(10) << G << ", " << std::setprecision(10) << dF << ", " << std::setprecision(10) << dxNorm << ", " << nBatches <<")\n"; 
				
				
				
				if (Progress.CurrentSteps % Progress.StepsPerPositionSave == 0)
				{
					Functor.SavePosition(false,Progress.CurrentSteps,Progress.UniquePositionSaves,x);
				}
				std::cout << "\t\t\t" << JSL::PrintCurrentTime();
			}
			else
			{
				if (batch == 0)
				{
					std::cout << "\t\tEpoch " << Progress.CurrentSteps + 1 << "   [";
					Progress.Hashes = 0;
				}
				int nHashes = (batch+1) * Progress.MaxHashes / nBatches;
				//~ std::cout << batch << "   "<<nHashes << std::endl;
				if (nHashes > Progress.Hashes)
				{
					std::string h = "";
					while (Progress.Hashes < nHashes)
					{
						h+="#";
						++Progress.Hashes;
					}
					std::cout << h << std::flush;
				}
			}
			
			
			
		}

		void SaveProgress(int n)
		{
			std::string saveFile = Progress.SaveLocation + "OptimiserProgress.txt";
			std::fstream file;
			int width = 20;
			if (Progress.BufferFileOpened == false)
			{
				file.open(saveFile,std::ios::out);
				std::vector<std::string> headers = {"Elapsed","Epoch","Batch", "nBatches","F","dF","GradNorm"};
				for (int i = 0; i < headers.size(); ++i)
				{
					file << std::setw(width) << headers[i] + ",";
				}
				file << "\n";
				Progress.BufferFileOpened = true;
			}
			else
			{
				file.open(saveFile,std::ios::app);
			}
			
			for (int i = 0; i < n; ++i)
			{
				int prec = 8l;
				file << std::setw(width) << std::setprecision(prec) << Buffer.Times[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Buffer.Epochs[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Buffer.MiniBatches[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Buffer.Batches[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Buffer.Fs[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Buffer.DFs[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Buffer.Gradnorms[i] << ",";
				
				file << "\n";
			}
			
			file.close();
		}

		std::string GetStatus()
		{
			std::string s = "";
			s += "Steps Taken: " + std::to_string(Progress.CurrentSteps) + " / " + std::to_string(HaltConditions.MaxSteps);
			s += "\nHalt conditions: ";
			std::vector<std::string> titles = {"Too many steps", "Reached Gradient Convergence", "Reached Step Convergence", "Reached Functional Convergence"};
			std::vector<bool> values = {Status.TooManySteps, Status.ReachedGradConvergence, Status.ReachedStepConvergence, Status.ReachedFunctionConvergence};
			for (int i = 0; i < titles.size(); ++i)
			{
				s +=  "\n\t" + titles[i] + ": ";
				std::string answer = "no";
				if (values[i])
				{
					answer = "yes";
				} 
				s += answer;
			}
			s+= "\n";
			return s;
		}
		
};



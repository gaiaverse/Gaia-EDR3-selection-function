#pragma once
#include <vector>
#include <math.h>
#include <ctime>
#include <iomanip>
#include <fstream>
#include "../GenericFunctions/timeCodes.h"
#include <iostream>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"
using Eigen::VectorXd;

template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
struct Conditions
{
	int MaxSteps;
	double StepSize;
	double xConvergence;
	double gConvergence;
	double fConvergence;
	int SaveSteps;
	int InitialStepMemory;
};
struct Statuses
{
	

	int CurrentSteps;
	
	bool TooManySteps;
	bool ReachedGradConvergence;
	bool ReachedStepConvergence;
	bool ReachedFunctionConvergence;
};
struct Progresser
{
	std::string ProgressDir;
	int BufferSize;
	bool HasSaved;
	int BufferPosition;
	std::chrono::time_point<std::chrono::system_clock> StartTime;

	std::chrono::time_point<std::chrono::system_clock> LastBufferTime;
	int BufferTime;

	std::vector<double> PastFs;
	std::vector<double> PastGradNorms;
	std::vector<double> PastDFs;
	std::vector<double> PastTimes;
	std::vector<int> PastMiniBatch;
	std::vector<int> PastEpoch;
	std::vector<int> PastBatchCount;
	
	int AnalysisSteps;
	int AnalysisMemorySize;
	std::vector<double> AnalysisMemory;
	double AnalysisMemoryGrad;
	int Hashes;
	int MaxHashes;
};

template<class T>
class Optimizer
{
	private:
		T & Functor;
		int Dimensions;
		
	public:
		bool Converged;
		Conditions Condition;
		Statuses Status;
		Progresser Progress;
		
		Optimizer<T>(int nRawParams, T& functor) : Functor(functor)
		
		{
			Dimensions = nRawParams;
			SetDefaults();	
		}
		
		
		
		void SetDefaults()
		{
			Converged = false;
			Condition.MaxSteps = 10000;
			Condition.xConvergence = 0;
			Condition.gConvergence = 1e-7;
			Condition.fConvergence = 0;
			Condition.StepSize = 0.02;
			Condition.SaveSteps = 5;
			Condition.InitialStepMemory= 10;
			Status.CurrentSteps = 0;
		
			Status.TooManySteps = false;
			Status.ReachedGradConvergence = false;
			Status.ReachedStepConvergence = false;
			Status.ReachedFunctionConvergence = false;
			
			Progress.BufferSize = 10;
			Progress.ProgressDir = "";
			Progress.AnalysisMemorySize = 5;
			Progress.MaxHashes = 20;
			
		}
		
		void InitialiseProgress()
		{
			Progress.HasSaved = false;
			Progress.BufferPosition = 0;
			Progress.StartTime = std::chrono::system_clock::now();
			Progress.LastBufferTime = std::chrono::system_clock::now();
			Progress.BufferTime = 300;
			int n = Progress.BufferSize;
			
			Progress.PastBatchCount = std::vector<int>(n,0);
			Progress.PastMiniBatch = std::vector<int>(n,0);
			Progress.PastFs = std::vector<double>(n,0);
			Progress.PastDFs = std::vector<double>(n,0);
			Progress.PastGradNorms = std::vector<double>(n,0);
			Progress.PastTimes = std::vector<double>(n,0);
			Progress.PastEpoch = std::vector<int>(n,0);
			Progress.AnalysisSteps = 0;
			Progress.AnalysisMemory = std::vector<double>(Progress.AnalysisMemorySize,0);
			Progress.AnalysisMemoryGrad = 0;
		}
		
		void Minimize(VectorXd & x)
		{
			if (x.size() != Dimensions)
			{
				std::cout << "OPTIMIZER ERROR: Initial position vector is not of the provided size." << std::endl;
				exit(2);
			}
			
			ADAM(x);
		}
		
		void Minimize(VectorXd & x, int nBatches)
		{
			if (x.size() != Dimensions)
			{
				std::cout << "OPTIMIZER ERROR: Initial position vector is not of the provided size." << std::endl;
				exit(2);
			}
			
			InitialiseProgress();
			ADABADAM(x,nBatches);
			
		}
		
	
		void ADABADAM(VectorXd &x,int nBatches)
		{
			int EffectiveBatches = nBatches;
			
			//initialise ADAM vectors
			VectorXd m = VectorXd::Zero(Dimensions);
			VectorXd v = VectorXd::Zero(Dimensions);
			VectorXd epochGradient = VectorXd::Zero(Dimensions);
			//~ VectorXd dx = VectorXd::Zero(Dimensions);
			
			//ADAM Variables
			double beta1 = 0.9;
			double beta2 = 0.999;
			double eps = 1e-8;
	
			double previousEpoch = 99999999;
			double previousMinibatch = 9999999;
			
			int epochs; 
			int t = 1;
			bool minimiseContinues = true;
			while (minimiseContinues)
			{
				epochs = Status.CurrentSteps + 1;
				double epochL = 0;
				epochGradient -= epochGradient;
				
				//prepare the pseudo-random batches
				std::vector<int> batchOrder;
				for (int i=0; i<EffectiveBatches; ++i) batchOrder.push_back(i);
				std::random_shuffle ( batchOrder.begin(), batchOrder.end() );
				
				
				for (int batches = 0; batches < EffectiveBatches; ++batches)
				{
					int currentBatch = batchOrder[batches];
							
					Functor.Calculate(x,currentBatch,EffectiveBatches);
					
					double b1Mod = 1.0/(1.0 - pow(beta1,t));
					double b2Mod = 1.0/(1.0 - pow(beta2,t));			
					
					//~ VectorXd grad = Eigen::Map<Eigen::VectorXd>(Functor.Gradient.data(),Functor.Gradient.size());
					//~ m = (  beta1 *m + (1.0-beta1)*grad);
					//~ v = ( beta2 * v + (1.0-beta2)* (VectorXd)( grad.array() * grad.array()) );
				
					//~ dx = b1Mod * m * Condition.StepSize;
	
					double gNorm = 0;
					for (int i = 0; i < Dimensions; ++i)
					{
						double g = Functor.Gradient[i];
						gNorm += g*g;
						m[i] = beta1 * m[i] + (1.0 - beta1)*g;
						v[i] = beta2 * v[i] + (1.0 - beta2) * (g*g);
						
						double dx_i = b1Mod * m[i] * Condition.StepSize /  (sqrt(v[i]*b2Mod) + eps);
						x[i] -= dx_i;
						epochGradient[i] += g;
					}	
					++t;
					
					epochL += Functor.Value;
					//epochGradient += grad;
					if (EffectiveBatches > 1)
					{
						double df_mini = Functor.Value - previousMinibatch;
						previousMinibatch = Functor.Value;
						double sqrtgNorm = sqrt(gNorm);
						UpdateProgress(batches,EffectiveBatches,Functor.Value,sqrtgNorm,df_mini);
						Progress.AnalysisMemoryGrad += sqrtgNorm;
					}
				}
				
				epochL/=EffectiveBatches;
				epochGradient *= 1.0/EffectiveBatches;
				double df = epochL - previousEpoch;
				previousEpoch = epochL;
				double epochGradNorm = epochGradient.norm();
				
				++Status.CurrentSteps;
				UpdateProgress(-1,EffectiveBatches,epochL,epochGradNorm,df);
				
				EffectiveBatches = UpdateBatchSize(df,EffectiveBatches,epochGradNorm);
				minimiseContinues = CheckContinues(epochGradient,df);
				
				if (minimiseContinues == false && EffectiveBatches > 1)
				{
					minimiseContinues = true;
					EffectiveBatches = 1;
				}
				
				
			}
			SaveProgress(Progress.BufferPosition);
			Functor.SavePosition(true,0);
		}
		
		void GradientTester(VectorXd &x)
		{
			double dx = 1e-6;
			
			
			std::cout << "Performing a gradient test at the following location: \n\t x = " << x.transpose() << std::endl;
			Functor.Calculate(x);
			double F = Functor.Value;
			std::vector<double> Grad = Functor.Gradient;
			for (int i = Dimensions - 1; i >= 0;--i)
			{
				VectorXd xH = x;
				xH[i] += dx;
				Functor.Calculate(xH);
				
				double measureGrad = (Functor.Value - F)/dx;
				double err = (measureGrad - Grad[i])/Grad[i];
				std::cout << "dF/dx_" <<  std::setprecision(10) << i << " = " << measureGrad << "   (theory value = " <<  std::setprecision(10) << Grad[i] << ", err = " << err << ")\n";
				
			}
			
		}
			
		bool CheckContinues(const VectorXd & dg, double df)
		{
		
			if (Status.CurrentSteps > Condition.MaxSteps)
			{
				//~ std::cout << "STEPS" << std::endl;
				Status.TooManySteps = true;
				return false;
			}
			//~ if (Condition.xConvergence > 0 && dx.norm() < Condition.xConvergence)
			//~ {
				//~ std::cout << "X " << dx.norm() << std::endl;
				//~ Status.ReachedStepConvergence = true;
				//~ Converged = true;
				//~ return false;
			//~ }
			if (Condition.gConvergence > 0 && dg.norm() < Condition.gConvergence)
			{
				//~ std::cout << "G " << dg.norm() << std::endl;
				Status.ReachedGradConvergence = true;
				Converged = true;
				return false;
			}
			if (Condition.fConvergence > 0 && abs(df) < Condition.fConvergence)
			{
				//~ std::cout << "DF " << abs(df) << std::endl;
				Status.ReachedFunctionConvergence = true;
				Converged = true;
				return false;
			}

			
			return true;
		} 

		int UpdateBatchSize(double df,int currentSize,double meanGNorm)
		{
		
			int analysisPos = Progress.AnalysisSteps % Progress.AnalysisMemorySize; 
			
			Progress.AnalysisMemory[analysisPos] = df;
			++Progress.AnalysisSteps;	
			Progress.AnalysisMemoryGrad /= currentSize;
			double newSize = currentSize;
			
			if (Progress.AnalysisSteps >= Progress.AnalysisMemorySize)
			{
			
				if (NeedsBatchReduction(meanGNorm))
				{
				
					Progress.AnalysisSteps = 0;
					newSize = currentSize / 2;
					if (newSize < 1)
					{
						newSize = 1;
					}
				}
			
			}		
			Progress.AnalysisMemoryGrad = 0;
			return newSize;
		}

		bool NeedsBatchReduction(double meanGNorm)
		{
			bool batchesAreAProblem = false;
			
			int N = Progress.AnalysisMemorySize;
			
			double sum = 0;
			int signChanges = 0;
	
			for (int i = 0; i <N; ++i)
			{
				double x = Progress.AnalysisMemory[i];

				sum += x;
				if (i > 0 && (sgn(x) != sgn(Progress.AnalysisMemory[i-1])) &&  abs(x) > 0)
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

			double gradThreshold = 3;

			if (3*meanGNorm < Progress.AnalysisMemoryGrad)
			{
				batchesAreAProblem = true;
			}

			return batchesAreAProblem;
		}

		void UpdateProgress(int batch, int nBatches,double F, double G, double dF)
		{
			int i = Progress.BufferPosition;
			Progress.PastMiniBatch[i] = batch;
			Progress.PastBatchCount[i] = nBatches;
			Progress.PastFs[i] = F;
			Progress.PastGradNorms[i] = G;
			Progress.PastDFs[i] = dF;
			Progress.PastEpoch[i] = Status.CurrentSteps + 1;
			auto time = std::chrono::system_clock::now();
			std::chrono::duration<double> diff = time - Progress.StartTime;
			
			Progress.PastTimes[i] = diff.count();
			
			++Progress.BufferPosition;
			
			std::chrono::duration<double> savediff = time - Progress.LastBufferTime;
			if (Progress.BufferPosition >= Progress.BufferSize || savediff.count() >= Progress.BufferTime)
			{
				SaveProgress(std::min(Progress.BufferSize,Progress.BufferPosition));
				Progress.BufferPosition = 0;
				Progress.LastBufferTime = time;
			}
			
			if (batch == -1)
			{
				if (nBatches == 1)
				{
					std::cout << "\t\tEpoch " << Status.CurrentSteps;
				}
				else
				{
					std::cout << "]";
				}
				std::cout << " complete, at Calculation Evaluation " << Functor.LoopID << "\n";
				std::cout << "\t\t\t(L,Gradnorm,dL,nBatch) = (" << std::setprecision(10) << F << ", " <<  std::setprecision(10) << G << ", " << std::setprecision(10) << dF << ", " << nBatches <<")\n"; 
				
				
				
				if (Status.CurrentSteps % Condition.SaveSteps == 0)
				{
					Functor.SavePosition(false,Status.CurrentSteps);
				}
				std::cout << "\t\t\t"; printTime();
			}
			else
			{
				if (batch == 0)
				{
					std::cout << "\t\tEpoch " << Status.CurrentSteps + 1 << "   [";
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
			std::string saveFile = Progress.ProgressDir + "OptimiserProgress.txt";
			std::fstream file;
			int width = 20;
			if (Progress.HasSaved == false)
			{
				file.open(saveFile,std::ios::out);
				std::vector<std::string> headers = {"Elapsed","Epoch","Batch", "nBatches","F","dF","GradNorm"};
				for (int i = 0; i < headers.size(); ++i)
				{
					file << std::setw(width) << headers[i] + ",";
				}
				file << "\n";
				Progress.HasSaved = true;
			}
			else
			{
				file.open(saveFile,std::ios::app);
			}
			
			for (int i = 0; i < n; ++i)
			{
				int prec = 8l;
				file << std::setw(width) << std::setprecision(prec) << Progress.PastTimes[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Progress.PastEpoch[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Progress.PastMiniBatch[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Progress.PastBatchCount[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Progress.PastFs[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Progress.PastDFs[i] << ",";
				file << std::setw(width) << std::setprecision(prec) << Progress.PastGradNorms[i] << ",";
				
				file << "\n";
			}
			
			file.close();
		}

		std::string GetStatus()
		{
			std::string s = "";
			s += "Steps Taken: " + std::to_string(Status.CurrentSteps) + " / " + std::to_string(Condition.MaxSteps);
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
		
		
				
		void BinaryWolfeDescent(VectorXd & x)
		{
			double prevF = 0;
			Functor.Calculate(x);
			bool minimiseContinues = true;
			double alpha;
			double alphaInit = 1;
			
			double c1_orig = 1e-2;
			double c2_orig = 0.9;
			

			while (minimiseContinues)
			{
				double c1 = c1_orig;
				double c2 = c2_orig;
			
				double OriginalValue = Functor.Value;
				VectorXd Grad = Functor.Gradient;
				double alpha = alphaInit;
				bool alphaFound = false;
				double gNorm = Grad.norm();
				VectorXd pk = -1*Grad; 
				double armijoValue = pk.dot(Grad);
				VectorXd dx;
				int alphaSteps = 0;
				
				
				while (!alphaFound)
				{
					dx = alpha * pk;
					VectorXd xHyp = x + dx;
					Functor.Calculate(xHyp);
					bool armijoSuccess = (Functor.Value <= OriginalValue);
					//~ bool curvatureSuccess = ( - pk.dot(Functor.Gradient) <= - c2* armijoValue);
					bool nanSuccess = ! (std::isnan(Functor.Value) || Functor.Gradient.hasNaN() );
					
					std::cout << "\t\t\tTrying alpha = " << alpha << " which gives |dx| = " << dx.norm() << " \n\t\t\t\tL = " << Functor.Value << " <=! " << OriginalValue - alpha*c1*armijoValue;
					std::cout << "\n\t\t\t\tGrad: " <<   Grad.dot(Functor.Gradient) << "<=! " << c2*armijoValue << "\n";
					
					if (armijoSuccess && nanSuccess)
					{
						alphaFound = true;
						x = xHyp;
					}
					else
					{
						//~ std::cout << Grad.transpose() << std::endl;
						++alphaSteps;
						alpha = alpha*0.5;
					}
					
					if (alphaSteps > 10)
					{
						alphaSteps = 0;
						std::cout << "Reducing convergence conditions" << std::endl;
						c1 = c1*0.5;
						c2 = c2 * 1.1;
					}
				}
				if (alpha == alphaInit)
				{
					alphaInit *= 1.3;
				}
				else
				{
					alphaInit = alpha;
				}
				++Status.CurrentSteps;				
				double df = Functor.Value - prevF;
				minimiseContinues = CheckContinues(Functor.Gradient,df);
				prevF = Functor.Value;
				if (Status.CurrentSteps % Condition.SaveSteps == 0)
				{
					Functor.SavePosition(false,Status.CurrentSteps);
				}
				std::cout << "\t\tStep " << Status.CurrentSteps << " Taken, at Calculation Evaluation " << Functor.LoopID << " (L,Gradnorm,df) = (" <<prevF << ", " <<  Functor.Gradient.norm() << ", " << df << ")\n\t\t\t"; printTime();
			}
			
			Functor.SavePosition(true);
		}
		
		void ADAM(VectorXd &x)
		{
			
			//initialise ADAM vectors
			VectorXd m = VectorXd::Zero(Dimensions);
			VectorXd v = VectorXd::Zero(Dimensions);
			
			double beta1 = 0.9;
			double beta2 = 0.999;
			double eps = 1e-10;
			bool minimiseContinues = true;
			double prevF = 0;
			VectorXd ones = VectorXd::Constant(Dimensions,1.0);
			while (minimiseContinues)
			{
				int t = Status.CurrentSteps + 1;
				Functor.Calculate(x);
				double b1Mod = 1.0/(1.0 - pow(beta1,t));
				double b2Mod = 1.0/(1.0 - pow(beta2,t));
				m = (  beta1 *m + (1.0-beta1)*Functor.Gradient);
				
				
				VectorXd gSq = Functor.Gradient.array() * Functor.Gradient.array(); 
				v = ( beta2 * v + (1.0-beta2)* ( gSq) );
			
				VectorXd dx = b1Mod * m * Condition.StepSize;

				for (int i = 0; i < Dimensions; ++i)
				{
					dx[i] /= (sqrt(v[i]*b2Mod) + eps);
				}

				x -= dx;
				
				double df = Functor.Value - prevF;
				minimiseContinues = CheckContinues(Functor.Gradient,df);
				prevF = Functor.Value;
				++Status.CurrentSteps;
				if (Status.CurrentSteps % Condition.SaveSteps == 0)
				{
					Functor.SavePosition(false);
				}
				std::cout << "\t\tStep " << Status.CurrentSteps << " Taken, at Calculation Evaluation " << Functor.LoopID << "\n";
				std::cout << "\t\t\t(L,Gradnorm,df) = (" << std::setprecision(10) << prevF << ", " <<  std::setprecision(10) <<Functor.Gradient.norm() << ", " << std::setprecision(10) <<df << ")\n"; 
				std::cout << "\t\t\t"; printTime();
				
				if (std::isnan(Functor.Value))
				{
					exit(10);
				}
			}
		}

};



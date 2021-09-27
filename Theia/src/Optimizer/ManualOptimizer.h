#pragma once
#include <vector>
#include <math.h>
#include <ctime>
#include <iomanip>
#include <fstream>
#include <iostream>

#include "OptimizerProperties.h";

namespace ADABADAM
{
	//! A type-ambivilant implementation of the sign(x) function with sign(0) == 0 \param val a member of any type where `<` is defined \returns -1, 0 or 1 
	template <typename T> int sgn(T val) {
	    return (T(0) < val) - (val < T(0));
	}
	
	//! Calculates the standard Euclidean norm of a provided vector \param x A vector on R^n \returns ||x|| = sqrt(x.x)
	double norm(const std::vector<double> &x)
	{
		double sum = 0;
		for (int i = 0; i < x.size(); ++i)
		{
			sum += x[i]*x[i];
		}
		return sum;
	}
	
	//! Calculates the pointwise difference off a vector \param x a vector on R^n \param y a vector also on R^n \returns A vector z such that z_i = x_i - y_i
	std::vector<double> vectorDiff(const std::vector<double> &x, const std::vector<double> &y)
	{
		if (x.size() != y.size())
		{
			std::cout << "Vectors of different sizes passed to ADABADAM::vectorDiff, exit!" << std::endl;
			exit(5);
		}
		
		std::vector<double> newVec = x;
		for (int i = 0; i < x.size(); ++i)
		{
			newVec[i] -= y[i];
		}
		return newVec;
	}
	
	//! Calculates the norm of the difference between two vectors. Quicker than calculating diff() then norm(), at the cost of some duplicate code \param x a vector on R^n \param y a vector also on R^n \returns ||x -y|| = sqrt( (x-y).(x-y) )
	double vectorDiffNorm(const std::vector<double> &x, const std::vector<double> &y)
	{
		if (x.size() != y.size())
		{
			std::cout << "Vectors of different sizes passed to ADABADAM::vectorDiffNorm, exit!" << std::endl;
			exit(5);
		}
		double sum = 0;
		for (int i = 0; i < x.size(); ++i)
		{
			double d = x[i] - y[i];
			sum += d*d;
		}
		return sum;
	}
	
	//! Generates a randomly shuffled array of numbers between 0 and n-1. Suitable for then randomly looping over an object. \param n The number of random elements \returns a randomly ordered vector of unique elements between 0 and n-1
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
	
	
	
	/*!
	 * An implementation of the The ADABADAM (Adaptive Batched, Adaptive Moment) Estimator. The object is templated against the #Functor object, which contains the function that is to be minimized. 
	*/
	template<class T>
	class Optimizer
	{
		private:
			/*!
			 * A function-like class which acts as the templating property of this class. In order to work with this optimizer, the chosen functor should have the following properties:
			 *  * ``void Calculate(std::vector<double> x, int batchID, int minibatch)`` function which populates the ``Value`` and ``Gradient`` members
			 *  * ``void SavePosition(bool finalSave,int saveStep,bool uniqueSave)`` function which (optionally) saves the last used value of x to file.
			 *  * ``double Value`` member, populated by Calculate()
			 *  * ``std::vector<double>> Gradient`` member, populated by Calculate()
			 * 
			 * See LikelihoodFunctor for an example of a valid Functor object for this class.
			*/
			T & Functor;
	
			/*!Part of the  \verbatim embed:rst:inline :doc:`optimiser-data` \endverbatim. See that page for more information*/
			MemoryBuffer Buffer;
			
			//! Sets the default values for several elements of Buffer, Status, Properties and HaltConditions, so that they aren't awfully uninitialised if the user forgets to set them. Called during the constructor function.
			void SetDefaults()
			{
				HaltConditions.MaxSteps = 10000;
				HaltConditions.PositionChangeThreshold = 1e-5;
				HaltConditions.GradientThreshold = 1e-6;
				HaltConditions.FunctionChangeThreshold = 1e-7;
				HaltConditions.SingleBatchStepThreshold = 40;
				

				Properties.MiniBatches = 1;
				Properties.StepSize = 0.05;			
				
				Properties.HarnessReleaseSteps = 5;
				Properties.MaxHarnessFactor = 100;
				
				Properties.MinibatchDownStep = 2;
				
				Buffer.Size = 30;
				Buffer.AnalysisSize = 15;
				Buffer.OverrideTime = 300;
				Progress.SaveLocation = "";
				Progress.MaxHashes = 32;
				Properties.StepsPerPositionSave = 1;
				
				Properties.adamBeta1 = 0.7;
				Properties.adamBeta2 = 0.99;
			}
			
			//! Sets the relevant values of Status, Progress and Buffer to their appropriate values for the beginning of an optimisation loop. Called during #Minimize(), so any changes made after construction are erased. 
			void Initialise()
			{
				Status.TooManySteps = false;
				Status.CarryingOnRegardless = false;
				Status.Converged = false;
				Status.ReachedFunctionConvergence = false;
				Status.ReachedGradConvergence = false;
				Status.ReachedStepConvergence = false;
				Status.ExternalDownStep = false;
				Status.ExternalTermination = false;
				Progress.Harness = 1.0/Properties.MaxHarnessFactor;
				Progress.BufferFileOpened = false;
				Progress.InitialSaveComplete = false;
				
				Progress.PreviousEpoch = 9999999;
				Progress.PreviousMinibatch = 999999;
				Buffer.Position = 0;
				Buffer.StartTime = std::chrono::system_clock::now();
				Buffer.LastSaveTime = std::chrono::system_clock::now();
		
				int n = Buffer.Size;
				
				Buffer.Batches = std::vector<int>(n,0);
				Buffer.MiniBatches = std::vector<int>(n,0);
				Buffer.Fs = std::vector<double>(n,0);
				Buffer.DXs = std::vector<double>(n,0);
				Buffer.Gradnorms = std::vector<double>(n,0);
				Buffer.Times = std::vector<double>(n,0);
				Buffer.Epochs = std::vector<int>(n,0);
				
				Buffer.AnalysisSteps = 0;
				Buffer.Analysis = std::vector<double>(Buffer.AnalysisSize,0);
				
				Progress.LearningRate = Properties.StepSize;
				Progress.Hashes = 0;
				Progress.CurrentSteps = 0;
				Progress.SlowdownTriggers = 0;
				Status.Continues = true;
				
				if (HaltConditions.UseExternalInstructions)
				{
					InitialiseExternalFiles();
				}
			}
				
			//! Executes the per-minibatch portion, asking #Functor to calculate a gradient, and then using that to generate a step, dx, thereby performing a single optimization step.
			void ADABADAM_Batch(int t, std::vector<double> &m, std::vector<double> &v, std::vector<double> &x, std::vector<double> & epochGradient, double & epochL, const int batches, const std::vector<int> & batchOrder)
			{
				const double eps = 1e-10;
				int currentBatch = batchOrder[batches];

				Functor.Calculate(x,currentBatch,Progress.CurrentMinibatches);
				//save initial position
				
				
				double b1Mod = 1.0/(1.0 - pow(Properties.adamBeta1,t));
				double b2Mod = 1.0/(1.0 - pow(Properties.adamBeta2,t));			
				
				double gNorm = 0;
				double dxNorm = 0;
				
				for (int i = 0; i < Properties.Dimensions; ++i)
				{
					double g = Functor.Gradient[i];
				
					gNorm += g*g;
					m[i] = Properties.adamBeta1 * m[i] + (1.0 - Properties.adamBeta1)*g;
					v[i] = Properties.adamBeta2* v[i] + (1.0 - Properties.adamBeta2) * (g*g);
					
					double effectiveRate = Progress.LearningRate * Progress.Harness;
					
					double dx_i = - b1Mod * m[i] /  ((sqrt(v[i]*b2Mod) + eps) ) * effectiveRate;
					
					dxNorm += dx_i * dx_i;
					x[i] += dx_i;
					epochGradient[i] += g;
				}	
				
				++t;
				
				epochL += Functor.Value;

				if (Progress.CurrentMinibatches > 1)
				{
					double df_mini = Functor.Value - Progress.PreviousMinibatch;
					Progress.PreviousMinibatch = Functor.Value;
					double sqrtgNorm = sqrt(gNorm);
					UpdateBuffer(batches,Progress.CurrentMinibatches,Functor.Value,sqrtgNorm,df_mini,sqrt(dxNorm),x);
				}
				if (!Progress.InitialSaveComplete)
				{
					Progress.PreviousEpoch = Functor.Value;
					Functor.SavePosition(false,0,Properties.UniquePositionSaves);
					Progress.InitialSaveComplete = true;
				}
				double harnessFactor = pow(Properties.MaxHarnessFactor, 1.0/(Properties.HarnessReleaseSteps * Progress.CurrentMinibatches));

				Progress.Harness = std::min(1.0,Progress.Harness * harnessFactor);
			};
			
			
			//! Executes the per-Epoch loop: generates a random order for the minibatches and then loops over them, calling ADABADAM_Batch(), before performing convergence and downstep checks.
			void ADABADAM_Epoch(int & t, std::vector<double> &m, std::vector<double> &v, std::vector<double> & epochGradient, std::vector<double> &x, std::vector<double> & oldX)
			{
				
				int epochs = Progress.CurrentSteps + 1;
				
				double epochL = 0;
				std::fill(epochGradient.begin(), epochGradient.end(),0.0);
								
				std::vector<int> batchOrder = randomShuffle(Progress.CurrentMinibatches);

		
				for (int batches = 0; batches < Progress.CurrentMinibatches; ++batches)
				{
					ADABADAM_Batch(t,m,v,x,epochGradient, epochL, batches,batchOrder);		
				}
				
				epochL/=Progress.CurrentMinibatches;
				double df = epochL - Progress.PreviousEpoch;
				Progress.PreviousEpoch = epochL;
				double epochGradNorm = norm(epochGradient)/sqrt(Progress.CurrentMinibatches);
				double epochDx = vectorDiffNorm(x,oldX);
				oldX = x;
				
				++Progress.CurrentSteps;
				CheckExternalFiles();
				UpdateBuffer(-1,Progress.CurrentMinibatches,epochL,epochGradNorm,df,epochDx,x);
				
				double newBatches = CheckMinibatches(df,Progress.CurrentMinibatches);
				
				bool justSlowed = false;
				if (newBatches < Progress.CurrentMinibatches)
				{
					Progress.CurrentMinibatches = newBatches;
					Progress.Harness = 1.0/Properties.MaxHarnessFactor;
					justSlowed = true;
				}
				
				
				if (df > 0 && justSlowed == false)
				{
					Progress.LearningRate *= 0.5;
					++Progress.SlowdownTriggers;
				}
				if (df < 0)
				{
					Progress.LearningRate *= 1.02;
					if (Progress.LearningRate > Properties.StepSize)
					{
						Progress.LearningRate = Properties.StepSize;
					}
				}
				
				
				CheckConvergence(norm(epochGradient),df,epochDx);
				if (Status.Continues == false && (Progress.CurrentMinibatches > 1 || Progress.Harness < 1) && !Status.ExternalTermination)
				{
					Status.Continues = true;
					Progress.CurrentMinibatches = std::max(1,Progress.CurrentMinibatches/2);
					if (Progress.Harness ==1)
					{
						Progress.Harness = 1.0/Properties.MaxHarnessFactor;
					}
				}
					
				
			}
				
			
			/*! The bulk of the ADABADAM machinery, excuting the main workflow. For a discussion about what this function does, see  \verbatim embed:rst:inline :ref:`adabadam-theory` \endverbatim. \param x The initial starting position to optimize towards. */
			void ADABADAM_Body(std::vector<double> &x)
			{
				Progress.CurrentMinibatches = Properties.MiniBatches;
	
				int n = x.size();
				Properties.Dimensions = n;
				
				//initialise ADAM vectors
				std::vector<double> m(n,0.0);
				std::vector<double> v(n,0.0);
				std::vector<double> epochGradient(n,0.0);
				std::vector<double> oldX = x;
				int t = 1;
	
		
				bool initSaved = false;
				while (Status.Continues && ~Status.CarryingOnRegardless)
				{
					ADABADAM_Epoch(t, m, v, epochGradient, x, oldX);
				}
				SaveBuffer(Buffer.Position);
				Functor.SavePosition(true,0,Properties.UniquePositionSaves);
				CleanExternalFiles();
			}
				
			//! Check the current progress of the optimizer against the StopConditions to see if convergence has been reached. Updates the values of Status if the conditions have been met. \param dg The norm of the Gradient \param df The change in the functional value across the epoch \param dx The norm of the total dx change across the epoch
			void CheckConvergence( double dg, double df, double dx)
			{
				if (Progress.CurrentSteps > HaltConditions.MaxSteps)
				{
					Status.TooManySteps = true;
					Status.Continues = false;
				}
				if (HaltConditions.PositionChangeThreshold > 0 && dx < HaltConditions.PositionChangeThreshold && Progress.Harness > 0.95)
				{
					Status.ReachedStepConvergence = true;
					Status.Converged = true;
					Status.Continues = false;
				}
				if (HaltConditions.GradientThreshold > 0 && dg < HaltConditions.GradientThreshold)
				{
					Status.ReachedGradConvergence = true;
					Status.Converged = true;
					Status.Continues = false;
				}
				if (HaltConditions.FunctionChangeThreshold > 0  && abs(df) < HaltConditions.FunctionChangeThreshold)
				{
					Status.ReachedFunctionConvergence = true;
					Status.Converged = true;
					Status.Continues = false;
				}
				if (Status.ExternalTermination == true)
				{
					Status.Continues = false;
					Status.Converged = false;
				}
				
				Status.CarryingOnRegardless = false;
				if (Status.Continues == false && Progress.EpochsSinceSingleBatch < HaltConditions.SingleBatchStepThreshold)
				{
					Status.CarryingOnRegardless = true || !Status.ExternalTermination;
				}
			} 
	
			//! Looks for a reason to reduce the number of minibatches, either due to NeedsBatchReduction() returning true, or because OptimizerStatus::ExternalDownStep is true. Divides the current minibatches by OptimizerProperties::MinibatchDownStep \param df The change in the functional value across the epoch \param currentSize the current number of minibatches used \returns the number of minibatches which will be used going forwards
			int CheckMinibatches(double df,int currentSize)
			{
			
				int analysisPos = Buffer.AnalysisSteps % Buffer.AnalysisSize; 
				
				Buffer.Analysis[analysisPos] = df;
				++Buffer.AnalysisSteps;	
	
				double newSize = currentSize;
				
				if (Status.ExternalDownStep || Buffer.AnalysisSteps >= Buffer.AnalysisSize)
				{
					if (NeedsBatchReduction() || Status.ExternalDownStep)
					{	
						Buffer.AnalysisSteps = 0;
						Progress.SlowdownTriggers = 0;
						newSize = currentSize / Properties.MinibatchDownStep;
						if (newSize < 1)
						{
							newSize = 1;
						}
					}
				
				}		
				if (currentSize == 1)
				{
					++Progress.EpochsSinceSingleBatch;
				}
				return newSize;
			}
	
			/*!
			 * Looks through the MemoryBuffer::Analysis vector to look for problematic trends which would indicate that a minibatch reduction is needed. These triggers are: 
			 * 1. More than 3 sign changes in df across the buffer's duration (i.e. oscillatory behaviour)
			 * 2. Mean df across buffer duration is positive (i.e. not minimizing)
			 * 3. More than 2 ProgressTracker::Harness triggers since the last minibatch downchange
			*/
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
				if ((mean > 0) || signChanges >= problematicSignChanges || Progress.SlowdownTriggers > 2)
				{
					batchesAreAProblem = true;
				}
	
	
				return batchesAreAProblem;
			}
	
			//! Inserts relevant values into the Buffer. If the buffer is full, or it was a long time (c.f. MemoryBuffer::OverrideTime) since the last save, it then calls SaveBuffer(). Also included in this are the progress bars and output to the terminal. 
			void UpdateBuffer(int batch, int nBatches,double F, double G, double dF,double dxNorm,const std::vector<double> & x)
			{
				int i = Buffer.Position;
				Buffer.MiniBatches[i] = batch;
				Buffer.Batches[i] = nBatches;
				Buffer.Fs[i] = F;
				Buffer.Gradnorms[i] = G;
				Buffer.DXs[i] = dxNorm;
				Buffer.Epochs[i] = Progress.CurrentSteps + 1;
				auto time = std::chrono::system_clock::now();
				std::chrono::duration<double> diff = time - Buffer.StartTime;
				
				Buffer.Times[i] = diff.count();
				
				++Buffer.Position;
				
				std::chrono::duration<double> savediff = time - Buffer.LastSaveTime;
				if (Buffer.Position >= Buffer.Size || savediff.count() >= Buffer.OverrideTime)
				{
					SaveBuffer(std::min(Buffer.Size,Buffer.Position));
					Buffer.Position = 0;
					Buffer.LastSaveTime = time;
				}
				
				if (batch == -1)
				{

					
					if (nBatches == 1)
					{
						std::cout << "\t\tEpoch " << Progress.CurrentSteps;
					}
					else
					{
						std::cout << "]";
					}
					std::cout << " complete\n";
					std::cout << "\t\t\t(L,Gradnorm,dL,|dx|,nBatch) = (" << std::setprecision(10) << F << ", " <<  std::setprecision(10) << G << ", " << std::setprecision(10) << dF << ", " << std::setprecision(10) << dxNorm << ", " << nBatches <<")\n"; 
					
					
					
					if (Progress.CurrentSteps % Properties.StepsPerPositionSave == 0)
					{
						Functor.SavePosition(false,Progress.CurrentSteps,Properties.UniquePositionSaves);
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
	
			//! Called by UpdateBuffer(), this saves the contents of the Buffer into the file given by the ProgressTracker::SaveLocation member of Progress.
			void SaveBuffer(int n)
			{
				std::string saveFile = Progress.SaveLocation + "OptimizerProgress.txt";
				std::fstream file;
				int width = 20;
				if (Progress.BufferFileOpened == false)
				{
					file.open(saveFile,std::ios::out);
					std::vector<std::string> headers = {"Elapsed","Epoch","Batch", "nBatches","F","dX","GradNorm"};
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
					file << std::setw(width) << std::setprecision(prec) << Buffer.DXs[i] << ",";
					file << std::setw(width) << std::setprecision(prec) << Buffer.Gradnorms[i] << ",";
					
					file << "\n";
				}
				
				file.close();
			}
	
			//! Creates and writes the default values (and explanatory text) to the StopConditions::DownStepFile and StopConditions::TerminationFile of HaltConditions.
			void InitialiseExternalFiles()
			{
				
				std::vector<std::string> names = {HaltConditions.TerminationFile, HaltConditions.DownStepFile};
				std::vector<std::string> purposes = {"terminate the optimisation process", "trigger a minibatch downstep"};
				std::fstream file;
				
				for (int i = 0; i < names.size(); ++i)
				{
					file.open(names[i],std::ios::out);
	
					file << "#### This file can be used to " << purposes[i] << ": set the value below to 1 to trigger the action at the end of the current epoch ####\n0";
					file.close();
				}
				
			}
			
			//! Reads the default values (and explanatory text) to the StopConditions::DownStepFile and StopConditions::TerminationFile of HaltConditions. If they are 1, sets the relevant flags OptimizerStatus::ExternalDownstep and OptimizerStatus::ExternalTermination respectively.
			void CheckExternalFiles()
			{
				if (HaltConditions.UseExternalInstructions)
				{
					std::vector<std::string> names = {HaltConditions.TerminationFile, HaltConditions.DownStepFile};
		
					//~ Status.ExternalTermination = false;
					//~ Status.ExternalDownStep = false;
					for (int i = 0; i < names.size(); ++i)
					{
						int j = 0;
						forLineIn(names[i],
							if (j != 0)
							{
								bool fileActivated = bool(std::stod(FILE_LINE));
								if (i == 0)
								{
									Status.ExternalTermination = fileActivated;
								}
								else
								{
									Status.ExternalDownStep = fileActivated;
								}
							}
							++j;
						);
					}
					InitialiseExternalFiles();
				}
			}
			
			//! After the optimizer has completed, delete the files created by InitialiseExternalFiles(). It's always good to clean up your mess.
			void CleanExternalFiles()
			{
				if (HaltConditions.UseExternalInstructions)
				{
					std::vector<std::string> names = {HaltConditions.TerminationFile, HaltConditions.DownStepFile};
					for (int i = 0; i < names.size(); ++i)
					{
						std::string command = "rm " + names[i];
						system(command.c_str());
					}
				}
			}
		public:
			
			/*!Part of the  \verbatim embed:rst:inline :doc:`optimiser-data` \endverbatim. See that page for more information*/
			OptimizerStatus Status;
			
			/*!Part of the  \verbatim embed:rst:inline :doc:`optimiser-data` \endverbatim. See that page for more information*/
			OptimizerProperties Properties;
			
			/*!Part of the  \verbatim embed:rst:inline :doc:`optimiser-data` \endverbatim. See that page for more information*/
			StopConditions HaltConditions;
			
			/*!Part of the  \verbatim embed:rst:inline :doc:`optimiser-data` \endverbatim. See that page for more information*/
			ProgressTracker Progress;
			
			
			//! Constructor class. Doesn't do much except initialise #Functor and call SetDefaults()
			Optimizer<T>(T& functor) : Functor(functor)
			{
				SetDefaults();	
			}
		
			//! Calls Initialise() to properly clean everything before launching into the meat of matters: actually performing the ADABADAM_Body() optimization.
			void Minimize(std::vector<double> & x)
			{
				Initialise();
				ADABADAM_Body(x);
			}
	
			//! Returns a nicely formatted string detailing the contents of #Status. Useful for outputting a debrief after optimization has stopped.
			std::string GetStatus()
			{
				std::string s = "";
				s += "Steps Taken: " + std::to_string(Progress.CurrentSteps) + " / " + std::to_string(HaltConditions.MaxSteps);
				s += "\nHalt conditions: ";
				std::vector<std::string> titles = {"Too many steps", "Reached Gradient Convergence", "Reached Step Convergence", "Reached Functional Convergence","External Quit"};
				std::vector<bool> values = {Status.TooManySteps, Status.ReachedGradConvergence, Status.ReachedStepConvergence, Status.ReachedFunctionConvergence,Status.ExternalTermination};
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
}


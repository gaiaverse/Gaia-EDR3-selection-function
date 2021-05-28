#pragma once
#include <vector>
#include <mpi.h>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"

//~ #include "Star.h"
#include "../Likelihood/LogLikelihoodPrior.h"

#include "../GenericFunctions/timeCodes.h"

using Eigen::VectorXd;

//DescentFunctor is a function-like class which acts as a wrapper for the gradient descent algorithm. 
//The overloaded operator () allows the class to be called as a function by LBFGs, but the classlike nature allows the function to 
//access data without needing to continually reload it into the function


class DescentFunctor
{
	private:
		int RunningID;

		//since the descent functor runs only once (i.e. on Root), we would still like root to use its CPU cycles to do some calculating, so we have a copy of 
		//the structures needed to do liklihood analysis stored within
		const std::vector<std::vector<Star>> &Data; 
		LogLikelihoodPrior L;

		
		std::chrono::time_point<std::chrono::system_clock> Start;

		//running values for the loglikelihood and gradient 

		
		std::string OutputDir;
		//to prevent double evaluations at the same point, prevlock saves current position against a threshold
		VectorXd PrevLock;
		const double lockLim = 1e-15;
		
		
		//Needlet stuff - has to be public
		int needletN;
		std::vector<int> needlet_u;
    	std::vector<int> needlet_v;
    	std::vector<double> needlet_w;
  
		int NStars;
		int StarsInLastBatch;
		//holder for transformed values
		std::vector<double> TransformedPosition;
		std::vector<double> TransformedGradient;
		std::vector<double> bVector;
		std::vector<double> mut_gaps;
		
		std::vector<bool> freezeOuts;
		std::vector<bool> freezeOuts_mag; 

		void ForwardTransform(const VectorXd &z);
		void BackwardTransform();		
		void ResetPosition();
		

	public:
		int LoopID;

		double Value;
		std::vector<double> Gradient;
	
	    DescentFunctor(int n,const std::vector<std::vector<Star>> & data, int nParams,std::string outdir, int nStars): Data(data), L(data, n)
	    {
				NStars = nStars;
				RunningID = n;
				LoopID = 0;
				Start = std::chrono::system_clock::now();
								
				PrevLock = VectorXd::Random(totalRawParams);
				
				TransformedPosition = std::vector<double>(totalTransformedParams,0);
				TransformedGradient = std::vector<double>(totalTransformedParams,0);
				OutputDir = outdir;
				
				std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
				int i = 0;
			    forLineVectorInFile(needlet_file,',',
			 
					if (i > 0)
					{
				        needlet_u.push_back(std::stoi(FILE_LINE_VECTOR[0]));
				        needlet_v.push_back(std::stoi(FILE_LINE_VECTOR[1]));
				        needlet_w.push_back(std::stod(FILE_LINE_VECTOR[2]));
					}
			        ++i;
			    );    
			    needletN = needlet_u.size();
			    bVector = std::vector<double>(Nm*Ns,0);
			    
			    
				Value = 0;
				Gradient = std::vector<double>(totalRawParams,0);
				
				
				mut_gaps = std::vector<double>(Nt,0);
				std::string gapFile = "../../ModelInputs/gaps_prior.dat";
				double timeFactor = (double)TotalScanningTime / Nt;
				int it = 0;
				bool inGap = false;
				int borderWidth = 0;
				int modifiedBorderWidth = borderWidth * timeFactor;
				bool inBorder= false;
				int trueTime = 0;
				int lastEnd = -9999;
				freezeOuts = std::vector<bool>(Nt,true);
				freezeOuts_mag = std::vector<bool>(Nt_m,true);
				forLineVectorInFile(gapFile,' ',
					
					int gapStart = std::stoi(FILE_LINE_VECTOR[0]);
					int gapEnd = std::stoi(FILE_LINE_VECTOR[1]);
					
					trueTime = floor(it * timeFactor);
					while (trueTime < gapEnd)
					{
						int leftDistance = std::min(abs(trueTime - gapStart),abs(trueTime - lastEnd));
						int rightDistance = abs(trueTime - gapEnd);
						
						bool inGap = (trueTime >= gapStart) && (trueTime <= gapEnd);
						
						bool nearGapEdge = (leftDistance < modifiedBorderWidth) || (rightDistance < modifiedBorderWidth);
						double insertValue = mut_normal;
						if (inGap)
						{
							insertValue = mut_gap;
							//~ freezeOuts[it] = true;

						}
						if (nearGapEdge)
						{
							insertValue = mut_border;
						}
					
						mut_gaps[it] = insertValue;
						//~ std::cout << "\t " <<it << "  " << trueTime << "   " << insertValue << "   " << leftDistance << "   " << rightDistance << std::endl;
						
						++it;
						trueTime = floor((double)it * timeFactor);
						
					}
					//~ std::cout << "Gap finished at it = " << it << " t = " << trueTime << std::endl;
					lastEnd = gapEnd;
				);
				
				while (it<Nt)
				{
					mut_gaps[it] = mut_normal;
					++it;
				}
			
		}
	    void DistributeCalculations(const VectorXd &y, int batchID, int effectiveBatches);
 
		
		void Calculate(const VectorXd &x, int batchID, int effectiveBatches);
		void Calculate(const VectorXd &x);
		void SavePosition(bool finalSave, int saveStep);
		void Unfreeze();
};



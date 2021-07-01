#pragma once
#include <vector>
#include <mpi.h>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#define EIGEN_MPL2_ONLY
#include "../libs/Eigen/Core"

//~ #include "Star.h"
#include "../Likelihood/LogLikelihoodPrior.h"
#include "../Main/GlobalConstants.h"
#include "../libs/JSL/JSL.h"
#include "../Likelihood/ProbabilityFunctions.h"

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

		
		//Needlet stuff - has to be public
		int needletN;
		std::vector<int> needlet_u;
    	std::vector<int> needlet_v;
    	std::vector<double> needlet_w;
  
		int NStars;
		int StarsInLastBatch;
		
		
		bool SpaceActive;
		bool TimeActive;
		bool HyperActive;
		int ActiveParams;
		//holder for transformed values
		std::vector<double> TransformedPosition;
		std::vector<double> TransformedGradient;
		std::vector<double> bVector;
		
		
		std::vector<bool> freezeOuts;
		//~ std::vector<bool> freezeOuts_mag; 

		void ForwardTransform(const VectorXd &z);
		void ForwardTransform_Spatial(const VectorXd &z);
		void ForwardTransform_Temporal(const VectorXd &z);
		void ForwardTransform_Hyper(const VectorXd &z);

		void BackwardTransform();		
		void BackwardTransform_Spatial();
		void BackwardTransform_Temporal();
		void BackwardTransform_Hyper();
		
		void ResetPosition();
		
		int MaxBatches;
	public:
		int LoopID;

		double Value;
		std::vector<double> Gradient;
		
		std::vector<double> mut_gaps;
		std::vector<double> FrozenSpace;
		std::vector<double> FrozenTime;
		std::vector<double> FrozenHypers;
		
	    DescentFunctor(int n,const std::vector<std::vector<Star>> & data, int nParams,std::string outdir, int nStars, int maxBatches, bool timeActive, bool spaceActive, bool hyperActive): Data(data), L(data, n)
	    {
			TimeActive = timeActive;
			SpaceActive = spaceActive;
			HyperActive = hyperActive;
			
			NStars = nStars;
			RunningID = n;
			LoopID = 0;
			Start = std::chrono::system_clock::now();
			MaxBatches = maxBatches;		
			
			TransformedPosition = std::vector<double>(totalTransformedParams,0);
			TransformedGradient = std::vector<double>(totalTransformedParams,0);
			OutputDir = outdir;
			
			std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
			int i = 0;
		    forLineVectorIn(needlet_file,',',
		 
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
			ActiveParams = Nt*TimeActive + Ns*Nm*SpaceActive + NHyper*HyperActive;
			FrozenTime = std::vector<double>(Nt,0.0);
			FrozenSpace = std::vector<double>(Nm*Nl,0.0);
			FrozenHypers = std::vector<double>(NHyper,0.0);
			Gradient = std::vector<double>(ActiveParams,0);
			
			
			mut_gaps = std::vector<double>(Nt,xtPriorNonGap);
			
			for (int i = 0; i < Nt; ++i)
			{
				if (GapList[i] == true)
				{
					mut_gaps[i] = xtPriorInsideGap;
				}
			}
		}
	    void DistributeCalculations(const VectorXd &y, int batchID, int effectiveBatches);
 
		
		void Calculate(const VectorXd &x, int batchID, int effectiveBatches);
		void Calculate(const VectorXd &x);
		void SavePosition(bool finalSave, int saveStep, bool uniqueSave, const VectorXd & x);
		void Unfreeze();
};



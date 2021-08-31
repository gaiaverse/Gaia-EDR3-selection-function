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

#include "EfficiencyVector.h"


/*!
 * A function-like class which acts as a wrapper for the function calls used by the gradient descent algorithm and the necessary message passing for communication with the workers. The overloaded operator () allows the class to be called as a function, but it can store data preventing the need for reloading the data or passing huge numbers of arguments.
*/
class LikelihoodFunctor
{
	private:

		//! Allows the calling process to calculate its own portion of the Likelihood, and as a guaranteed priveleged worker, also can call the associated priors.
		LogLikelihoodPrior L;
		
		//! The maximum number of minibatches which are used (and which therefore define the structure of #LogLikelihoodPrior::Data)
		int MaxBatches;
		
	public:
		//! The total number of calls that have been passed to DistributeCalculations()
		unsigned int LoopID;
		
		//! The current proposed efficiency vector
		EfficiencyVector Efficiency;
		
		//! A required member for templating values of the ADABADAM::Optimizer class. After Calculate() is called, holds  current value of the Functor. *Note that we perform an explicit renormalisation of this value, dividing through by the number of stars used to calculate it (hence it is the approximate per-star Value). In addition, as ADABADAM attempts to minimise functions, whilst the Likelihood returns a value we wish to* **maximise**, *so this is the negation of that value*.
		double Value;
		
		//! A required member for templating values of the ADABADAM::Optimizer class. After Calculate() is called, holds  current value of the Gradient. *Note that we perform an explicit renormalisation of this value, dividing through by the number of stars used to calculate it (hence it is the approximate per-star Value). In addition, as ADABADAM attempts to minimise functions, whilst the Likelihood returns a value we wish to* **maximise**, *so this value is equal to the product of -1 and* EfficiencyVector::RawGradient.
		std::vector<double> Gradient;
		
	
		//! Constructor class -- initialises #L and #MaxBatches. \param data The data, arranged according to the minibatching schedule. \param maxBatches the original number of minibatches which are used
	    LikelihoodFunctor(const std::vector<std::vector<Star>> & data, int maxBatches): L(data), MaxBatches(maxBatches)
		{}
		
		/*! A separate chunk of almost-constructor, initialising the remaning components + getting the starting position. An argument could be made that this should be part of the constructor, but it is kept separate so that you can re-run the optimisation with a different starting position without destroying the object, you can just reinitialise it. 
		 *\param loadPosition The place for #Efficiency to check for valid starting positions
		 * \param outdir The output directory where savefiles are located
		 * \returns The EfficiencyVector::RawPosition component of #Efficiency,  the default starting point for the optimizer
		*/
		std::vector<double> Initialise(std::string loadPosition, std::string outdir)
		{
			LoopID = 0;
			Efficiency = EfficiencyVector(loadPosition,outdir);
			Gradient = std::vector<double>(totalRawParams,0.0);
			Value = 0;
			return Efficiency.RawPosition;
		}
		

		/*!
		 * A required function call for templating values of the ADABADAM::Optimizer class, and the main work loop of the class. Sends MPI messages to the workers and instructs them to calculate the next value of their LogLikelihood object. Then uses MPI to collect these objects, call the associated Priors and then populate the #Value and #Gradient objects for the optimizer to use.  
		 * \param x The new value of the EfficiencyVector::RawPosition, the current proposed operating efficiency
		 * \param batchID The current (randomised) batchID from which to draw the population of stars used to calculate the likelihood. Can be 0 <= batchID < effectiveBatches
		 * \param effectiveBatches The current number of minibatches used by the optimizer
		 * \returns No explicit return, but populates the #Value and #Gradient objects
		*/
		void Calculate(const std::vector<double> &x, int batchID, int effectiveBatches);
		
		//! An overloading alias in case minibatching is ever disabled: calls Calculate(x,0,1)
		void Calculate(const std::vector<double> &x);
		
		/*!
		 *  A required function call for templating values of the ADABADAM::Optimizer class. Saves the current value of Efficiency to file via the EfficiencyVector::Save() call.
		*/
		void SavePosition(bool finalSave, int saveStep, bool uniqueSave);
		
};



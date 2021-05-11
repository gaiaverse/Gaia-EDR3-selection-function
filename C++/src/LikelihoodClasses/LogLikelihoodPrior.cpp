#include "LogLikelihoodPrior.h"


void LogLikelihoodPrior::PriorCalculate(const std::vector<double> & x, int effectiveBatchID, int effectiveBatches)
{

	Calculate(x,effectiveBatchID,effectiveBatches);
	
	
	
	//~ int n = totalRawParams;
	//~ for (int i = 0; i < n ; ++i)
	//~ {
		//~ double d = x[i];
		//~ Value -= 0.5 * d * d / effectiveBatches;
		
		//~ Gradient[i] -= d / effectiveBatches;
	//~ }
	
}



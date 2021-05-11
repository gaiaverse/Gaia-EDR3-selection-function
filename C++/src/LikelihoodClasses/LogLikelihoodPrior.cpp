#include "LogLikelihoodPrior.h"


void LogLikelihoodPrior::PriorCalculate(const std::vector<double> & x, int effectiveBatchID, int effectiveBatches)
{
	Calculate(x,effectiveBatchID,effectiveBatches);
	
	
	int n = Nt + Nm*Ns;
	for (int i = 0; i < n ; ++i)
	{
		double d = x[i];
		Value -= 0.5 * d * d / effectiveBatches;
		
		Gradient[i] -= d / effectiveBatches;
	}
}



#include "LogLikelihoodPrior.h"


void LogLikelihoodPrior::Prior(const Eigen::VectorXd& RawParams, double * currentValue, Eigen::VectorXd * currentGradient, int effectiveBatches)
{
	
	int n = Nt + Nm*Ns;
	for (int i = 0; i < n ; ++i)
	{
		double d = RawParams[i];
		currentValue[0] -= 0.5 * d * d / effectiveBatches;
		
		currentGradient[0][i] -= d / effectiveBatches;
	}
}


void LogLikelihoodPrior::MakeCovarianceMatrix()
{
	Eigen::Matrix<double, Nm, Nm> Kg;
	for (int i = 0; i < Nm; i++) 
	{
		for (int j = 0; j < i; j++) 
		{
			Kg(i,j) = Kg(j,i) = exp(-pow(i - j,2)/(2.0*lm*lm));
		}
		Kg(i,i) = 1.0 + SingularityPreventer;
	}

	//decompose to make CholeskyKg
	CholeskyKg = Kg.llt().matrixL();
	Kg_decomposed = true;
}

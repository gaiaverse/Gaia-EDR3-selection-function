#include "LogLikelihoodPrior.h"


void LogLikelihoodPrior::Prior(const Eigen::VectorXd& RawParams, double * currentValue, std::vector<double> * currentGradient, int effectiveBatches)
{
	int n = totalRawParams;
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

	std::vector<double> max_in_row = std::vector<double>(Nm,0);
	for (int i = 0; i < Nm; i++) 
	{
		for (int j = 0; j <= i; j++) 
		{
			max_in_row[i] += std::max(0.0,abs(CholeskyKg(i,j))-max_in_row[i]);
		}
	}

	choleskyN = 0;
	for (int i = 0; i < Nm; i++) 
	{
		for (int j = 0; j <= i; j++) 
		{
			if (abs(CholeskyKg(i,j)) > cholesky_tol * max_in_row[i])
			{
				choleskyN += 1;
				cholesky_u.push_back(i);
				cholesky_v.push_back(j);
				cholesky_w.push_back(CholeskyKg(i,j));
			}
		}
	}

	Kg_decomposed = true;
}

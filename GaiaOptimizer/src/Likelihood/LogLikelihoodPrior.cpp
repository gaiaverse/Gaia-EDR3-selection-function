#include "LogLikelihoodPrior.h"


void LogLikelihoodPrior::Prior(const Eigen::VectorXd& RawParams, double * currentValue, std::vector<double> * currentGradient, int effectiveBatches, bool space, bool time, bool hyper)
{
	int n = rawNonHyperParams;
	int N = n;
	if (useHyperPrior == true)
	{
		N = totalRawParams;
	}
	
	if (time)
	{
		for (int i = 0; i < Nt; ++i)
		{

			double d = RawParams[i];
			currentValue[0] -= 0.5 * d * d / effectiveBatches;
			
			currentGradient[0][i] -= d / effectiveBatches;
		}
	}
	if (space)
	{
		int spaceOffset = 0;
		if (time)
		{
			spaceOffset += Nt;
		}
		for (int i = 0; i < Nm*Ns; ++i)
		{
			double d = RawParams[spaceOffset + i];
			currentValue[0] -= 0.5 * d * d / effectiveBatches;
			
			currentGradient[0][spaceOffset + i] -= d / effectiveBatches;
		}
	}
	
	
	if (hyper && useHyperPrior)
	{
		int hyperOffset = 0;
		if (time)
		{
			hyperOffset += Nt;
		}
		if (space)
		{
			hyperOffset += Ns*Nm;
		}
		
		for (int i = 0; i < hyperOrder + 1;++i)
		{
			for (int j = 0;j < NVariancePops; ++j)
			{
				double mean = 0;
				if (i > 0)
				{
					mean = log(1e-3);
				}
				double d = (RawParams[hyperOffset + i*NVariancePops + j] - mean);
				currentValue[0] -= 0.5 * d * d / effectiveBatches;
				
				currentGradient[0][hyperOffset + i*NVariancePops + j] -= d / effectiveBatches;
			}
		}
		
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

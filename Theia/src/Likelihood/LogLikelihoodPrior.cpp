#include "LogLikelihoodPrior.h"

struct F_dF
{
	double F;
	double dF;
	F_dF()
	{
		F = 0;
		dF = 0;
	}
	F_dF(double f, double df) : F(f), dF(df){};
};


F_dF Normal(double x, double mu, double sigma)
{
	double d = (x - mu);
	double s2 = sigma*sigma;
	double v = -0.5 * d*d/(s2);
	double dv = - d/s2;
	
	return F_dF(v,dv);
}

F_dF StudentT(double x, double mu, double nu)
{
	double div = 2;
	double d = (x - mu)/div;
	double v = - (nu + 1)/2 * log(1 + d*d/nu);
	double dv = - (nu + 1) * d/ ( nu * (1 + d*d/nu) * div);
	
	return F_dF(v,dv); 
}

F_dF GapEnforcer(double x)
{
	double v = gapPriorAlpha * x - (gapPriorBeta + gapPriorAlpha)* log(1.0 + exp(x));
	double dv = (gapPriorAlpha - gapPriorBeta* exp(x))/(1 + exp(x));
	return F_dF(v,dv);
}

void LogLikelihoodPrior::RawPrior(const std::vector<double>& RawParams, double * currentValue, std::vector<double> * currentGradient, int effectiveBatches)
{
	
	
	for (int i = 0; i < Nt; ++i)
	{
		F_dF p = StudentT(RawParams[i],0,studentNu);
		currentValue[0] += p.F / effectiveBatches;
		currentGradient[0][i] += p.dF / effectiveBatches;
	}
	

	for (int i = 0; i < Nm*Ns; ++i)
	{
		double d = RawParams[Nt + i];
		F_dF p = Normal(RawParams[Nt + i],0,1);
		currentValue[0] += p.F / effectiveBatches;
		currentGradient[0][Nt + i] += p.dF / effectiveBatches;
	}
	
	
	if ( useHyperPrior)
	{
		for (int i = 1; i < hyperOrder+1; ++i)
		{
			for (int j = 0; j < NVariancePops; ++j)
			{
				int index = rawNonHyperParams + i * NVariancePops + j;
				double d = RawParams[index];
				F_dF p = Normal(d,0,1);
				currentValue[0] += p.F / effectiveBatches;
				double old = currentGradient[0][index];
				currentGradient[0][index] += p.dF / effectiveBatches;
				
			}
		}
	}		
}

void LogLikelihoodPrior::TransformPrior(const std::vector<double> & TransformPosition, double * currentValue, std::vector<double> & TransformGradient, int effectiveBatches)
{
	
	for (int i = 0; i < Nt; ++i)
	{
		if (BufferedGapList[i])
		{
			F_dF p = GapEnforcer(TransformPosition[i]);
			currentValue[0] += p.F / effectiveBatches;
			TransformGradient[i] += p.dF / effectiveBatches;
		}
	}

}


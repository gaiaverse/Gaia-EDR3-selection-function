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

void LogLikelihoodPrior::RawPrior(const Eigen::VectorXd& RawParams, double * currentValue, std::vector<double> * currentGradient, int effectiveBatches, bool space, bool time, bool hyper)
{
	
	int spaceOffset = 0;
	int hyperOffset = 0;
	if (time)
	{
		spaceOffset += Nt;
		hyperOffset += Nt;
	}
	if (space)
	{
		hyperOffset += Ns*Nm;
	}
	
	if (time)
	{
		for (int i = 0; i < Nt; ++i)
		{
			F_dF p = StudentT(RawParams[i],0,studentNu);
			currentValue[0] += p.F / effectiveBatches;
			currentGradient[0][i] += p.dF / effectiveBatches;
		}
	}
	if (space)
	{
		for (int i = 0; i < Nm*Ns; ++i)
		{
			double d = RawParams[spaceOffset + i];
			F_dF p = Normal(RawParams[spaceOffset + i],0,1);
			currentValue[0] += p.F / effectiveBatches;
			currentGradient[0][spaceOffset + i] += p.dF / effectiveBatches;
		}
	}
	
	if (hyper && useHyperPrior)
	{	
		for (int i = 1; i < hyperOrder+1; ++i)
		{
			for (int j = 0; j < NVariancePops; ++j)
			{
				int index = hyperOffset + i * NVariancePops + j;
				double d = RawParams[index];
				F_dF p = Normal(d,0,1);
				currentValue[0] += p.F / effectiveBatches;
				double old = currentGradient[0][index];
				currentGradient[0][index] += p.dF / effectiveBatches;
				
			}
		}
			
	}
}

void LogLikelihoodPrior::TransformPrior(const std::vector<double> & TransformPosition, double * currentValue, std::vector<double> & TransformGradient, int effectiveBatches,bool space, bool time, bool hyper)
{
	if (time)
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

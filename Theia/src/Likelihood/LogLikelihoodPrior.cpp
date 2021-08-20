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

double LogLikelihoodPrior::RawPrior(EfficiencyVector &x, int effectiveBatches)
{
	
	double L = 0;
	for (int i = 0; i < Nt; ++i)
	{
		double d = x.Access(x.Raw,x.Temporal,x.Position,i);
		F_dF p = StudentT(d,0,studentNu);
		L += p.F / effectiveBatches;
		x.Increment(x.Raw,x.Temporal,x.Gradient,i,p.dF / effectiveBatches);
	}
	

	for (int i = 0; i < Nm*Ns; ++i)
	{
		double d = x.Access(x.Raw,x.Spatial,x.Position,i);
		
		F_dF p = Normal(d,0,1);
		L += p.F / effectiveBatches;
		x.Increment(x.Raw,x.Spatial,x.Gradient,i, p.dF / effectiveBatches);
	}
	
	
	//~ if ( useHyperPrior)
	//~ {
		//~ for (int i = 0; i < hyperFractionOffset; ++i)
		//~ {
			//~ double d = x.Access(x.Raw,x.Hyper,x.Position,i);
			//~ F_dF p = Normal(d,0,1);
			//~ L += p.F / effectiveBatches;
			//~ x.Assign(x.Raw,x.Hyper,x.Gradient,i, p.dF / effectiveBatches);
		//~ }
	//~ }		
	return L;
}

double LogLikelihoodPrior::TransformPrior(EfficiencyVector&x, int effectiveBatches)
{
	double L = 0;
	for (int i = 0; i < Nt; ++i)
	{
		if (BufferedGapList[i])
		{
			double d = x.Access(x.Transformed,x.Temporal,x.Position,i);
			F_dF p = GapEnforcer(d);
			L += p.F / effectiveBatches;
			x.Increment(x.Transformed,x.Temporal,x.Gradient,i,p.dF / effectiveBatches);
		}
	}
	return L;
	
}


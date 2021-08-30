#include "LogLikelihoodPrior.h"


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
		if (GapList[i])
		{
			double d = x.Access(x.Transformed,x.Temporal,x.Position,i);
			F_dF p = TemporalBetaPrior(d,gapPriorAlpha,gapPriorBeta);
			L += p.F / effectiveBatches;
			x.Increment(x.Transformed,x.Temporal,x.Gradient,i,p.dF / effectiveBatches);
		}
	}
	return L;
	
}


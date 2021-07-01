#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers


void DescentFunctor::ResetPosition()
{
	Value = 0;
	std::fill(TransformedPosition.begin(), TransformedPosition.end(),xmPrior);
	std::fill(TransformedGradient.begin(), TransformedGradient.end(),0);
	std::fill(Gradient.begin(), Gradient.end(),0);
}

void DescentFunctor::SavePosition(bool finalSave,int saveStep,bool uniqueSave,const VectorXd  & x)
{
	std::string transBase = OutputDir + "/";
	std::string intBase = OutputDir + "/";
	if (finalSave)
	{
		transBase += "FinalPosition_";
		intBase = transBase;
		ForwardTransform(x);
	}
	else
	{
		intBase += TempDirName + "/TempPosition";
		transBase = intBase;
		if (uniqueSave)
		{
			transBase += std::to_string(saveStep);
		}
		intBase += "_";
		transBase += "_";
	}
	
	
	std::fstream rawfile;
	rawfile.open(intBase + "InternalParameters.dat",std::ios::out);
	
	for (int i = 0; i < x.size(); ++i)
	{
		rawfile << std::setprecision(10) << x[i] << "\n";
	}
	std::fstream transfile;
	transfile.open(transBase + "TransformedParameters.dat",std::ios::out);

	
	for (int i = 0; i < totalTransformedParams; ++i)
	{
		transfile << std::setprecision(10) << TransformedPosition[i] << "\n";
		

	}

	rawfile.close();
	transfile.close();
}



void DescentFunctor::ForwardTransform(const VectorXd &z)
{
	ResetPosition();
	//check that cholesky decomposition has occurred, if not execute it now
	if (L.Kg_decomposed == false)
	{
		L.MakeCovarianceMatrix();
	}
	
	ForwardTransform_Temporal(z);
	
	ForwardTransform_Spatial(z);
	
	ForwardTransform_Hyper(z);
}

void DescentFunctor::ForwardTransform_Spatial(const VectorXd &z)
{
	if (SpaceActive)
	{
		
		int spaceOffset = Nt;
		if (TimeActive == false)
		{
			spaceOffset = 0;
		}
		
		std::fill(bVector.begin(), bVector.end(),0);
		for (int s = 0; s < Ns; ++s)
		{
			for (int i = 0; i < L.choleskyN; ++i)
			{
				bVector[s*Nm+L.cholesky_u[i]] += L.cholesky_w[i] * z[spaceOffset+s*Nm+L.cholesky_v[i]];
			}
		}
	
		// yml
		for (int i = 0; i < needletN; ++i)
		{
			for (int m = 0; m < Nm; ++m)
			{
				TransformedPosition[Nt+needlet_u[i]*Nm+m] += needlet_w[i]*bVector[needlet_v[i]*Nm+m];
			}
		}
	}
	else
	{
		for (int i = 0; i < Nm*Nl; ++i)
		{
			TransformedPosition[Nt + i] = FrozenSpace[i];
		}
	}
}

void DescentFunctor::ForwardTransform_Temporal(const VectorXd &z)
{
	if (TimeActive)
	{
		double u = exp(-1.0/lt);
		double ua = sqrt(1.0-u*u);
		double previous = z[Nt-1]; // First case is trivial
		TransformedPosition[Nt-1] = mut_gaps[Nt-1] + sigmat * previous;
		for (int i = Nt - 2; i >= 0; i--) 
		{
	    	previous = ua * z[i] + u * previous;
	    	TransformedPosition[i] = mut_gaps[i] + sigmat * previous;
		}
	}
	else
	{
		for (int i = 0; i < Nt; ++i)
		{
			TransformedPosition[i] = FrozenTime[i];
		}
	}
}

void DescentFunctor::ForwardTransform_Hyper(const VectorXd &z)
{
	// [[ Nt + Nl*Nm + (zeroth order weightings) + (first order weightings) + ... +(pop fractions) + popSum ]
	if (HyperActive)
	{
		int hyperStart = 0;
		if (TimeActive)
		{
			hyperStart+= Nt;
		}
		if (SpaceActive)
		{
			hyperStart += Nm*Ns;
		}
		
		int offset = hyperFractionOffset;
		double X = VerySmallLog;
		
		for (int i = 0; i < NVariancePops; ++i)
		{
			X = log_add_exp(X,z[hyperStart + offset + i]);
			
		}
		
		for (int i = 0; i < NHyper; ++i)
		{
			double normalisation = 0;
			if (i >= offset)
			{
				normalisation = X;
			}
			double v = exp(z[hyperStart + i] - normalisation);
			TransformedPosition[transformedNonHyperParams + i] = v;
		}	
	}
	else
	{
		for (int i = 0; i < NHyper; ++i)
		{
			TransformedPosition[transformedNonHyperParams + i] = FrozenHypers[i];
		}
	}
}

void DescentFunctor::BackwardTransform()
{
	BackwardTransform_Temporal();
	
	BackwardTransform_Spatial();
	
	BackwardTransform_Hyper();

}

void DescentFunctor::BackwardTransform_Temporal()
{
	if (TimeActive)
	{
		// Backward transformation
		double u = exp(-1.0/lt);
	    
		double ua = 1.0/sqrt(1.0-u*u);
		double ub = -u*ua;
		double sigmata = sigmat/ua;
	
		if (Nt > 1)
		{
		    Gradient[0] = sigmata * TransformedGradient[0];
		    for (int i = 1; i < Nt-1; i++) {
		        Gradient[i] = u * Gradient[i-1] + sigmata * TransformedGradient[i];
		    }
		    Gradient[Nt-1] = -ub * Gradient[Nt-2] + sigmat * TransformedGradient[Nt-1];
		    
		 
		}
		else
		{
			Gradient[0] = TransformedGradient[0]*sigmat;
		}
	}
}

void DescentFunctor::BackwardTransform_Spatial()
{
	if (SpaceActive)
	{
		int insertOffset = 0;
		if (TimeActive)
		{
			insertOffset += Nt;
		}
		// yml
		std::fill(bVector.begin(), bVector.end(),0);
		for (int i = 0; i < needletN; ++i)
		{
			for (int m = 0; m < Nm; ++m)
			{
				bVector[needlet_v[i]*Nm+m] += needlet_w[i]*TransformedGradient[Nt+needlet_u[i]*Nm+m];
			}
		}
		// bms = Lmnzns
		for (int s = 0; s < Ns; ++s)
		{
			for (int i = 0; i < L.choleskyN; ++i)
			{
				Gradient[insertOffset+s*Nm+L.cholesky_v[i]] += L.cholesky_w[i] * bVector[s*Nm+L.cholesky_u[i]];
			}
		}
	}
}

void DescentFunctor::BackwardTransform_Hyper()
{
	if (HyperActive)
	{
		int insertOffset = 0;
		if (TimeActive)
		{
			insertOffset+=Nt;
		}
		if (SpaceActive)
		{
			insertOffset += Ns*Nm;
		}
		
		for (int i = 0; i < hyperFractionOffset; ++i)
		{
			double x = TransformedPosition[transformedNonHyperParams + i];
			double df = TransformedGradient[transformedNonHyperParams + i];
	
			Gradient[insertOffset + i] = x* df;		
		}
		
		for (int i = 0; i < NVariancePops; ++i)
		{
			double xi = TransformedPosition[transformedNonHyperParams + hyperFractionOffset+ i];
			double sum = 0;
			
			for (int j = 0; j <  NVariancePops; ++j)
			{
				double xj = TransformedPosition[transformedNonHyperParams + hyperFractionOffset+ j];
				double dfj = TransformedGradient[transformedNonHyperParams + hyperFractionOffset+ j];
				
				double ijTerm = 0;
				if (i==j)
				{
					ijTerm = 1;
				}
				sum += dfj * xi*(ijTerm -xj);
			}
			Gradient[insertOffset + hyperFractionOffset + i] = sum;
		}
		
	}
	
}

void DescentFunctor::DistributeCalculations(const VectorXd &RawPosition, int batchID, int effectiveBatches)
{
	const int n =  totalTransformedParams;
	
	//circuitBreaker signal to workers, telling them to initiate another loop
	int circuitBreaker = batchID;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	MPI_Bcast(&effectiveBatches, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	
	//Transform then broadcast the vector to workers
	ForwardTransform(RawPosition);
	
	MPI_Bcast(&TransformedPosition[0], n, MPI_DOUBLE, RunningID, MPI_COMM_WORLD);
	
	
	L.Calculate(TransformedPosition,batchID,effectiveBatches,MaxBatches);
	
	//collect values
	double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	double Lsum = 0;
	int stars = L.StarsUsed;
	int totalStarsUsed = 0;
	MPI_Reduce(&stars, &totalStarsUsed, 1,MPI_INT, MPI_SUM, RunningID,MPI_COMM_WORLD);

	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &TransformedGradient[0], n,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	totalStarsUsed = std::max(1,totalStarsUsed);
	
	
	L.TransformPrior(TransformedPosition,&Lsum,TransformedGradient, effectiveBatches,SpaceActive,TimeActive,HyperActive);
	
	BackwardTransform();
	
	L.RawPrior(RawPosition,&Lsum,&Gradient,effectiveBatches,SpaceActive,TimeActive,HyperActive);
	
	Value = Lsum;

	StarsInLastBatch = totalStarsUsed;

	++LoopID;
	
	//negative sign for maximisation problem + normalise to number of stars
	
	for (int i = 0; i < Gradient.size(); ++i)
	{
		Gradient[i] = -Gradient[i] / StarsInLastBatch;
	}
	
	Value = -Value/StarsInLastBatch;
}

void DescentFunctor::Unfreeze()
{
	freezeOuts = std::vector<bool>(Nt,false);
}

void DescentFunctor::Calculate(const VectorXd & x)
{	
	Calculate(x,0,1);
}

void DescentFunctor::Calculate(const VectorXd &x, int batchID, int effectiveBatches)
{
	DistributeCalculations(x,batchID,effectiveBatches);
}

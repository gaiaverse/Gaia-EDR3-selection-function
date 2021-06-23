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
	
	for (int i = 0; i < totalRawParams; ++i)
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
	
	// Forward transformation
	double u = exp(-1.0/lt);
	double ua = sqrt(1.0-u*u);
	double previous = z[Nt-1]; // First case is trivial
	TransformedPosition[Nt-1] = mut_gaps[Nt-1] + sigmat * previous;
	for (int i = Nt - 2; i >= 0; i--) 
	{
    	previous = ua * z[i] + u * previous;
    	TransformedPosition[i] = mut_gaps[i] + sigmat * previous;
	}
	
	//~ std::cout << "Finished new forward bit" << std::endl;
	// bms = Lmnzns
	
	std::fill(bVector.begin(), bVector.end(),0);
	for (int s = 0; s < Ns; ++s)
	{
		for (int i = 0; i < L.choleskyN; ++i)
		{
			bVector[s*Nm+L.cholesky_u[i]] += L.cholesky_w[i] * z[Nt+s*Nm+L.cholesky_v[i]];
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
	
	//hyper params
	
	
	// [[ Nt + Nl*Nm + (zeroth order weightings) + (first order weightings) + ... +(pop fractions) + popSum ]
	

	double popSum = 0;
	int offset = hyperFractionOffset;
	for (int i = 0; i < NHyper; ++i)
	{

		double v = exp(z[rawNonHyperParams + i]);
		TransformedPosition[transformedNonHyperParams + i] = v;
		
		if (i >=offset)
		{
			popSum += v;
		}
	}
	
	for (int i = 0; i < NVariancePops; ++i)
	{
		TransformedPosition[transformedNonHyperParams + offset + i] /= popSum;
	}

}

void DescentFunctor::BackwardTransform()
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
			Gradient[Nt+s*Nm+L.cholesky_v[i]] += L.cholesky_w[i] * bVector[s*Nm+L.cholesky_u[i]];
		}
	}
	
	for (int i = 0; i < NHyper; ++i)
	{
		double x = TransformedPosition[transformedNonHyperParams + i];
		double df = TransformedGradient[transformedNonHyperParams + i];
		Gradient[rawNonHyperParams + i] = x * df;
		
		if (i >= hyperFractionOffset)
		{
			Gradient[rawNonHyperParams + i] *= (1.0 - x);
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
	
	
	
	for (int i = 0; i < Nt;++i)
	{
		if (freezeOuts[i])
		{
			TransformedGradient[i] = 0;
		}
	}
	
	
	
	BackwardTransform();
	
	L.Prior(RawPosition,&Lsum,&Gradient,effectiveBatches);
	
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
	//~ freezeOuts_mag = std::vector<bool>(Nt_m,false);
}

void DescentFunctor::Calculate(const VectorXd & x)
{	
	Calculate(x,0,1);
}

void DescentFunctor::Calculate(const VectorXd &x, int batchID, int effectiveBatches)
{
	
	DistributeCalculations(x,batchID,effectiveBatches);
	PrevLock = x;	
}

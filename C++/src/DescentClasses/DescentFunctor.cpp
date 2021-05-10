#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers

//~ void checkNan(const std::vector<double> & y, std::string origin)
//~ {
	//~ if(y.hasNaN())
	//~ {
		//~ ERROR(2,"One or more values sourced from " + origin +  "  was a NaN.");
	//~ }	
//~ }

void DescentFunctor::ResetPosition()
{
	Value = 0;
	std::fill(TransformedPosition.begin(), TransformedPosition.end(),0);
	std::fill(TransformedGradient.begin(), TransformedGradient.end(),0);
	std::fill(Gradient.begin(), Gradient.end(),0);
}

void DescentFunctor::SavePosition(bool finalSave)
{
	std::string fileBase = OutputDir + "/";
	if (finalSave)
	{
		fileBase += "FinalPosition_";
		ForwardTransform(PrevLock);
	}
	else
	{
		fileBase += TempDirName + "/TempPosition";
		if (SaveAllTemps)
		{
			fileBase += std::to_string(LoopID);
		}
		fileBase += "_";
	}
	
	
	std::fstream rawfile;
	rawfile.open(fileBase + "InternalParameters.dat",std::ios::out);
	
	for (int i = 0; i < totalRawParams; ++i)
	{
		rawfile << PrevLock[i] << "\n";
	}
	std::fstream transfile;
	transfile.open(fileBase + "TransformedParameters.dat",std::ios::out);

	
	for (int i = 0; i < totalTransformedParams; ++i)
	{
		transfile << TransformedPosition[i] << "\n";
		

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
	double ua = 1.0/sqrt(1.0-u*u);
	double ub = -u*ua;
	double previous = z[Nt-1]; // First case is trivial
	TransformedPosition[Nt-1] = mut + sigmat * previous;
	for (int i = Nt - 2; i >= 0; i--) 
	{
    	previous = (z[i] - ub * previous) / ua;
    	TransformedPosition[i] = mut + sigmat * previous;
	}

	// bms = Lmnzns
	for (int s = 0; s < Ns; ++s)
	{
		for (int m = 0; m < Nm; ++m)
		{
			bVector[s*Nm+m] = 0;
			for (int n = 0; n <= m; ++n)
			{
				bVector[s*Nm+m] += L.CholeskyKg(m,n) * z[Nt+s*Nm+n];
			}
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
		for (int m = 0; m < Nm; ++m)
		{
			for (int n = 0; n <=m; ++n)
			{
				Gradient[Nt+s*Nm+n] += L.CholeskyKg(m,n)*bVector[s*Nm+m];
			}
		}
	}

}

void DescentFunctor::DistributeCalculations(const VectorXd &RawPosition, int batchID, int effectiveBatches)
{
	const int n =  Nt+Nm*Nl;
	
	//circuitBreaker signal to workers, telling them to initiate another loop
	int circuitBreaker = batchID;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	MPI_Bcast(&effectiveBatches, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	
	//Transform then broadcast the vector to workers
	ForwardTransform(RawPosition);
	MPI_Bcast(&TransformedPosition[0], n, MPI_DOUBLE, RunningID, MPI_COMM_WORLD);
	L.Calculate(TransformedPosition,batchID,effectiveBatches);
	
	
	//collect values
	const double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	double Lsum = 0;
	int stars = L.StarsUsed;
	int totalStarsUsed = 0;
	
	MPI_Reduce(&stars, &totalStarsUsed, 1,MPI_INT, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &TransformedGradient[0], n,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	
	
	BackwardTransform();
	L.Prior(RawPosition,&Lsum,&Gradient,effectiveBatches);
	Value = Lsum;

	StarsInLastBatch = totalStarsUsed;

	//~ checkNan(Gradient,"Gradient Calculation");
	++LoopID;
	
	//negative sign for maximisation problem + normalise to number of stars
	for (int i = 0; i < Gradient.size(); ++i)
	{
		Gradient[i] = -Gradient[i]/StarsInLastBatch;
	}
	Value = -Value/StarsInLastBatch;
}

void DescentFunctor::Calculate(const VectorXd & x)
{
	if (N_SGD_Batches != 1)
	{
		ERROR(-9, "You have called somehow ended up calling a single-batch function when your data is distributed into multiple batches - you should go and cry now");
	}
	
	Calculate(x,0,1);
}

void DescentFunctor::Calculate(const VectorXd &x, int batchID, int effectiveBatches)
{
	
	DistributeCalculations(x,batchID,effectiveBatches);
	PrevLock = x;	
}

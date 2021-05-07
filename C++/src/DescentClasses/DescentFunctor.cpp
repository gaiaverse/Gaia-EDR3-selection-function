#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers

void checkNan(const VectorXd & y, std::string origin)
{
	if(y.hasNaN())
	{
		ERROR(2,"One or more values sourced from " + origin +  "  was a NaN.");
	}	
}

void DescentFunctor::ResetPosition()
{
	for (int i =0; i < totalTransformedParams; ++i)
	{
		TransformedPosition[i] = 0;
		TransformedGradient[i] = 0;
	}
	
	for (int i = 0; i < totalRawParams; ++i)
	{
		CurrentGradient[i] = 0;
	}
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
	
	
	VectorXd b = VectorXd::Zero(Nm*Ns);
	
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
			for (int n = 0; n < Nm; ++n)
			{
				b[s*Nm+m] += L.CholeskyKg(m,n) * z[Nt+s*Nm+n];
			}
		}
	}

	// yml
	for (int i = 0; i < needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			TransformedPosition[Nt+needlet_u[i]*Nm+m] += needlet_w[i]*b[needlet_v[i]*Nm+m];
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
	VectorXd bGradient = VectorXd::Zero(Ns*Nm);

	if (Nt > 1)
	{
	    CurrentGradient[0] = sigmata * TransformedGradient[0];
	    for (int i = 1; i < Nt-1; i++) {
	        CurrentGradient[i] = u * CurrentGradient[i-1] + sigmata * TransformedGradient[i];
	    }
	    CurrentGradient[Nt-1] = -ub * CurrentGradient[Nt-2] + sigmat * TransformedGradient[Nt-1];
	}
	else
	{
		CurrentGradient[0] = TransformedGradient[0]*sigmat;
	}
	// yml
	for (int i = 0; i < needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			bGradient[needlet_v[i]*Nm+m] += needlet_w[i]*TransformedGradient[Nt+needlet_u[i]*Nm+m];
		}
	}

	// bms = Lmnzns
	for (int s = 0; s < Ns; ++s)
	{
		for (int m = 0; m < Nm; ++m)
		{
			for (int n = 0; n < Nm; ++n)
			{
				
				CurrentGradient[Nt+s*Nm+n] += L.CholeskyKg(m,n)*bGradient[s*Nm+m];
			}
		}
	}

}


void DescentFunctor::DistributeCalculations(const VectorXd &RawPosition)
{
	const int n =  Nt+Nm*Nl;
	
	//circuitBreaker signal to workers, telling them to initiate another loop
	int circuitBreaker = 1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);

	//Transform then broadcast the vector to workers
	ForwardTransform(RawPosition);
	MPI_Bcast(&TransformedPosition[0], n, MPI_DOUBLE, RunningID, MPI_COMM_WORLD);
	L.Calculate(TransformedPosition);
	
	
	//collect values
	double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	double Lsum = 0;			
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &TransformedGradient[0], n,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);

	
	BackwardTransform();
	L.Prior(RawPosition,&Lsum,&CurrentGradient);
	CurrentValue = Lsum;

	
	checkNan(CurrentGradient,"Gradient Calculation");
	++LoopID;
}

void DescentFunctor::Calculate(const VectorXd &x)
{
	
	DistributeCalculations(x);
	PrevLock = x;
	
	//negative sign for maximisation problem + normalise to number of stars
	for (int i = 0; i < x.size(); ++i)
	{
		Gradient[i] = -CurrentGradient[i]/NStars;
	}
	Value = -CurrentValue/NStars;
	
}
double DescentFunctor::value(const VectorXd &y)
{
	VectorXd diff = (y - PrevLock);
	checkNan(y," Position (value call) ");
	double key = diff.norm();

	
	if (key > lockLim)
	{
		DistributeCalculations(y);
		PrevLock = y;
	}

	return -CurrentValue;
}
void DescentFunctor::gradient(const VectorXd &y, VectorXd &grad)
{
	checkNan(y," Position (grad call) ");
	VectorXd diff = (y - PrevLock);
	double key = diff.norm();

	if (key > lockLim)
	{
		DistributeCalculations(y);
		PrevLock = y;
	}

	//negative sign for maximisation problem
	for (int i = 0; i < y.size(); ++i)
	{
		grad[i] = -CurrentGradient[i];
	}
}


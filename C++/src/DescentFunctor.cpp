#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers

void checkNan(const VectorXd & y, std::string origin)
{
	if(y.hasNaN())
	{
		std::cout << "\n\n\nERROR \n One or more values sourced from " << origin << "  was a NaN. Cannot handle. Goodbye";
		exit(1); 
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
	if (finalSave)
	{
		std::cout << "\tConverged roots : ";
	}
	for (int i = 0; i < totalTransformedParams; ++i)
	{
		transfile << TransformedPosition[i] << "\n";
		if (finalSave)
		{
			std::cout << TransformedPosition[i] << ",\t";
		}
	}
	std::cout << "\n";
	
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
	for (int i = 0; i < L.needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			TransformedPosition[Nt+L.needlet_u[i]*Nm+m] += L.needlet_w[i]*b[L.needlet_v[i]*Nm+m];
		}
	}
}

void DescentFunctor::BackwardTransform(bool print)
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
	for (int i = 0; i < L.needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			bGradient[L.needlet_v[i]*Nm+m] += L.needlet_w[i]*TransformedGradient[Nt+L.needlet_u[i]*Nm+m];
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


void DescentFunctor::TestGradient(const VectorXd y)
{
	VectorXd ySave = y;
	VectorXd oldGrad = CurrentGradient;
	double oldVal = CurrentValue;
	

	double delta = 1e-5;
	double trueL = CurrentValue;
	
	VectorXd gradOverride = CurrentGradient;
	VectorXd yNudge;
	for (int i = 0; i < y.size(); ++i)
	{
		yNudge = y;
		double shove = std::max(delta*delta,delta*y[i]);
		yNudge[i] += shove;

		
		DistributeCalculations(yNudge,false);
		
		double testGrad = (CurrentValue - trueL)/(shove);
		//std::cout << "\t\t\t\tManual Gradient: " << i << ":  " << testGrad << ", \t Exact Gradient: " << oldGrad[i] <<  "\n";
		gradOverride[i] = testGrad;
	}
	
	CurrentValue = oldVal;
	CurrentGradient = oldGrad;
	bool overrideOn = true;
	bool unsafe = std::isnan(CurrentValue);
	if (overrideOn & !(unsafe))
	{
		CurrentGradient = gradOverride;
		std::cout << "\t\t\t\t\tGradient by manual, new value: " << CurrentGradient.transpose() << std::endl;
	}
}


void DescentFunctor::DistributeCalculations(const TVector &y, bool printOn)
{
	
	if (printOn)
	{
		
		
		//~ std::cout << "\n\t\tEntering calculation iteration " << LoopID<< ". "; printTime();
		//~ std::cout << "Current position: " << y.transpose() << std::endl;
		++LoopID;
	}
	const int n =  Nt+Nm*Nl;
	
	//std::cout << "\tCalculation distribution " << LoopID << " begun" << std::endl;
	const VectorXd RawPosition = y;

	
	//circuitBreaker signal to workers, telling them to initiate another loop
	int circuitBreaker = 1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);



	ForwardTransform(RawPosition);
		

	//send position vector to workers
	MPI_Bcast(&TransformedPosition[0], n, MPI_DOUBLE, RunningID, MPI_COMM_WORLD);

	
	// IMPORTANT!
	// The Root also executes their own calculations of L on their provided data 
	L.Calculate(TransformedPosition);
	double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	
	//collect values
	double Lsum = 0;
	
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &TransformedGradient[0], n,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);

	
	if (LoopID % SaveSteps == 0)
	{
		SavePosition(false);
	}
	
	
	BackwardTransform(printOn);
	
	
	L.Prior(RawPosition,&Lsum,&CurrentGradient);
	
	CurrentValue = Lsum;
	
	//~ if (printOn)
	//~ {
		//~ std::cout << "\t\t\tCurrent position: " << y.transpose() << "\n \t\t\tCurrent Gradient: " << CurrentGradient.transpose() << std::endl;
		//~ std::cout << "\n\t\t\tTransformed position: " << TransformedPosition.transpose() << "\n";
		//~ std::cout << "\t\t\tTransformed Gradient: " << TransformedGradient.transpose() << std::endl;
		//~ std::cout << "\t\t\tL = " << CurrentValue << "   Gradnorm = " <<  CurrentGradient.norm() << std::endl;
		
	//~ }
	checkNan(CurrentGradient,"Gradient Calculation");
}


double DescentFunctor::value(const TVector &y)
{
	VectorXd diff = (y - PrevLock);
	checkNan(y," Position (value call) ");
	double key = diff.norm();

	
	if (key > lockLim)
	{
		DistributeCalculations(y,true);
		PrevLock = y;
		
		//TestGradient(y);
	}


	return -CurrentValue;
}
void DescentFunctor::gradient(const TVector &y, TVector &grad)
{
	checkNan(y," Position (grad call) ");
	VectorXd diff = (y - PrevLock);
	double key = diff.norm();

	if (key > lockLim)
	{
		DistributeCalculations(y,true);
		PrevLock = y;
		
		
		//TestGradient(y);
	}

	//negative sign for maximisation problem
	for (int i = 0; i < y.size(); ++i)
	{
		grad[i] = -CurrentGradient[i];
	}

}


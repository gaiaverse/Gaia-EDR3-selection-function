#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers


void DescentFunctor::ResetPosition()
{
	int n = Nt+Nm*Nl;
	
	for (int i =0; i < n; ++i)
	{
		TransformedPosition[n] = 0;
		TransformedGradient[n] = 0;
	}
	
	for (int i = 0; i < CurrentGradient.size(); ++i)
	{
		CurrentGradient[i] = 0;
	}
	
	
}
void DescentFunctor::ForwardTransform(VectorXd &z)
{
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
	for (int i = Nt - 2; i >= 0; i--) {
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
				//b[s*Nm+m] += L.CholeskyKg[m,n] * z[Nt+s*Nm+n];
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

void DescentFunctor::BackwardTransform()
{
	// Backward transformation
	double u = exp(-1.0/lt);
	double ua = 1.0/sqrt(1.0-u*u);
	double ub = -u*ua;
	VectorXd bGradient = VectorXd::Zero(Ns*Nm);
	CurrentGradient[Nt-1] = TransformedGradient[Nt-1]/sigmat;
	for (int i = Nt - 2; i >= 0; i--) {
		CurrentGradient[i] = (ua * TransformedGradient[i] + ub * TransformedGradient[i+1]) / sigmat;
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
				//CurrentGradient[Nt+s*Nm+n] += L.CholeskyKg[m,n]*bGradient[s*Nm+m];
			}
		}
	}

	
}

void DescentFunctor::DistributeCalculations(const TVector &y)
{
	std::cout << "\t\tDescent functor has been activated" << std::endl;
	const int n =  Nt+Nm*Nl;
	ResetPosition();
	//std::cout << "\tCalculation distribution " << LoopID << " begun" << std::endl;
	VectorXd RawPosition = y;

	
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

	++LoopID;
	
	
	BackwardTransform();
	
	L.Prior(RawPosition,&CurrentValue,&CurrentGradient);
}


void DescentFunctor::ExamineInterestVectors(Eigen::VectorXd& position)
{
	//use this vector to pluck variables of interest into the comparison vector
	std::vector<int> interestIDs = {0,1,2,3,4};
	
	std::vector<std::string> interestNames = {"log_lt", "log_lg", "log_sigma2", "log_m","log_tau2"};
	if (LoopID == 0)
	{
		std::cout << "\nInitial guess parameters: \n\t";
		InterestVectors.resize(interestIDs.size());
		for (int i = 0; i < interestIDs.size(); ++i)
		{
			InterestVectors[i] = position[i];
			std::cout << interestNames[i] << " = " << position[i] << "   ";
		}
		std::cout << std::endl;
	}
	else
	{
		auto checkpoint = std::chrono::system_clock::now();
		std::string duration = formatDuration(Start,checkpoint);
		std::cout << "\n\n Loop " << LoopID << ": updating parameters. Current runtime: " + duration;
		std::cout << " \n\t";
		for (int i = 0; i < interestIDs.size(); ++i)
		{
			double delta = (position[i] - InterestVectors[i])/InterestVectors[i]*100;
			InterestVectors[i] = position[i];
			std::string val;
			if (!std::isnan(delta) && !std::isinf(delta))
			{
				val = std::to_string(delta);
				if (delta <  0.0)
				{
					val = val;
				}
				else
				{
					val = "+" + val;
				}
				
				val += "%";
			}
			else
			{
				val = "NaN";
			}
			std::cout << interestNames[i] << " = " << InterestVectors[i] << "(" << val << ")";
		}
		std::cout << "\n";
	}
}



double DescentFunctor::value(const TVector &y)
{
	VectorXd diff = (y - PrevLock);
	double key = diff.norm();

	std::cout.precision(10);
	if (key > lockLim)
	{

		DistributeCalculations(y);
		PrevLock = y;
	}

	return -CurrentValue;
}
void DescentFunctor::gradient(const TVector &y, TVector &grad)
{

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

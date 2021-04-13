#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers


void DescentFunctor::DistributeCalculations(const TVector &y)
{
	std::cout << "\t\tDescent functor has been activated" << std::endl;
	
	//std::cout << "\tCalculation distribution " << LoopID << " begun" << std::endl;
	VectorXd z = y;
	x = VectorXd::Zero(Nt+Nm*Nl);
	b = VectorXd::Zero(Nm*Ns);
	// Somehow need a Cholesky matrix here

	//std::cout << "Position:\n" << x.transpose() << std::endl;
	//ExamineInterestVectors(x);
	//circuitBreaker signal to workers, telling them to initiate another loop
	int n = x.size();
	int circuitBreaker = 1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);

	// Forward transformation
	double u = exp(-1.0/lt);
	double ua = 1.0/sqrt(1.0-u*u);
	double ub = -u*ua;
	double previous = z[Nt-1]; // First case is trivial
	x[Nt-1] = mut + sigmat * previous;
	for (int i = Nt - 2; i >= 0; i--) {
    	previous = (z[i] - ub * previous) / ua;
    	x[i] = mut + sigmat * previous;
	}

	// bms = Lmnzns
	for (int s = 0; s < Ns; ++s)
	{
		for (int m = 0; m < Nm; ++m)
		{
			for (int n = 0; n < Nm; ++n)
			{
				b[s*Nm+m] += L[m,n] * z[Nt+s*Nm+n];
			}
		}
	}

	// yml
	for (int i = 0; i < needlet_n; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			x[Nt+needlet_u[i]*Nm+m] += needlet_w[i]*b[needlet_v[i]*Nm+m];
		}
	}

	// x[0:Nm] = mu, x[Nm:Nm+Nt] = zt, x[Nm+Nt:] = zml
	
	
	//send position vector to workers
	MPI_Bcast(&x[0], n, MPI_DOUBLE, RunningID, MPI_COMM_WORLD);

	//initialise an empty array for summing in the output results of GradL and L
	std::vector<double> empty(n,0.0);
	std::vector<double> output(n,0.0);
	double emptyS = 0;
	
	
	// IMPORTANT!
	// The Root also executes their own calculations of L on their provided data 
	L.Calculate(x);
	double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	
	//collect values
	RawGradient = VectorXd::Zero(n);
	double Lsum = 0;
	
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &RawGradient[0], n,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);

	++LoopID;
	
	

	CurrentValue = Lsum;
	CurrentGradient = VectorXd::Zero(Nt+Ns*Nm);
	bGradient = VectorXd::Zero(Ns*Nm);


	// Backward transformation

	CurrentGradient[Nt-1] = RawGradient[Nt-1]/sigmat;
	for (int i = Nt - 2; i >= 0; i--) {
		CurrentGradient[i] = (a * RawGradient[i] + b * RawGradient[i+1]) / sigmat;
	}

	// yml
	for (int i = 0; i < needlet_n; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			bGradient[needlet_v[i]*Nm+m] += needlet_w[i]*RawGradient[Nt+needlet_u[i]*Nm+m];
		}
	}

	// bms = Lmnzns
	for (int s = 0; s < Ns; ++s)
	{
		for (int m = 0; m < Nm; ++m)
		{
			for (int n = 0; n < Nm; ++n)
			{
				CurrentGradient[Nt+s*Nm+n] += L[m,n]*bGradient[s*Nm+m];
			}
		}
	}



	//std::cout << "Gradient:\n" << CurrentGradient.transpose() << std::endl;
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

#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers

double DescentFunctor::operator()(Eigen::VectorXd& x, Eigen::VectorXd& grad)
{
	
	//circuitBreaker signal to workers, telling them to initiate another loop
	int n = x.size();
	int circuitBreaker = 1;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	
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
	double Lsum = 0;
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &grad[0], n,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);

	return Lsum;
}

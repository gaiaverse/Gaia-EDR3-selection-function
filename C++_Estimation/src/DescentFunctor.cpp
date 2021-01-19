#include "DescentFunctor.h"

double DescentFunctor::operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad)
{

	//chance that this copies it + wastes a bunch of time...
	std::vector<double> xVec;
	xVec.resize(x.size());
	VectorXd::Map(&xVec[0], x.size()) = x;
	
	//send vector size so can be computed
	int n = xVec.size();
	MPI_Bcast(&n, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	
	//send position vector to workers
	MPI_Bcast(&xVec[0], x.size(), MPI_DOUBLE, RunningID, MPI_COMM_WORLD);

	
	std::vector<double> empty(n,0.0);
	std::vector<double> output(n,0.0);
	double emptyS = 0;
	
	L.Calculate(xVec);
	double l = L.Value;
	double Lsum = 0;
	MPI_Barrier(MPI_COMM_WORLD);
	
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &output[0], n,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);

	

	
	for (int i = 0; i < n; ++i)
	{
		
		grad[i] = output[i];

	}

	return Lsum;
}

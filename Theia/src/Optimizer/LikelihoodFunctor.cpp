#include "LikelihoodFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers


void LikelihoodFunctor::SavePosition(bool finalSave,int saveStep,bool uniqueSave)
{
	Efficiency.Save(finalSave,saveStep,uniqueSave);
}

void LikelihoodFunctor::Calculate(const std::vector<double> & x)
{	
	Calculate(x,0,1);
}

void LikelihoodFunctor::Calculate(const std::vector<double> &x, int batchID, int effectiveBatches)
{
	const int n =  totalTransformedParams;
	
	//circuitBreaker signal to workers, telling them to initiate another loop

	int breaker = batchID;
	MPI_Bcast(&breaker, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	MPI_Bcast(&effectiveBatches, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	//Transform then broadcast the vector to workers
	Efficiency.ForwardTransform(x);
	MPI_Bcast(&Efficiency.TransformedPosition[0], n, MPI_DOUBLE, RootID, MPI_COMM_WORLD);


	L.Calculate(Efficiency,batchID,effectiveBatches,MaxBatches);

	//collect values
	double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	double Lsum = 0;
	int stars = L.StarsUsed;
	int totalStarsUsed = 0;
	MPI_Reduce(&stars, &totalStarsUsed, 1,MPI_INT, MPI_SUM, RootID,MPI_COMM_WORLD);

	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	
	std::vector<double> gradientCatcher(totalTransformedParams,0.0);
	MPI_Reduce(&L.Gradient[0], &gradientCatcher[0], n,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	

	totalStarsUsed = std::max(1,totalStarsUsed);
	
	//Efficiency.TransformedGradient= gradientCatcher;
	
	
	Lsum += L.TransformPrior(Efficiency, effectiveBatches);

	Efficiency.BackwardTransform();
	
	
	//Lsum += L.RawPrior(Efficiency,effectiveBatches);

	Value = Lsum;

	int StarsInLastBatch = totalStarsUsed;

	++LoopID;
	
	//negative sign for maximisation problem + normalise to number of stars
	for (int i = 0; i < Efficiency.RawGradient.size(); ++i)
	{
		Gradient[i] = -Efficiency.RawGradient[i] / totalStarsUsed;
	}
	Value = -Value/totalStarsUsed;
}

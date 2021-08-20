#include "DescentFunctor.h"

//the overloaded () operator executes the full evaluation of L and GradL at a given position x
//this function is executed only once (by root), and distributes the task to the remaining workers


void DescentFunctor::ResetPosition()
{
	Value = 0;
}

void DescentFunctor::SavePosition(bool finalSave,int saveStep,bool uniqueSave,const VectorXd  & x)
{
	std::string transBase = OutputDir + "/";
	std::string intBase = OutputDir + "/";
	Efficiency.ForwardTransform(x);
	if (finalSave)
	{
		transBase += "FinalPosition_";
		intBase = transBase;
		
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
		transfile << std::setprecision(10) << Efficiency.TransformedPosition[i] << "\n";
		

	}

	rawfile.close();
	transfile.close();
	

}


void DescentFunctor::DistributeCalculations(const VectorXd &inputPosition, int batchID, int effectiveBatches)
{
	const int n =  totalTransformedParams;
	
	//circuitBreaker signal to workers, telling them to initiate another loop
	MPI_Bcast(&batchID, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	MPI_Bcast(&effectiveBatches, 1, MPI_INT, RootID, MPI_COMM_WORLD);
	
	//Transform then broadcast the vector to workers
	Efficiency.ForwardTransform(inputPosition);
	MPI_Bcast(&Efficiency.TransformedPosition[0], n, MPI_DOUBLE, RootID, MPI_COMM_WORLD);
	
	L.Calculate(Efficiency.TransformedPosition,batchID,effectiveBatches,MaxBatches);
	
	//collect values
	double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	double Lsum = 0;
	int stars = L.StarsUsed;
	int totalStarsUsed = 0;
	MPI_Reduce(&stars, &totalStarsUsed, 1,MPI_INT, MPI_SUM, RootID,MPI_COMM_WORLD);

	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &Efficiency.TransformedGradient[0], n,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	totalStarsUsed = std::max(1,totalStarsUsed);
	
	L.TransformPrior(Efficiency.TransformedPosition,&Lsum,Efficiency.TransformedGradient, effectiveBatches);
	
	Efficiency.BackwardTransform();
	
	L.RawPrior(Efficiency.RawPosition,&Lsum,&Efficiency.RawGradient,effectiveBatches);
	
	Value = Lsum;

	StarsInLastBatch = totalStarsUsed;

	++LoopID;
	
	//negative sign for maximisation problem + normalise to number of stars
	for (int i = 0; i < Efficiency.RawGradient.size(); ++i)
	{
		Gradient[i] = -Efficiency.RawGradient[i] / StarsInLastBatch;
	}
	Value = -Value/StarsInLastBatch;
}


void DescentFunctor::Calculate(const VectorXd & x)
{	
	Calculate(x,0,1);
}

void DescentFunctor::Calculate(const VectorXd &x, int batchID, int effectiveBatches)
{
	DistributeCalculations(x,batchID,effectiveBatches);
}

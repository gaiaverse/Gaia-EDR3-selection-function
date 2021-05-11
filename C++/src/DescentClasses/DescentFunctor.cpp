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
	std::fill(Gradient.begin(), Gradient.end(),0);
}

void DescentFunctor::SavePosition(bool finalSave)
{
	std::string fileBase = OutputDir + "/";
	if (finalSave)
	{
		fileBase += "FinalPosition_";
		//ForwardTransform(PrevLock);
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
		//rawfile << PrevLock[i] << "\n";
	}
	std::fstream transfile;
	transfile.open(fileBase + "TransformedParameters.dat",std::ios::out);

	
	for (int i = 0; i < totalTransformedParams; ++i)
	{
	//	transfile << TransformedPosition[i] << "\n";
		

	}

	rawfile.close();
	transfile.close();
}



void DescentFunctor::DistributeCalculations(const VectorXd &inputPosition, int batchID, int effectiveBatches)
{
	ResetPosition();
	
	std::vector<double> RawPosition(inputPosition.data(), inputPosition.data() + inputPosition.size() );
	
	const int n =  totalRawParams;
	
	//circuitBreaker signal to workers, telling them to initiate another loop
	int circuitBreaker = batchID;
	MPI_Bcast(&circuitBreaker, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	MPI_Bcast(&effectiveBatches, 1, MPI_INT, RunningID, MPI_COMM_WORLD);
	
	//Transform then broadcast the vector to workers

	MPI_Bcast(&RawPosition[0], n, MPI_DOUBLE, RunningID, MPI_COMM_WORLD);
	L.PriorCalculate(RawPosition,batchID,effectiveBatches);
	
	
	//collect values
	const double l = L.Value; //as with the workers, have to store here temporarily for a reason I don't understand. It breaks if you MPI_Reduce(&L.Value), so learn from my mistake
	double Lsum = 0;
	int stars = L.StarsUsed;
	int totalStarsUsed = 0;
	
	MPI_Reduce(&stars, &totalStarsUsed, 1,MPI_INT, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&l, &Lsum, 1,MPI_DOUBLE, MPI_SUM, RunningID,MPI_COMM_WORLD);
	MPI_Reduce(&L.Gradient[0], &Gradient[0], n,MPI_DOUBLE, MPI_SUM, RootID,MPI_COMM_WORLD);
	

	Value = Lsum;
	StarsInLastBatch = totalStarsUsed;

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

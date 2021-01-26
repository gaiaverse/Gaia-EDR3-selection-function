#include "Liklihood.h"

Liklihood::Liklihood(const std::vector<Star> &data, int nPoints, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = std::vector<double>(nPoints,0.0);
	//NVists = 5
}


void Liklihood::Calculate(Eigen::VectorXd& x)
{

	Value = 0;
	Gradient[0] = 0;
	Gradient[1] = 0;
	Gradient[2] = 0;
	Gradient[3] = 0;	
	

	if (ID == 0)
	{
		Prior();
	}

	for (int i = 0; i < Data.size(); ++i)
	{

	}

	
}


void Liklihood::Prior()
{
	//Value += whatever
	//Gradient[i] += etc
}

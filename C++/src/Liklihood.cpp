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


void Liklihood::Prior(Eigen::VectorXd& params)
{
    
    // Unpack parameters
    double lt = exp(params(0));
    double lg = exp(params(1));
    double sigma2 = exp(params(2));
    double m = exp(params(3));
    double tau2 = exp(params(4));
    
    VectorXd mu = params.segment(Nh, Ng);
    VectorXd x = params.segment(Nh+Ng, Ng*Nt);
    
	//Value += whatever
	//Gradient[i] += etc
}

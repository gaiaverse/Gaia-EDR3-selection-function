#include "GlobalVariables.h"
Eigen::VectorXd initialisedVector(int n)
{
	 VectorXd x = VectorXd::Zero(n);
    
    
    //~ //initialisation of hyperhyperparameters
    std::vector<double> hyperhyper = {0,0,0,2,0};
    for (int i = 0; i < Nh; ++i)
    {
		x[i] = hyperhyper[i];
	}

	//initialisation of boring old hyperparameters
	for (int i = 0; i < Ng; ++i)
	{
		x[Nh + i] = sin(2 * M_PI * (double)i/5) * exp(-(double)i/30);
	}
	
	return x;
}

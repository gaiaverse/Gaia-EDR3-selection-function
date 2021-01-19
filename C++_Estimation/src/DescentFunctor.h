#pragma once
#include <vector>
#include <mpi.h>
#include "libs/Eigen/Core"
#include "customClasses.h"
#include "Liklihood.h"
using Eigen::VectorXd;



 
class DescentFunctor
{
private:
	int RunningID;
	const std::vector<Star> &Data; 
	Liklihood L;
	
public:
    DescentFunctor(int n,const std::vector<Star> & data, int nParams) : RunningID(n), Data(data), L(data,nParams,n)
    {
			
	}
    double operator()(const Eigen::VectorXd& x, Eigen::VectorXd& grad);
 
};

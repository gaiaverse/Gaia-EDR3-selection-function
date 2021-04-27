#include "GlobalVariables.h"


Eigen::VectorXd initialisedVector(int n,bool print)
{
	 VectorXd x = VectorXd::Random(n);
    
   if (print)
   {
	   std::cout << "\tInitial position set to: (";
	   for (int i = 0; i < n; ++i)
	   {
		   std::cout << x[i];
		   if (i < n-1)
		   {
			   std::cout << ", ";
		   } 
	   }
	   std::cout << ")\n";
   }
   
	return x;
}

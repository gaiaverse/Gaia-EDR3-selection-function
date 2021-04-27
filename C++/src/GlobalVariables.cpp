#include "GlobalVariables.h"


Eigen::VectorXd initialisedVector(int n)
{
	 VectorXd x = VectorXd::Random(n);
    
   
	GlobalLog(1,
	   std::cout << "\tPosition Vetcor initialised to: (";
	   for (int i = 0; i < n; ++i)
	   {
		   std::cout << x[i];
		   if (i < n-1)
		   {
			   std::cout << ", ";
		   } 
	   }
	   std::cout << ")\n";
	);
   
	return x;
}

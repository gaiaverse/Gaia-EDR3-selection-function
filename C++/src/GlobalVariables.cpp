#include "GlobalVariables.h"


Eigen::VectorXd initialisedVector(int n)
{
	 VectorXd x = initialisationBounds*VectorXd::Random(n);
    
   
	//~ GlobalLog(0,
	   //~ std::cout << "\tPosition Vector initialised to: (";
	   //~ for (int i = 0; i < n; ++i)
	   //~ {
		   //~ std::cout << x[i];
		   //~ if (i < n-1)
		   //~ {
			   //~ std::cout << ", ";
		   //~ } 
	   //~ }
	   //~ std::cout << ")\n";
	//~ );
   
	return x;
}
void PrintStatus(std::string location)
{
	std::vector<std::string> properties = {"Nt","Nm","healpix_order","needlet_order","Nl","Ns","totalRawParams","totalTransformedParams","mu_t","sigma_t","l_m","l_t","initialisationBounds"};
	std::vector<double> vals = {Nt,Nm,healpix_order, needlet_order, Nl,Ns,totalRawParams,totalTransformedParams,mut,sigmat,lm,lt,initialisationBounds};
	
	std::fstream file;
	file.open(location + "/Optimiser_Properties.dat",std::ios::out);
	int w = 15;
	
	for (int i = 0; i < 2; ++i)
	{
		for (int j = 0; j < properties.size(); ++j)
		{
			//~ file << std::setw(w) << std::left;
			if (i ==0)
			{
				file << properties[j];
			}
			else
			{
				file << vals[j];
			}
			if (j < properties.size() - 1)
			{
				file << ", ";
			}
		}
		file << "\n";
	}
	file.close();
}

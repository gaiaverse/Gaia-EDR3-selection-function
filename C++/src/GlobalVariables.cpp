#include "GlobalVariables.h"


Eigen::VectorXd initialisedVector(int n, bool loadIn, std::string loadLocation)
{
	VectorXd x;
	if (loadIn)
	{
		
		VectorXd possibleX = VectorXd::Zero(n);
		
		int i = 0;
		forLineInFile(loadLocation,
			possibleX[i] = (std::stod(FILE_LINE));
			++i;
			if (i > n)
			{
				std::cout << "Internal quit " << std::endl;
				ERROR(100,"Asked to load in start position from file, but it was the wrong length");
			}
		);
		
		if (i < n);
		{
			std::cout << "End quit, i = " << i << std::endl;
			ERROR(100,"Asked to load in start position from file, but it was the wrong length");
		}
		
		x = possibleX;
	}
	else
	{
		x = initialisationBounds*VectorXd::Random(n);
    }
   
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

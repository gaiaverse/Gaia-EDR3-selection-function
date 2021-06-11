#include "GlobalVariables.h"


Eigen::VectorXd initialisedVector(int n, bool loadIn, std::string loadLocation)
{
	VectorXd x;
	if (loadIn)
	{
		
		VectorXd possibleX = VectorXd::Zero(n);
		
		int i = 0;
		forLineIn(loadLocation,
			possibleX[i] = (std::stod(FILE_LINE));
			++i;
			if (i > n)
			{
				std::cout << "Internal quit " << std::endl;
				ERROR(100,"Asked to load in start position from file, but it was the wrong length");
			}
		);
		
		bool wentUnder = (i < n);
		if (wentUnder)
		{
			std::cout << "End quit, i = " << i << "  n = " << n << "  " << wentUnder<< std::endl;
			ERROR(100,"Asked to load in start position from file, but it was the wrong length");
		}
		
		x = possibleX;
	}
	else
	{
		x = initialisationBounds*VectorXd::Random(n);
		
		//initialise spatial part properly
		
		Eigen::Matrix<double, Nm, Nm> Kg;
		for (int i = 0; i < Nm; i++) 
		{
			for (int j = 0; j < i; j++) 
			{
				Kg(i,j) = Kg(j,i) = exp(-pow(i - j,2)/(2.0*lm*lm));
			}
			Kg(i,i) = 1.0 + SingularityPreventer;
		}
		
		//decompose to make CholeskyKg
		Eigen::Matrix<double, Nm, Nm> CholeskyKg = Kg.llt().matrixL();
		Eigen::Matrix<double, Nm, Nm>  Inverted = CholeskyKg.inverse();
		
		Eigen::VectorXd mums = VectorXd::Constant(Nm,mum_init - mum_prior);
		
		mums = Inverted * mums;
		
		for (int i = 0; i < Nm; ++i)
		{
			x[Nt+i] += mums[i];
		}
	}
   
	return x;
}
void PrintStatus(std::string location)
{
	std::vector<std::string> properties = {"Nt","Nm","healpix_order","needlet_order","Nl","Ns","totalRawParams","totalTransformedParams","mu_t","sigma_t","l_m","l_t","initialisationBounds"};
	std::vector<double> vals = {(double)Nt,(double)Nm,(double)healpix_order, (double)needlet_order, (double)Nl,(double)Ns,(double)totalRawParams,(double)totalTransformedParams,(double)mut_normal,(double)sigmat,(double)lm,(double)lt,(double)initialisationBounds};
	
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
				file << std::setprecision(10) <<vals[j];
				
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

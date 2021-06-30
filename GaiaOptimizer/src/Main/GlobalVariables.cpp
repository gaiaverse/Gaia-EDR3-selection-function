#include "GlobalVariables.h"
#include "GlobalConstants.h"

Eigen::VectorXd initialisedVector(int n, std::string loadLocation)
{
	VectorXd x;
	bool loadIn = !(loadLocation == "__null_location__");
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
		
		Eigen::VectorXd mums = VectorXd::Constant(Nm,xmInitialised - xmPrior);
		
		mums = Inverted * mums;
		
		for (int i = 0; i < Nm; ++i)
		{
			x[Nt+i] += mums[i];
		}
		
		
		//hyper parameters
		for (int i = 0; i < hyperOrder+ 1; ++i)
		{			
			double prior = pow(10.0,-(1+2*i));
			if ( i > 0)
			{
				prior = pow(prior,1.0/i);
			}
			else
			{
				prior = 10;
			}
			prior = log(prior);
			
			for (int j = 0;j < NVariancePops; ++j)
			{
				int index = rawNonHyperParams + i*NVariancePops + j;
				x[index] += prior;
			}
			
		}
		for (int i = 0; i < NVariancePops; ++i)
		{
			x[rawNonHyperParams + hyperFractionOffset + i] = pow(-0.65822,i);
		}
	}

	return x;
}
void PrintStatus(std::string location)
{
	std::vector<std::string> properties = {"Nt","Nm","healpix_order","needlet_order","Nl","Ns","hyperOrder","NVariancePopulations","totalRawParams","totalTransformedParams","mu_t","sigma_t","l_m","l_t","xt_Prior_normal","xt_Prior_gaps","xm_Prior"};
	std::vector<double> vals = {(double)Nt,(double)Nm,(double)healpix_order, (double)needlet_order, (double)Nl,(double)Ns,(double)hyperOrder,(double)NVariancePops,(double)totalRawParams,(double)totalTransformedParams,(double)xtPriorNonGap,(double)sigmat,(double)lm,(double)lt,(double)xtPriorNonGap,(double)xtPriorInsideGap,double(xmPrior)};
	
	std::fstream file;
	file.open(location + "/OptimiserProperties.dat",std::ios::out);
	int w = 15;
	
	for (int i = 0; i < properties.size(); ++i)
	{
		file << properties[i] << " = " << std::setprecision(10) << vals[i] << "\n";
	}
	file.close();
}

std::vector<bool> AssembleGapList()
{
	std::string gapFile = "../../ModelInputs/gaps_prior.dat";
	double timeFactor = (double)TotalScanningTime / Nt;
	int it = 0;
	bool inGap = false;
	int borderWidth = 0;
	int modifiedBorderWidth = borderWidth * timeFactor;
	bool inBorder= false;
	int trueTime = 0;
	int lastEnd = -9999;

	std::vector<bool> list(Nt,false);
	forLineVectorIn(gapFile,' ',
		
		int gapStart = std::stoi(FILE_LINE_VECTOR[0]);
		int gapEnd = std::stoi(FILE_LINE_VECTOR[1]);
		
		trueTime = floor(it * timeFactor);
		while (trueTime < gapEnd)
		{
			int leftDistance = std::min(abs(trueTime - gapStart),abs(trueTime - lastEnd));
			int rightDistance = abs(trueTime - gapEnd);
			
			bool inGap = (trueTime >= gapStart) && (trueTime <= gapEnd);
			
			bool nearGapEdge = (leftDistance < modifiedBorderWidth) || (rightDistance < modifiedBorderWidth);
			bool insertValue = false;;
			if (inGap)
			{
				insertValue = true;
				
			}
			//~ if (nearGapEdge)
			//~ {
				//~ insertValue = fa;
			//~ }
		
			list[it] = insertValue;
			
			++it;
			trueTime = floor((double)it * timeFactor);
			
		}
		lastEnd = gapEnd;
	);
	
	while (it<Nt)
	{
		list[it] = false;
		++it;
	}
		
	return list;
	
}
std::vector<bool> GapList = AssembleGapList();

#include "GlobalVariables.h"
#include "GlobalConstants.h"

void PrintStatus(std::string location)
{
	std::vector<std::string> properties = {"Nt","Nm","healpix_order","needlet_order","Nl","Ns","hyperOrder","NVariancePopulations","totalRawParams","totalTransformedParams","mu_t","sigma_t","l_m","l_t","xm_Prior"};
	std::vector<double> vals = {(double)Nt,(double)Nm,(double)healpix_order, (double)needlet_order, (double)Nl,(double)Ns,(double)hyperOrder,(double)NVariancePops,(double)totalRawParams,(double)totalTransformedParams,(double)xtPrior,(double)sigmat,(double)lm,(double)lt,double(xmPrior)};
	
	std::fstream file;
	file.open(location + "/OptimiserProperties.dat",std::ios::out);
	int w = 15;
	
	for (int i = 0; i < properties.size(); ++i)
	{
		file << properties[i] << " = " << std::setprecision(10) << vals[i] << "\n";
	}
	file.close();
}

std::vector<bool> AssembleGapList(bool buffered)
{
	std::string gapFile = "../../ModelInputs/gaps_prior.dat";
	double timeFactor = (double)TotalScanningTime/Nt;
	int it = 0;
	bool inGap = false;
	double borderWidthRevs = 0.03;
	int modifiedBorderWidth = borderWidthRevs*2160;
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

			if (buffered && nearGapEdge)
			{
				insertValue = false;
			}
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

//define this here so that the extern in the header file is assigned on startup
std::vector<bool> GapList = AssembleGapList(false);


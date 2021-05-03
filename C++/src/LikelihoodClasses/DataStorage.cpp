#include "DataStorage.h"

LikelihoodData::LikelihoodData(const std::vector<Star> &data, int id) : Stars(data)
{
	//initialise the frequently overwritten vectors
	
	int nBig = NumberLargerThanMaxObservations;
	
	pmf_forward = std::vector<std::vector<double>>(nBig,std::vector<double>(nBig,0));
	pmf_backward =  std::vector<std::vector<double>>(nBig,std::vector<double>(nBig,0));
	subpmf =  std::vector<std::vector<double>>(3,std::vector<double>(nBig,0));
	pt = std::vector<double>(nBig,0);
	pml = std::vector<double>(nBig,0);
	p = std::vector<double>(nBig,0);
	dfdp = std::vector<double>(nBig,0);
	
	NStars = Stars.size();
	
    //read in the healpix stuff
    std::string healpix_fov_file = "../../ModelInputs/scanninglaw_to_healpix_"+std::to_string(healpix_order)+".csv";
	healpix_fov_1 = std::vector<int>(TotalScanningTime,0);
	healpix_fov_2 = std::vector<int>(TotalScanningTime,0);
	int i = 0;
	forLineVectorInFile(healpix_fov_file,',',
		if (i > 0)
		{
	        healpix_fov_1[i] = std::stoi(FILE_LINE_VECTOR[1]);
	        healpix_fov_2[i] = std::stoi(FILE_LINE_VECTOR[2]);  
	        
	    }
	     ++i;
    );
    
    
    //prepare the temporal down stepping vector
    time_mapping = std::vector<int>(TotalScanningTime,0);
	
    double time_ratio = 1;
    if (Nt < TotalScanningTime)
    {
		time_ratio = (double)Nt/TotalScanningTime;
	}

	for (int i = 0; i < TotalScanningTime; ++i)
	{
		time_mapping[i] = floor(time_ratio*i);
	}
}

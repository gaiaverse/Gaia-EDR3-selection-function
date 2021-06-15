#include "DataStorage.h"

LikelihoodData::LikelihoodData(const std::vector<std::vector<Star>> &data, int id) : Stars(data)
{
	//initialise the frequently overwritten vectors
	ID = id;
	int nBig = NumberLargerThanMaxObservations;
	
	pmf_forward = std::vector<std::vector<double>>(nBig,std::vector<double>(nBig,0));
	pmf_backward =  std::vector<std::vector<double>>(nBig,std::vector<double>(nBig,0));
	subpmf =  std::vector<std::vector<double>>(3,std::vector<double>(nBig,0));
	pt = std::vector<double>(nBig,0);
	//~ ptm = std::vector<double>(nBig,0);
	pml = std::vector<double>(nBig,0);
	p = std::vector<double>(nBig,0);
	grad_elu_xml1 = std::vector<double>(nBig,0);
	grad_elu_xml2 = std::vector<double>(nBig,0);
	dfdp = std::vector<double>(nBig,0);
	
	//~ NStars = 0;
	//~ for (int i = 0; i < N_SGD_Batches; ++i)
	//~ {
		//~ NStars +=	Stars[i].size();
	//~ }
    //read in the healpix stuff
    std::string healpix_fov_file = "../../ModelInputs/scanninglaw_to_healpix_"+std::to_string(healpix_order)+".csv";
	healpix_fov_1 = std::vector<int>(TotalScanningTime,0);
	healpix_fov_2 = std::vector<int>(TotalScanningTime,0);
	int i = 0;
	forLineVectorIn(healpix_fov_file,',',
		if (i > 0)
		{
	        healpix_fov_1[i] = std::stoi(FILE_LINE_VECTOR[1]);
	        healpix_fov_2[i] = std::stoi(FILE_LINE_VECTOR[2]);  
	        
	    }
	     ++i;
    );
    
    
    //prepare the temporal down stepping vector
    time_mapping = std::vector<int>(TotalScanningTime,0);
	//~ magtime_mapping = std::vector<int>(TotalScanningTime,0);
    
    double time_ratio = 1;
    double magtime_ratio = 1;
    if (Nt < TotalScanningTime)
    {
		time_ratio = (double)Nt/TotalScanningTime;
		//~ magtime_ratio = (double)Nt_m/TotalScanningTime;
	}

	for (int i = 0; i < TotalScanningTime; ++i)
	{
		time_mapping[i] = std::min(Nt-1,(int)round(time_ratio*i));
		//~ magtime_mapping[i] = std::min(Nt_m-1,(int)round(magtime_ratio*i));
	}
	
	Mode = NormalApproximation;
	
	for(int i = 0; i < VariancePopulations.size(); ++i)
	{
		VariancePopulation p = VariancePopulation(VariancePopulationFractions[i],VarianceBaselines[i],VarianceLinears[i],VarianceQuadratics[i]);
		VariancePopulations.push_back(p);	
	}
}

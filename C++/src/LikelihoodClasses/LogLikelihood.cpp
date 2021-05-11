#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<std::vector<Star>> &data, int id): Data(data,id)
{
	Value = 0.0;
	StarsUsed = 0;
	Gradient = std::vector<double>(totalRawParams,0.0);
	TransformedPosition = std::vector<double>(totalTransformedParams,0.0);
	TransformedGradient = std::vector<double>(totalTransformedParams,0.0);
	
	MakeCovarianceMatrix();
}

void LogLikelihood::MakeCovarianceMatrix()
{
	std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
	int i = 0;
    forLineVectorInFile(needlet_file,',',
		if (i > 0)
		{
	        needlet_u.push_back(std::stoi(FILE_LINE_VECTOR[0]));
	        needlet_v.push_back(std::stoi(FILE_LINE_VECTOR[1]));
	        needlet_w.push_back(std::stod(FILE_LINE_VECTOR[2]));
		}
        ++i;
    );    
    needletN = needlet_u.size();
    bVector = std::vector<double>(Nm*Ns,0);
	
	
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
	CholeskyKg = Kg.llt().matrixL();

}


void LogLikelihood::ForwardTransform(const std::vector<double> &z)
{
	
	// Forward transformation
	double u = exp(-1.0/lt);
	double ua = 1.0/sqrt(1.0-u*u);
	double ub = -u*ua;
	double previous = z[Nt-1]; // First case is trivial
	TransformedPosition[Nt-1] = mut + sigmat * previous;
	for (int i = Nt - 2; i >= 0; i--) 
	{
    	previous = (z[i] - ub * previous) / ua;
    	TransformedPosition[i] = mut + sigmat * previous;
	}

	// bms = Lmnzns
	for (int s = 0; s < Ns; ++s)
	{
		for (int m = 0; m < Nm; ++m)
		{
			bVector[s*Nm+m] = 0;
			for (int n = 0; n <= m; ++n)
			{
				bVector[s*Nm+m] += CholeskyKg(m,n) * z[Nt+s*Nm+n];
			}
		}
	}

	// yml
	for (int i = 0; i < needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			TransformedPosition[Nt+needlet_u[i]*Nm+m] += needlet_w[i]*bVector[needlet_v[i]*Nm+m];
		}
	}
}


void LogLikelihood::BackwardTransform()
{
	// Backward transformation
	double u = exp(-1.0/lt);
    
	double ua = 1.0/sqrt(1.0-u*u);
	double ub = -u*ua;
	double sigmata = sigmat/ua;

	if (Nt > 1)
	{
	    Gradient[0] = sigmata * TransformedGradient[0];
	    for (int i = 1; i < Nt-1; i++) {
	        Gradient[i] = u * Gradient[i-1] + sigmata * TransformedGradient[i];
	    }
	    Gradient[Nt-1] = -ub * Gradient[Nt-2] + sigmat * TransformedGradient[Nt-1];
	}
	else
	{
		Gradient[0] = TransformedGradient[0]*sigmat;
	}
	// yml
	
	std::fill(bVector.begin(), bVector.end(),0);
	for (int i = 0; i < needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			bVector[needlet_v[i]*Nm+m] += needlet_w[i]*TransformedGradient[Nt+needlet_u[i]*Nm+m];
		}
	}

	// bms = Lmnzns
	for (int s = 0; s < Ns; ++s)
	{
		for (int m = 0; m < Nm; ++m)
		{
			for (int n = 0; n <=m; ++n)
			{
				Gradient[Nt+s*Nm+n] += CholeskyKg(m,n)*bVector[s*Nm+m];
			}
		}
	}

}


void LogLikelihood::Calculate(const std::vector<double> & x, int effectiveBatchID, int effectiveBatches)
{

	Reset();	

	ForwardTransform(x);

	int realBatchesPerEffective = N_SGD_Batches / effectiveBatches;
	
	int start = effectiveBatchID * realBatchesPerEffective;
	int end = N_SGD_Batches;
	if (effectiveBatchID < effectiveBatches-1)
	{
		end = (effectiveBatchID+1) * realBatchesPerEffective;
	}
	StarsUsed = 0;
	for (int i = start; i < end; ++i)
	{
		int n = Data.Stars[i].size();
		for (int j = 0; j < n; ++j)
		{
			PerStarContribution(i,j,x);
		}
		StarsUsed += n;
	}
	
	BackwardTransform();
}

void LogLikelihood::Reset()
{
	Value = 0;
	std::fill(Gradient.begin(),Gradient.end(),0.0);
	std::fill(TransformedGradient.begin(),TransformedGradient.end(),0.0);
	std::fill(TransformedPosition.begin(),TransformedPosition.end(),0.0);
}


void LogLikelihood::PerStarContribution(int batchId, int starID, const std::vector<double> & x)
{

	const Star * candidate = &Data.Stars[batchId][starID];

	GeneratePs(candidate,x);
	
	GenerateContribution(candidate);
		
	AssignGradients(candidate);
}

void LogLikelihood::GeneratePs(const Star * candidate, const std::vector<double> & x)
{
	int n = candidate->nVisit;
	//generate p vectors
	for (int i = 0; i < n; ++i)
	{
		int t= candidate->TimeSeries[i];

		double xt = x[Data.time_mapping[t]];
		int idx1 = Nt + Data.healpix_fov_1[t] * Nm + candidate->gBin;
		int idx2 = Nt + Data.healpix_fov_2[t] * Nm + candidate->gBin;
		double xml1 = x[idx1];
		double xml2 = x[idx2];
		
		Data.pt[i] = sigmoid(xt);
		Data.pml[i] = sigmoid(xml1 + xml2);
		Data.p[i] = Data.pt[i] *Data.pml[i];
	}
}

void LogLikelihood::GenerateContribution(const Star * candidate)
{
	int n = candidate->nVisit;
	int k = candidate->nMeasure;

	// lots of probability black magic stuff in this function
	// Ask Douglas for help!
	
	poisson_binomial_pmf_forward(Data.p,n,Data.pmf_forward);
	
	
	double correction = 1;
	for (int i = 0; i < PipelineMinVisits; ++i)
	{
		correction -= Data.pmf_forward[n-1][i];
	}
	if (correction < VerySmallNumber)
	{
		correction = 0;
		for (int i = n; i >= PipelineMinVisits; --i)
		{
			correction += Data.pmf_forward[n-1][i];
		}
	}
	
	double likelihood = Data.pmf_forward[n-1][k];
	double log_likelihood = log(likelihood);
	double log_correction = log(correction);
	
	double contribution = log_likelihood - log_correction;



	bool floatingPointTruncated = (Data.pmf_forward[n-1][0] < VerySmallNumber)  || ( Data.pmf_forward[n-1][n] < VerySmallNumber);
	bool somethingVeryBad = std::isnan(contribution) || std::isinf(contribution);

	if (floatingPointTruncated || somethingVeryBad)
	{
		GenerateExactContribution(candidate);
		return;
	}
	Value += contribution;


	//everything else is to calculate the gradients....
	poisson_binomial_pmf_backward(Data.p,n,Data.pmf_backward);
	poisson_binomial_subpmf(PipelineMinVisits-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[0]);
	
	bool measuredAtLeastOnce = false;
	bool missedAtLeastOnce = false;
	if (k > 0)
	{
		poisson_binomial_subpmf(k-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[1]);
		measuredAtLeastOnce = true;
	}
	if (k < n)
	{
		poisson_binomial_subpmf(k,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[2]);
		missedAtLeastOnce = true;
	}
	
	//plonk the gradients into the vector to be used in the assignment functions
	for (int i = 0; i < n; ++i)
	{
		double dfdp_i =  - Data.subpmf[0][i]/correction;
		
		if (measuredAtLeastOnce)
		{
			dfdp_i += Data.subpmf[1][i]/likelihood;
		}
		if (missedAtLeastOnce)
		{
			dfdp_i -= Data.subpmf[2][i]/likelihood;
		}
		
		if (std::isnan(dfdp_i) || std::isinf(dfdp_i))
		{
			GenerateExactContribution(candidate);
			return;
		}
		Data.dfdp[i] =  dfdp_i;
	}
	
	
}

void LogLikelihood::GenerateExactContribution(const Star * candidate)
{
	int n = candidate->nVisit;
	int k = candidate->nMeasure;

	// Even more black magic happening here - structurally the same as above, but with logs happening all over the shop
	// Ask Douglas for help!
	
	poisson_binomial_lpmf_forward(Data.p,n,Data.pmf_forward);
	poisson_binomial_lpmf_backward(Data.p,n,Data.pmf_backward);
	poisson_binomial_sublpmf(PipelineMinVisits-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[0]);
	
	double log_likelihood = Data.pmf_forward[n-1][k];
	double likelihood = exp(log_likelihood);
	
	double log_correction = 0;
	for (int i = 0; i < PipelineMinVisits; ++i)
	{
		log_correction += log1p(-exp(Data.pmf_forward[n-1][i]-log_correction));
	}
	bool triggeredLoopInvert = false;
	if (std::isinf(log_correction) || std::isnan(log_correction))
	{
		triggeredLoopInvert = true;
		log_correction = VerySmallLog;
		for (int i = n; i >= PipelineMinVisits; --i)
		{
			log_correction = log_add_exp(log_correction,Data.pmf_forward[n-1][i]);
		}
	}
	double correction = exp(log_correction);
	
	double contribution = log_likelihood - log_correction;
	

	Value += contribution;
	

	bool measuredAtLeastOnce = false;
	bool missedAtLeastOnce = false;
	if (k > 0)
	{
		poisson_binomial_sublpmf(k-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[1]);
		measuredAtLeastOnce = true;
	}
	if (k < n)
	{
		poisson_binomial_sublpmf(k,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[2]);
		missedAtLeastOnce = true;
	}
	
	bool dfdpEmergency = false;
	for (int i = 0; i < n; ++i)
	{
		double dfdp_i =   - exp(Data.subpmf[0][i] - log_correction);
		if (measuredAtLeastOnce)
		{
				dfdp_i += exp(Data.subpmf[1][i]-log_likelihood);
		}
		if (missedAtLeastOnce)
		{
			dfdp_i -= exp(Data.subpmf[2][i]-log_likelihood);
		}
		
		
		if (std::isnan(dfdp_i) || std::isinf(dfdp_i))
		{
			dfdpEmergency = true;
		}
		Data.dfdp[i] = dfdp_i;
	}
	
	if (std::isnan(contribution) || std::isinf(contribution) || dfdpEmergency)
	{
		std::cout << "\n\n Error detected! NaN found in Value calculation on core " << Data.ID << "\n";
		std::cout << "n = " << n << "k = " << k << std::endl;
		std::cout << "likelihood = " << likelihood << "  correction = " << correction << "\n";
		std::cout << "triggered loop inversion: " << triggeredLoopInvert << "\n";
		std::cout << "log_likelihood = " << log_likelihood << "  log_correction = " << log_correction << "\n";
		std::cout << "p = (";
		for (int i = 0; i < n; ++i)
		{
			std::cout << Data.p[i] << ", ";
		}
		std::cout << ")\n\nl_pmf = (";
		for (int i = 0; i <= n; ++i)
		{
			double v = Data.pmf_forward[n-1][i];

			std::cout << v << ",";
		}
		
		std::cout << "\nsub_pmf_0 \t\tsub_pmf_1\t\tsub_pmf_2";
		for (int i = 0; i < n; ++i)
		{
			std::cout << Data.subpmf[0][i] << "\t\t" << Data.subpmf[1][i] << "\t\t" << Data.subpmf[2][i] << "\n";
		}
		
		std::cout << "\n\ndfdp = (";
		for (int i = 0; i < n; ++i)
		{
			std::cout << Data.dfdp[i] << ", ";
		}
		std::cout << ")\n\n";
		ERROR(100, "See above output");
	}
}

void LogLikelihood::AssignGradients(const Star * candidate)
{
	int n = candidate->nVisit;
	for (int i = 0; i < n; ++i)
	{		
		double dFdP_p = Data.dfdp[i] * Data.p[i];

		int t = candidate->TimeSeries[i];

		int offset = Nt + candidate->gBin;
		double mlGrad = dFdP_p * (1.0 - Data.pml[i]);	
		int index1 = offset +  Data.healpix_fov_1[t] * Nm;
		int index2 = offset +  Data.healpix_fov_2[t] * Nm;
		
		TransformedGradient[Data.time_mapping[t]] += dFdP_p * (1.0 - Data.pt[i]);
		TransformedGradient[index1] += mlGrad;
		TransformedGradient[index2] += mlGrad;
	}
}

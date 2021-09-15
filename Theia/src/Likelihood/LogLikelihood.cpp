#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<std::vector<Star>> &data): Data(data)
{
	Value = 0.0;
	StarsUsed = 0;
	Gradient = std::vector<double>(totalTransformedParams,0.0);
}

void LogLikelihood::Calculate(const EfficiencyVector & x, int effectiveBatchID, int effectiveBatches, int maxBatches)
{

	
	Reset();	
	Data.GeneratePopulations(x);
	int realBatchesPerEffective = maxBatches / effectiveBatches;
	
	int start = effectiveBatchID * realBatchesPerEffective;
	int end = maxBatches;
	if (effectiveBatchID < effectiveBatches-1)
	{
		end = (effectiveBatchID+1) * realBatchesPerEffective;
	}
	StarsUsed = 0;
	
	std::string p = "I am preparing to calculate batches ";
	for (int i = start; i < end; ++i)
	{
		p += std::to_string(i) + ", ";
	}
	p += "\nMy vectors have the following sizes:";
	p+= "  Data[i] " + std::to_string(Data.Stars[start].size());
	p+= "  dfdp(N) " + std::to_string(Data.dfdp_constantN.size());
	p += "  hypergradient " + std::to_string(Data.hypergradient.size());
	p +=  "  pt " + std::to_string(Data.pt.size());
	p += "   pml  " + std::to_string(Data.pml.size());
	p += "   p " + std::to_string(Data.p.size() );
	p += "   elu_grad  " + std::to_string(Data.grad_elu_xml1.size() ) + "  " + std::to_string(Data.grad_elu_xml2.size() );
	p += "  Varpop  " + std::to_string(Data.VariancePopulations.size() );
	p += "   x.transform  " + std::to_string(x.TransformedGradient.size() ) + "   " + std::to_string(x.TransformedPosition.size() );
	p += "   x.raw  " + std::to_string(x.RawGradient.size() ) + "   " + std::to_string(x.RawPosition.size() );
	p += "   healpix  " + std::to_string(Data.healpix_fov_1.size() ) + "   " + std::to_string(Data.healpix_fov_2.size() );
	std::cout << p << std::endl;
	
	for (int i = start; i < end; ++i)
	{
		int n = Data.Stars[i].size();
		for (int j = 0; j < n; ++j)
		{
			PerStarContribution(i,j,x);
		}
		StarsUsed += n;
	}
}

void LogLikelihood::Reset()
{
	Value = 0;
	std::fill(Gradient.begin(),Gradient.end(),0.0);
}



void LogLikelihood::PerStarContribution(int batchId, int starID, const EfficiencyVector & x)
{
	const Star candidate = Data.Stars[batchId][starID];

	GeneratePs(candidate,x);
	
	GenerateContribution(candidate);
	AssignGradients(candidate);
}


void LogLikelihood::GeneratePs(const Star candidate, const EfficiencyVector & x)
{
	int n = candidate.nVisit;
	//generate p vectors

	for (int i = 0; i < n; ++i)
	{
		int t= candidate.TimeSeries[i];
		int T= Data.time_mapping[t];
		
		//double time_multiplier = time_ratio * t - T;
		//double xt = (1.0-time_multiplier) * x[T] + time_multiplier * x[T+1];
		double xt = x.Access(x.Transformed,x.Temporal,x.Position,T);

		int idx1 = Nt + Data.healpix_fov_1[t] * Nm + candidate.gBin;
		int idx2 = Nt + Data.healpix_fov_2[t] * Nm + candidate.gBin;
		
		double xlm1 =  x.Access(x.Transformed,x.Spatial,x.Position,Data.healpix_fov_1[t], candidate.gBin);
		double xlm2 =  x.Access(x.Transformed,x.Spatial,x.Position,Data.healpix_fov_2[t], candidate.gBin);
		
		double elu_xml1 = elu(xlm1);
		double elu_xml2 = elu(xlm2);
		
		double pt = sigmoid(xt);
		Data.pt[i] = pt;

        Data.grad_elu_xml1[i] = elu_grad(xlm1, elu_xml1);
        Data.grad_elu_xml2[i] = elu_grad(xlm2, elu_xml2);
		Data.pml[i] = exp(-spatialAddingPrefactor * (elu_xml1 + elu_xml2));
		Data.p[i] = Data.pt[i] *Data.pml[i];
	}
}

void LogLikelihood::GenerateContribution(const Star candidate)
{
	// lots of probability black magic stuff in this function
	// Ask Douglas for help!
	switch(Data.Mode)
	{
		case PoissonBinomial:
		{
			PoissonContribution(candidate);
			break;
		}
		case NormalApproximation:
		{
			NormalContribution(candidate);
			break;
		}
	}
}
void LogLikelihood::NormalContribution(const Star candidate)
{
	int n = candidate.nVisit;
	int k = candidate.nMeasure;
	Value += poisson_binomial_normal_lpmf(k, n, Data);
}

void LogLikelihood::PoissonContribution(const Star candidate)
{
	int n = candidate.nVisit;
	int k = candidate.nMeasure;

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
		ExactPoissonContribution(candidate);
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
			ExactPoissonContribution(candidate);
			return;
		}
		Data.dfdp_constantN[i] =  dfdp_i;
	}
	
}

void LogLikelihood::ExactPoissonContribution(const Star candidate)
{
	int n = candidate.nVisit;
	int k = candidate.nMeasure;

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
		Data.dfdp_constantN[i] = dfdp_i;
	}
	
}
void LogLikelihood::AssignGradients(const Star candidate)
{
	int n = candidate.nVisit;
	for (int i = 0; i < n; ++i)
	{		
		int t= candidate.TimeSeries[i];
		int T= Data.time_mapping[t];
		
		int offset = Nt + candidate.gBin;
		int index1 = offset +  Data.healpix_fov_1[t] * Nm;
		int index2 = offset +  Data.healpix_fov_2[t] * Nm;
		
		double dndp_time, dndp_space;
		switch (ScalingMode)
		{
			case (NScaling):
			{
				dndp_time = 0;
				dndp_space = 0;
				break;
			}
			case (MScaling):
			{
				dndp_time = Data.pml[i];
				dndp_space = Data.pt[i];
				break;
			}
			
			case (ActiveNScaling):
			{
				dndp_time = 1.0;
				dndp_space = 0.0;
			}
		}
		double dfdP_time = Data.pml[i] * Data.dfdp_constantN[i] + Data.dfdN_constantP * dndp_time;
		double dfdP_space = Data.pt[i] * Data.dfdp_constantN[i] + Data.dfdN_constantP * dndp_space;
		
		
		
		Gradient[T] += Data.pt[i] * (1.0 - Data.pt[i]) * dfdP_time;
		Gradient[index1] -= spatialAddingPrefactor * Data.grad_elu_xml1[i] * Data.pml[i] * dfdP_space;
		Gradient[index2] -= spatialAddingPrefactor * Data.grad_elu_xml2[i] * Data.pml[i] * dfdP_space;
	}
	
	for (int i = 0; i < NHyper; ++i)
	{
		Gradient[transformedNonHyperParams + i] += Data.hypergradient[i];
	}
	
}

#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<std::vector<Star>> &data, int id): Data(data,id)
{
	Value = 0.0;
	StarsUsed = 0;
	Gradient = std::vector<double>(totalTransformedParams,0.0);
}

void LogLikelihood::Calculate(const std::vector<double> & x, int effectiveBatchID, int effectiveBatches, int maxBatches)
{

	Reset();	


	int realBatchesPerEffective = maxBatches / effectiveBatches;

	
	
	int start = effectiveBatchID * realBatchesPerEffective;
	int end = maxBatches;
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
}

void LogLikelihood::Reset()
{
	Value = 0;
	std::fill(Gradient.begin(),Gradient.end(),0.0);
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
		int T= Data.time_mapping[t];
		
		//double time_multiplier = time_ratio * t - T;
		//double xt = (1.0-time_multiplier) * x[T] + time_multiplier * x[T+1];
		double xt = x[T];
		double xtm = 100000.0;
		if (Nt_m > 0)
		{
			int Tm = Data.magtime_mapping[t];
			xtm = x[Nt + Nl*Nm + Tm * Nm + candidate->gBin];
		}
		int idx1 = Nt + Data.healpix_fov_1[t] * Nm + candidate->gBin;
		int idx2 = Nt + Data.healpix_fov_2[t] * Nm + candidate->gBin;
		double elu_xml1 = elu(x[idx1]);
		double elu_xml2 = elu(x[idx2]);
		
		Data.pt[i] = sigmoid(xt);
		Data.ptm[i] = sigmoid(xtm);
        Data.grad_elu_xml1[i] = elu_grad(x[idx1], elu_xml1);
        Data.grad_elu_xml2[i] = elu_grad(x[idx2], elu_xml2);
		Data.pml[i] = exp(-density_alpha * (elu_xml1 + elu_xml2));
		Data.p[i] = Data.pt[i] *Data.pml[i] * Data.ptm[i];
	}
}

void LogLikelihood::GenerateContribution(const Star * candidate)
{
	int n = candidate->nVisit;
	int k = candidate->nMeasure;

	// lots of probability black magic stuff in this function
	// Ask Douglas for help!
	double contribution;
	poisson_binomial_normal_lpmf(k, Data.p, n, contribution, Data.dfdp);
	Value += contribution;
	
}

void LogLikelihood::AssignGradients(const Star * candidate)
{
	int n = candidate->nVisit;
	for (int i = 0; i < n; ++i)
	{		
		double dFdP_p = Data.dfdp[i] * Data.p[i];

		int t= candidate->TimeSeries[i];
		int T= Data.time_mapping[t];
		
		//double time_multiplier = time_ratio * t - T;

		int offset = Nt + candidate->gBin;
		int index1 = offset +  Data.healpix_fov_1[t] * Nm;
		int index2 = offset +  Data.healpix_fov_2[t] * Nm;
		
		//Gradient[T] += dFdP_p * (1.0 - Data.pt[i]) * (1.0 - time_multiplier);
		//Gradient[T+1] += dFdP_p * (1.0 - Data.pt[i]) * time_multiplier;
		Gradient[T] += dFdP_p * (1.0 - Data.pt[i]);

		Gradient[index1] -= density_alpha * Data.grad_elu_xml1[i] * dFdP_p;
		Gradient[index2] -= density_alpha * Data.grad_elu_xml2[i] * dFdP_p;
		
		if (Nt_m > 0)
		{
			int Tm = Data.magtime_mapping[t];
			Gradient[offset+Nm*(Nl + Tm)] += dFdP_p * (1.0 - Data.ptm[i]); 
		}
	}
}

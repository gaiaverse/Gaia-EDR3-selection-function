#include "LogLikelihood.h"

LogLikelihood::LogLikelihood(const std::vector<Star> &data, int id): Data(data,id)
{

	Value = 0.0;
	Gradient = Eigen::VectorXd::Zero(totalTransformedParams);
}

void LogLikelihood::Calculate(Eigen::VectorXd& x)
{

	Reset();	

	for (int i = 0; i < Data.NStars; ++i)
	{
		PerStarContribution(i,x);
	}
}

void LogLikelihood::Reset()
{
	Value = 0;
	for (int i = 0; i < Gradient.size(); ++i)
	{
		Gradient[i] = 0;
	}
}


void LogLikelihood::PerStarContribution(int star, Eigen::VectorXd& x)
{

	const Star * candidate = &Data.Stars[star];

	GeneratePs(candidate,x);
	
	GenerateContribution(candidate);
		
	AssignGradients(candidate);
}

void LogLikelihood::GeneratePs(const Star * candidate, Eigen::VectorXd & x)
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
	
	double zeroMeasureKiller = 0;
	double nMeasureKiller = 0;
	if (k > 0)
	{
		poisson_binomial_subpmf(k-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[1]);
		zeroMeasureKiller = 1;
	}
	if (k < n)
	{
		poisson_binomial_subpmf(k,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[2]);
		nMeasureKiller = 1;
	}
	
	//plonk the gradients into the vector to be used in the assignment functions
	for (int i = 0; i < n; ++i)
	{
		Data.dfdp[i] =  (Data.subpmf[1][i]*zeroMeasureKiller-Data.subpmf[2][i]*nMeasureKiller)/likelihood - Data.subpmf[0][i]/correction;
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
			log_correction += log1p(exp(Data.pmf_forward[n-1][i] - log_correction));
		}
	}
	double correction = exp(log_correction);
	
	double contribution = log_likelihood - log_correction;
	
	
	if (std::isnan(contribution) || std::isinf(contribution))
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
		for (int i = 0; i < n; ++i)
		{
			double v = Data.pmf_forward[n-1][i];

			std::cout << v << ",";
		}
		std::cout << ") \n\n\n I will now repeatt he calculation which triggered this error....\n\nForward loop:";
		log_correction = 0;
		for (int i = 0; i < PipelineMinVisits; ++i)
		{
			double c = log1p(-exp(Data.pmf_forward[n-1][i]-log_correction));
			std::cout << i << "   pmf = " << Data.pmf_forward[n-1][i] << "  dc = " << c;
			log_correction += c;
			std::cout << "   l_c=" << log_correction << "\n";
		}
		
		if (std::isnan(log_correction) || std::isinf(log_correction) )
		{
				std::cout << "\nBackwards loop: \n";
				log_correction = VerySmallLog;
				for (int i = n; i >= PipelineMinVisits; --i)
				{
					double c = log1p(exp(Data.pmf_forward[n-1][i] - log_correction));
					std::cout << i << "   pmf = " << Data.pmf_forward[n-1][i] << "  dc = " << c;
					log_correction += c;
					std::cout << "   l_c=" << log_correction << "\n";
				}
			
		}
		ERROR(100, "See above output");
	}
	
	

	Value += contribution;
	

	double zeroMeasureKiller = 0;
	double nMeasureKiller = 0;
	if (k > 0)
	{
		poisson_binomial_sublpmf(k-1,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[1]);
		zeroMeasureKiller = 1;
	}
	if (k < n)
	{
		poisson_binomial_sublpmf(k,n,Data.pmf_forward,Data.pmf_backward,Data.subpmf[2]);
		nMeasureKiller = 1;
	}
	
	for (int i = 0; i < n; ++i)
	{
		Data.dfdp[i] =  exp(Data.subpmf[1][i]-log_likelihood)*zeroMeasureKiller - exp(Data.subpmf[2][i]-log_likelihood)*nMeasureKiller - exp(Data.subpmf[0][i] - log_correction);
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
		
		Gradient[Data.time_mapping[t]] += dFdP_p * (1.0 - Data.pt[i]);
		Gradient[index1] += mlGrad;
		Gradient[index2] += mlGrad;
	}
}

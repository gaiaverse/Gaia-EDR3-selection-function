#include "EfficiencyVector.h"


EfficiencyVector::EfficiencyVector(std::string load_location)
{
	bool loadIn = !(load_location == "__null_location__");
	RawPosition = std::vector<double>(totalRawParams,0.0);
	RawGradient= std::vector<double>(totalRawParams,0.0);
	TransformedPosition= std::vector<double>(totalTransformedParams,0.0);
	TransformedGradient= std::vector<double>(totalTransformedParams,0.0);
	if (loadIn)
	{
		LoadVector(load_location);
	}
	else
	{
		GenerateVector();	
	}
	LoadNeedlets();
	LoadCholesky();
}


double EfficiencyVector::Access(VectorMode mode, VectorComponent component, VectorType type, int i)
{	
	int index;
	switch (component)
	{
		case Temporal:
			index = i;
			break;
		case Spatial:
			index = Nt + i;
			break;
		case Hyper:
			index = rawNonHyperParams + i;
			if (mode == Transformed)
			{
				index = transformedNonHyperParams + i;
			}
			break;
	}

	double v;
	
	if (mode == Raw)
	{
		if (type == Position)
		{
			v = RawPosition[index];
		}
		if (type == Gradient)
		{
			v = RawGradient[index];
		}
	}
	if (mode == Transformed)
	{
		if (type == Position)
		{
			v = TransformedPosition[index];
		}
		if (type == Gradient)
		{
			v = TransformedGradient[index];
		}
	}
	return v;
}

double EfficiencyVector::Access(VectorMode mode, VectorComponent component, VectorType type, int sl, int m)
{
	switch (component)
	{
		case Temporal:
			ERROR(2, "Incorrect access mode used for temporal component of the EfficiencyVector");
			break;
		case Spatial:
			Access(mode, component,type, sl*Nm + m);
			break;
		case Hyper:
			ERROR(2, "Incorrect access mode used for hyper component of the EfficiencyVector");
			break;
	}
}

void EfficiencyVector::Assign(VectorMode mode, VectorComponent component, VectorType type, int i, double newValue)
{	
	int index;
	switch (component)
	{
		case Temporal:
			index = i;
			break;
		case Spatial:
			index = Nt + i;
			break;
		case Hyper:
			index = rawNonHyperParams + i;
			if (mode == Transformed)
			{
				index = transformedNonHyperParams + i;
			}
			break;
	}
	
	if (mode == Raw)
	{
		if (type == Position)
		{
			RawPosition[index] = newValue;
		}
		if (type == Gradient)
		{
			RawGradient[index] = newValue;
		}
	}
	if (mode == Transformed)
	{
		if (type == Position)
		{
			TransformedPosition[index] = newValue;
		}
		if (type == Gradient)
		{
			TransformedGradient[index] = newValue;
		}
	}
	
}

void EfficiencyVector::Assign(VectorMode mode, VectorComponent component, VectorType type,int sl, int m, double newValue)
{
	switch (component)
	{
		case Temporal:
			ERROR(2, "Incorrect access mode used for temporal component of the EfficiencyVector");
			break;
		case Spatial:
			Assign(mode, component,type,  sl*Nm + m, newValue);
			break;
		case Hyper:
			ERROR(2, "Incorrect access mode used for hyper component of the EfficiencyVector");
			break;
	}
}

void EfficiencyVector::Increment(VectorMode mode, VectorComponent component, VectorType type,int i, double value)
{	
	int index;
	switch (component)
	{
		case Temporal:
			index = i;
			break;
		case Spatial:
			index = Nt + i;
			break;
		case Hyper:
			index = rawNonHyperParams + i;
			if (mode == Transformed)
			{
				index = transformedNonHyperParams + i;
			}
			break;
	}

	
	if (mode == Raw)
	{
		if (type == Position)
		{
			RawPosition[index] += value;
		}
		if (type == Gradient)
		{
			RawGradient[index] += value;
		}
	}
	if (mode == Transformed)
	{
		if (type == Position)
		{
			TransformedPosition[index] += value;
		}
		if (type == Gradient)
		{
			TransformedGradient[index] += value;
		}
	}
}

void EfficiencyVector::Increment(VectorMode mode, VectorComponent component, VectorType type,int sl, int m, double value)
{
	switch (component)
	{
		case Temporal:
			ERROR(2, "Incorrect access mode used for temporal component of the EfficiencyVector");
			break;
		case Spatial:
			Increment(mode, component,type,  sl*Nm + m, value);
			break;
		case Hyper:
			ERROR(2, "Incorrect access mode used for hyper component of the EfficiencyVector");
			break;
	}
}

void EfficiencyVector::LoadVector(std::string load_location)
{
	int i = 0;
	forLineIn(load_location,
		if (i >= totalRawParams)
		{
			ERROR(100,"Asked to load in start position from file, but it was the wrong length");
		}
		
		double value; 
		try
		{
			value = std::stod(FILE_LINE);
		}
		catch(const std::exception& e)
		{
			value = 0;
			//sometimes the stod fails for values <1e-300, so just default these to zero
		}
		
		RawPosition[i] = value;
		++i;
	);
	
	bool wentUnder = (i < totalRawParams);
	if (wentUnder)
	{
		ERROR(100,"Asked to load in start position from file, but it was the wrong length");
	}
	
}

void EfficiencyVector::GenerateVector()
{
	//firstly seed a bunch of random numbers
	for (int i = 0; i < totalRawParams; ++i)
	{
		RawPosition[i] = 2*initialisationBounds*((double)std::rand()/RAND_MAX -0.5);
	}
	
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
		Increment(Raw,Spatial,Position,0,i,mums[i]);
	}
}

void EfficiencyVector::LoadNeedlets()
{
	std::string needlet_file = "../../ModelInputs/needlets_"+std::to_string(healpix_order)+"_"+std::to_string(needlet_order)+".csv";
	int i = 0;
    forLineVectorIn(needlet_file,',',
 
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
}
	
void EfficiencyVector::LoadCholesky()
{
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

	std::vector<double> max_in_row = std::vector<double>(Nm,0);
	for (int i = 0; i < Nm; i++) 
	{
		for (int j = 0; j <= i; j++) 
		{
			max_in_row[i] += std::max(0.0,abs(CholeskyKg(i,j))-max_in_row[i]);
		}
	}

	choleskyN = 0;
	for (int i = 0; i < Nm; i++) 
	{
		for (int j = 0; j <= i; j++) 
		{
			if (abs(CholeskyKg(i,j)) > cholesky_tol * max_in_row[i])
			{
				choleskyN += 1;
				cholesky_u.push_back(i);
				cholesky_v.push_back(j);
				cholesky_w.push_back(CholeskyKg(i,j));
			}
		}
	}
	
}

void EfficiencyVector::Reset()
{
	std::fill(TransformedPosition.begin(), TransformedPosition.end(),xmPrior);
	std::fill(TransformedGradient.begin(), TransformedGradient.end(),0);
	std::fill(RawGradient.begin(), RawGradient.end(),0);
}

void EfficiencyVector::ForwardTransform(const VectorXd &z)
{
	Reset();

	RawPosition = std::vector<double>(&z[0],z.data() + z.cols()*z.rows());

	ForwardTransform_Temporal();
	
	ForwardTransform_Spatial();
	
	ForwardTransform_Hyper();
}

void EfficiencyVector::ForwardTransform_Spatial()
{
	std::fill(bVector.begin(), bVector.end(),0);
	for (int s = 0; s < Ns; ++s)
	{
		for (int i = 0; i < choleskyN; ++i)
		{
			bVector[s*Nm+cholesky_u[i]] += cholesky_w[i] * Access(Raw,Spatial,Position,s,cholesky_v[i]);
		}
	}

	// yml
	for (int i = 0; i < needletN; ++i)
	{
		for (int m = 0; m < Nm; ++m)
		{
			Increment(Transformed,Spatial,Position, needlet_u[i],m, needlet_w[i]*bVector[needlet_v[i]*Nm+m]);
		}
	}
}

void EfficiencyVector::ForwardTransform_Temporal()
{

	double u = exp(-1.0/lt);
	double ua = sqrt(1.0-u*u);
	double previous = Access(Raw,Temporal,Position,Nt-1); // First case is trivial
	Assign(Transformed,Temporal,Position,Nt-1, xtPrior + sigmat * previous);
	for (int i = Nt - 2; i >= 0; i--) 
	{
    	previous = ua * Access(Raw,Temporal,Position,i) + u * previous;
    	Assign(Transformed, Temporal, Position,i, xtPrior + sigmat * previous);
	}
}

void EfficiencyVector::ForwardTransform_Hyper()
{
	// [[ Nt + Nl*Nm + (zeroth order weightings) + (first order weightings) + ... +(pop fractions) + popSum ]
	int offset = hyperFractionOffset;
	double X = VerySmallLog;
	
	for (int i = 0; i < NVariancePops; ++i)
	{
		X = log_add_exp(X,Access(Raw,Hyper,Position,offset + i));
		
	}
		
	for (int i = 0; i < NHyper; ++i)
	{
		double v = Access(Raw,Hyper,Position,i);
		if (i >= offset)
		{
			v = exp(v - X);
		}
		Assign(Transformed,Hyper,Position,i,v);
	}	
}

void EfficiencyVector::BackwardTransform()
{
	BackwardTransform_Temporal();
	
	//~ BackwardTransform_Spatial();
	
	//~ BackwardTransform_Hyper();
}

void EfficiencyVector::BackwardTransform_Temporal()
{
	double u = exp(-1.0/lt);
	double ua = 1.0/sqrt(1.0-u*u);
	double ub = -u*ua;
	double sigmata = sigmat/ua;
	
	if (Nt > 1)
	{
	    double previous = sigmata * Access(Transformed,Temporal,Gradient,0);
	    Assign(Raw,Temporal,Gradient,0,previous);
	    for (int i = 1; i < Nt-1; i++) 
	    {
			previous = u * previous + sigmata * Access(Transformed,Temporal,Gradient,i);
	        Assign(Raw,Temporal,Gradient,i,previous);
	    }
	    previous = -ub * previous + sigmat * Access(Transformed,Temporal,Gradient,Nt-1);
	    Assign(Raw,Temporal,Gradient,Nt-1, previous);
	}
	else
	{
		Assign(Raw,Temporal,Gradient,0,sigmat*Access(Transformed,Temporal,Gradient,0));
	}
}

//~ void DescentFunctor::BackwardTransform_Spatial()
//~ {
	//~ if (SpaceActive)
	//~ {
		//~ int insertOffset = 0;
		//~ if (TimeActive)
		//~ {
			//~ insertOffset += Nt;
		//~ }
		//~ // yml
		//~ std::fill(bVector.begin(), bVector.end(),0);
		//~ for (int i = 0; i < needletN; ++i)
		//~ {
			//~ for (int m = 0; m < Nm; ++m)
			//~ {
				//~ bVector[needlet_v[i]*Nm+m] += needlet_w[i]*TransformedGradient[Nt+needlet_u[i]*Nm+m];
			//~ }
		//~ }
		//~ // bms = Lmnzns
		//~ for (int s = 0; s < Ns; ++s)
		//~ {
			//~ for (int i = 0; i < L.choleskyN; ++i)
			//~ {
				//~ Gradient[insertOffset+s*Nm+L.cholesky_v[i]] += L.cholesky_w[i] * bVector[s*Nm+L.cholesky_u[i]];
			//~ }
		//~ }
	//~ }
//~ }

//~ void DescentFunctor::BackwardTransform_Hyper()
//~ {
	//~ if (HyperActive)
	//~ {
		//~ int insertOffset = 0;
		//~ if (TimeActive)
		//~ {
			//~ insertOffset+=Nt;
		//~ }
		//~ if (SpaceActive)
		//~ {
			//~ insertOffset += Ns*Nm;
		//~ }
		
		//~ for (int i = 0; i < hyperFractionOffset; ++i)
		//~ {
			//~ double x = TransformedPosition[transformedNonHyperParams + i];
			//~ double df = TransformedGradient[transformedNonHyperParams + i];

			
			//~ Gradient[insertOffset + i] = df;		
		//~ }
		
		//~ for (int i = 0; i < NVariancePops; ++i)
		//~ {
			//~ double xi = TransformedPosition[transformedNonHyperParams + hyperFractionOffset+ i];
			//~ double sum = 0;
			
			//~ for (int j = 0; j <  NVariancePops; ++j)
			//~ {
				//~ double xj = TransformedPosition[transformedNonHyperParams + hyperFractionOffset+ j];
				//~ double dfj = TransformedGradient[transformedNonHyperParams + hyperFractionOffset+ j];
				
				//~ double ijTerm = 0;
				//~ if (i==j)
				//~ {
					//~ ijTerm = 1;
				//~ }
				//~ sum += dfj * xi*(ijTerm -xj);
			//~ }
			//~ Gradient[insertOffset + hyperFractionOffset + i] = sum;
		//~ }	
	//~ }
//~ }	

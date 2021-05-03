#include "ManualOptimizer.h"


TestFunctor::TestFunctor(int n)
{
	Dimensions = n;
	Value = 0;
	Gradient = VectorXd::Zero(n);
}

void TestFunctor::Calculate(const VectorXd & x)
{
	Value = 0;
	
	for (int i = 0; i < Dimensions; ++i)
	{
		
		double d = x[i] - i;
		int n = 4;
		double v = pow(d,n);

		Gradient[i] = 10*pow(d,n-1);

		Value += v;
	}
}
void TestFunctor::operator () (const VectorXd & x)
{
	Calculate(x);
}


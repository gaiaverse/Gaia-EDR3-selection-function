#include "ManualOptimizer.h"

TestFunctor::TestFunctor()
{
}
TestFunctor::TestFunctor(int n)
{
	Dimensions = n;
	Value = 0;
	LoopID = 0;
	Gradient = VectorXd::Zero(n);
}

void TestFunctor::Calculate(const VectorXd & x)
{
	Value = 0;
	
	double a = 3;
	double b = 100;
	
	double u = x[0];
	double v = x[1];
	
	double au = a-u;
	double vu = v - u*u;
	Value = au * au + b*vu*vu;
	
	Gradient[0] = -2*au - 4 * b * vu*u;
	Gradient[1] = 2*b*vu;
	
	++LoopID;
}
void TestFunctor::operator () (const VectorXd & x)
{
	Calculate(x);
}

void TestFunctor::SavePosition(bool finalSave)
{
	
}

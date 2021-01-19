#include "Liklihood.h"

Liklihood::Liklihood(const std::vector<Star> &data, int nPoints, int id): Data(data)
{
	ID = id;
	Value = 0.0;
	Gradient = std::vector<double>(nPoints,0.0);

}


void Liklihood::Calculate(Eigen::VectorXd& x)
{
	double a = x[0];
	double b = x[1];
	double c = x[2];
	
	Value = 0;
	Gradient[0] = 0;
	Gradient[1] = 0;
	Gradient[2] = 0;
	
	

	for (int i = 0; i < Data.size(); ++i)
	{
		double arg = Data[i].z - Data[i].x*a - Data[i].y*b - c;
		
		double errSq = Data[i].err * Data[i].err;
		
		Value += 0.5 * arg*arg / errSq;
		Gradient[0] += -Data[i].x * arg / errSq;
		Gradient[1] += -Data[i].y * arg / errSq;
		Gradient[2] += -1 * arg / errSq;
	}

	
}

#include "MiscFunctions.h"


double  sigmoid(double x)
{
   // return 0.5*(1.0+tanh(0.5*x));
    if (x < 0.0)
    {
        double a = exp(x);
        return a / (1.0 + a); 
    }
    else
    {
        return 1.0 / (1.0 + exp(-x));
    }
}

/*
double  elu(double x)
{
	return exp(-x);
}

double  elu_grad(double x, double elu_x)
{
	return -elu_x;
}
*/

double  elu(double x)
{
    if (x < density_cut)
    {
        return (1.0 + density_cut - x)*expm_density_cut; 
    }
    else
    {
        return exp(-x);
    }
}

double  elu_grad(double x, double elu_x)
{
    if (x < density_cut)
    {
        return -expm_density_cut; 
    }
    else
    {
        return -elu_x;
    }
}


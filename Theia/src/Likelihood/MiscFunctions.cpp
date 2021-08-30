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

double  elu(double x)
{
    if (x < elu_transitionPoint)
    {
        return (1.0 + elu_transitionPoint - x)*exp_elu_transitionPoint; 
    }
    else
    {
        return exp(-x);
    }
}

double  elu_grad(double x, double elu_x)
{
    if (x < elu_transitionPoint)
    {
        return -exp_elu_transitionPoint; 
    }
    else
    {
        return -elu_x;
    }
}

double  log_add_exp(double a, double b)
{
    if (a > b)
    {
        return a + log1p(exp(b-a));
    }
    else
    {
        return b + log1p(exp(a-b));
    }
}

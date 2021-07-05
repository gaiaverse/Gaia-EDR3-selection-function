#pragma once
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
// Implements an expit sigmoid via the tanh method
double  sigmoid(double x);
double  elu(double x);
double  elu_grad(double x, double elu_x);


const double one_over_root2 = 1.0/sqrt(2.0);
const double one_over_root2pi = 1.0/sqrt(2.0*M_PI);
const double root2_over_pi = sqrt(2.0/M_PI);


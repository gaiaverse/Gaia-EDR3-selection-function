#pragma once
#include "../GenericFunctions/FileHandler.h"

#include "../GlobalVariables.h"



//function stolen from a git repo somewhere, see the implementation for more details
void  poisson_binomial_pmf_forward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_pmf_backward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_subpmf(int m, int probslen, std::vector<std::vector<double>> & pmf_forward, std::vector<std::vector<double>> & pmf_backward, std::vector<double> & result);

// Log-versions
double  log_add_exp(double a, double b);
void  poisson_binomial_lpmf_forward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_lpmf_backward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_sublpmf(int m, int probslen, std::vector<std::vector<double>> & lpmf_forward, std::vector<std::vector<double>> & lpmf_backward, std::vector<double> & result);

// Implements an expit sigmoid via the tanh method
double  sigmoid(double x);
double  elu(double x);
double  elu_grad(double x, double elu_x)

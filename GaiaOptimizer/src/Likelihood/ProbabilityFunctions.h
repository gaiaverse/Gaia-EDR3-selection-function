#pragma once
#include "../GenericFunctions/FileHandler.h"

#include "../Main/GlobalVariables.h"
#include <cmath>



//function stolen from a git repo somewhere, see the implementation for more details
void  poisson_binomial_pmf_forward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_pmf_backward(std::vector<double> &  probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_subpmf(int m, int probslen, std::vector<std::vector<double>> & pmf_forward, std::vector<std::vector<double>> & pmf_backward, std::vector<double> & result);

// Log-versions
double  log_add_exp(double a, double b);
void  poisson_binomial_lpmf_forward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_lpmf_backward(std::vector<double> & probs, int probslen, std::vector<std::vector<double>> & result);
void  poisson_binomial_sublpmf(int m, int probslen, std::vector<std::vector<double>> & lpmf_forward, std::vector<std::vector<double>> & lpmf_backward, std::vector<double> & result);

// Normal approximation to the Poisson-Binomial
void logphi(double z, double& f, double& df);
double poisson_binomial_normal_lpmf(int k, std::vector<double> & probs, int probslen, double& value, std::vector<std::vector<double>> & gradient);

// Implements an expit sigmoid via the tanh method
double  sigmoid(double x);
double  elu(double x);
double  elu_grad(double x, double elu_x);

std::vector<double> logphi_c{ 0.00048204, -0.00142906, 0.0013200243174, 0.0009461589032,
       -0.0045563339802, 0.00556964649138, 0.00125993961762116,
       -0.01621575378835404, 0.02629651521057465, -0.001829764677455021,
       2.0*(1.0-M_PI/3.0), (4.0-M_PI)/3.0, 1.0, 1.0 };
std::vector<double> logphi_r{1.2753666447299659525, 5.019049726784267463450,
        6.1602098531096305441, 7.409740605964741794425,
        2.9788656263939928886};
std::vector<double> logphi_q{2.260528520767326969592,  9.3960340162350541504,
       12.048951927855129036034, 17.081440747466004316,
        9.608965327192787870698,  3.3690752069827527677};

double one_over_root2 = 1.0/sqrt(2.0);
double one_over_root2pi = 1.0/sqrt(2.0*M_PI);
double root2_over_pi = np.sqrt(2.0/M_pi);

#pragma once
#include <vector>
#include "libs/Eigen/Core"
using Eigen::VectorXd;

const int Nh = 2; // number of hyper-hyper-parameters
const int Ng = 214;//35; // number of magnitude bins
const int Nt = 8967691; // number of time bins
const double SingularityPreventer = 10e-13;

const double mu_mean = -3.0;
const double mu_variance = 1.0;
const double lg = 0.1;

const int PipelineMinVisits = 5; 

#define FILEGAP << ", " << 

Eigen::VectorXd initialisedVector(int n);

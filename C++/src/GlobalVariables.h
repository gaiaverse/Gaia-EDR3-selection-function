#pragma once
#include <vector>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "libs/Eigen/Core"
using Eigen::VectorXd;

const int Nh = 2; // number of hyper-hyper-parameters
const int Ng = 81;//35; // number of magnitude bins
const int Nt = 896769	; // number of time bins
const double SingularityPreventer = 10e-13;

const double mu_mean = -3.0;
const double mu_variance = 1.0;
const double lg = 0.1;

const int PipelineMinVisits = 5; 

#define FILEGAP << ", " << 

Eigen::VectorXd initialisedVector(int n);

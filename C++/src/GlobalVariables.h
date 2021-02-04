#pragma once
#include <vector>
#include "libs/Eigen/Core"
using Eigen::VectorXd;

const int Nh = 2; // number of hyper-hyper-parameters
const int Ng = 81;//35; // number of magnitude bins
const int Nt = 1000;//8967691; // number of time bins
const double SingularityPreventer = 10e-13;
//~ const std::vector<double> magnitudes = { 2.5, 10.5, 16.5, 17.5, 18.1, 18.3, 18.5, 18.7, 18.9, 19.05, 19.15, 19.25, 19.35, 19.45, 19.55, 19.65, 19.75, 19.85, 19.95, 20.05, 20.15, 20.25, 20.35, 20.45, 20.55, 20.65, 20.75, 20.85, 20.95, 21.05, 21.15, 21.25, 21.35, 21.45, 23.25}; // magniitude bins

#define FILEGAP << ", "<< 

Eigen::VectorXd initialisedVector(int n);

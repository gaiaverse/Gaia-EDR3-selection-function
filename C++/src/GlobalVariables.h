#pragma once
#include <vector>
#include "libs/Eigen/Core"
using Eigen::VectorXd;

const int Nm = 81;//35; // number of magnitude bins
const int Nt = 8967691; // number of time bins
const int healpix_order = 0; // order of healpix map, can be any integer >= 0
const int needlet_order = -1; // maximum order of needlets used, can be any integer >= -1

int healpix_nside = 2^healpix_order;
int Nl = 12*healpix_nside^2;
int Ns = 1;
for (int i = 0; i <= needlet_order; i++) {
    Ns += 12*2^(2*i);
}

const double SingularityPreventer = 10e-13;

const double mut = 3.0;
const double sigmat = 3.0;
const double lg = 0.3;
const double lt = 100.0;

const int PipelineMinVisits = 5; 

#define FILEGAP << ", " << 

Eigen::VectorXd initialisedVector(int n);

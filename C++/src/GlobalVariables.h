#pragma once
#include <vector>
#include "libs/Eigen/Core"
using Eigen::VectorXd;

const int Nm = 81;//35; // number of magnitude bins
const int Nt = 8967691; // number of time bins
const int healpix_order = 0; // order of healpix map, can be any integer >= 0
const int needlet_order = -1; // maximum order of needlets used, can be any integer >= -1

const int healpix_nside = pow(2,healpix_order);
const int Nl = 12*pow(healpix_nside,2);
const int Ns = pow(4,needlet_order+2) - 3;

const int totalRawParams = Nt + Nm*Ns;
const int totalTransformedParams = Nt + Nm*Nl;

const double SingularityPreventer = 10e-13;

const double mut = 3.0;
const double sigmat = 3.0;
const double lm = 0.3;
const double lt = 100.0;

const int PipelineMinVisits = 5; 

#define FILEGAP << ", " << 

Eigen::VectorXd initialisedVector(int n);

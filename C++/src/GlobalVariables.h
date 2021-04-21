#pragma once
#include <vector>
#include "libs/Eigen/Core"
#include <string.h>
using Eigen::VectorXd;

const int Nm = 1;//35; // number of magnitude bins
const int Nt = 2; // number of time bins, coarse, feel free to change
const int TotalScanningTime = 8967691; // number of time bins, must be 8967691, do not change!
const int healpix_order = 0; // order of healpix map, can be any integer >= 0
const int needlet_order = -1; // maximum order of needlets used, can be any integ*needler >= -1

const int healpix_nside = pow(2,healpix_order);
const int Nl = 12*pow(healpix_nside,2);
const int Ns = pow(4,needlet_order+2) - 3;

const int totalRawParams = Nt + Nm*Ns;
const int totalTransformedParams = Nt + Nm*Nl;

const double SingularityPreventer = 1e-18;

const double mut = 0;
const double sigmat = 3;
const double lm = 3;
const double lt = 0.01;

const int SaveSteps = 10000;
const int PipelineMinVisits = 5; 

const bool SaveAllTemps = false;
const std::string TempDirName = "TempPositions";
#define FILEGAP << ", " << 


Eigen::VectorXd initialisedVector(int n);

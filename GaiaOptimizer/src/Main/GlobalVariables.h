#pragma once
#include <vector>
#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include "../libs/Eigen/Core"
#include "../libs/Eigen/Householder"
#include "../libs/Eigen/QR"
#include "../libs/Eigen/LU"
#include <fstream>
#define EIGEN_MPL2_ONLY
#include <string.h>
#include <iostream>
#include <fstream>
#include "../libs/JSL/JSL.h"
#include <iomanip>
using Eigen::VectorXd;

//Optimiser + data properties

const int DataLoadCount = 1e5;	//set to a value > 0, this truncates any datafile readin to that many lines

const std::string TempDirName = "TempPositions";

//temporal, spatial and magnitude resolution

const int Nt = 3e5; // number of time bins, coarse, feel free to change
const int healpix_order = 2; // order of healpix map, can be any integer >= 0
const int needlet_order = 1; // maximum order of needlets used, can be any integ*needler >= -1
const int Nm = 213; // number of magnitude bins
const int magOffset = 0; //offset of loaded files from 0.csv (assuming default file/bin association)


//temporal and magnitude lengthscales + variances

const double sigmat = 5;
const double lm = 3;
const double lt = 100;

//prior and initialisation values
const double xmPrior = -3;
const double xmInitialised = 5;

const double xtPriorNonGap = 5;
const double xtPriorBorderCase =5;
const double xtPriorInsideGap = 5;

const double initialisationBounds = 0.1;

//normal approximation variances

//m-scaling models
const std::vector<double> VariancePopulationFractions = {0.499501,0.407572,0.0929271};
const std::vector<double> VarianceBaselines = {0.0138608,0.00445439,0.0736989};
const std::vector<double> VarianceLinears = {0.003143,0.0034633,0.0211255};
const std::vector<double> VarianceQuadratics = {0,0,0};

//n-gaps-scaling models
//~ const std::vector<double> VariancePopulationFractions = {0.595329,0.4033,0.00131};
//~ const std::vector<double> VarianceBaselines = {0.003599,0.0024822,0.599044};
//~ const std::vector<double> VarianceLinears = {6.95e-5,0.01267,0.03866};
//~ const std::vector<double> VarianceQuadratics = {1.6e-10,3.04e-7,0.000346};




Eigen::VectorXd initialisedVector(int n,std::string loadLocation);


/// OUTPUT STUFF

#define GLOBAL_LOGGING
//#define GLOBAL_DEBUGGING

const int GlobalLoggingLevel = 1;
const int GlobalDebuggingLevel = 8;


void PrintStatus(std::string location);




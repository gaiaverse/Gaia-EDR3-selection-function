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

const int DataLoadCount = 1e6;	//set to a value > 0, this truncates any datafile readin to that many lines

const std::string TempDirName = "TempPositions";

//temporal, spatial and magnitude resolution

const int Nt = 1e6; // number of time bins, coarse, feel free to change
const int healpix_order = 6; // order of healpix map, can be any integer >= 0
const int needlet_order = 5; // maximum order of needlets used, can be any integ*needler >= -1
const int Nm = 213; // number of magnitude bins
const int magOffset = 0; //offset of loaded files from 0.csv (assuming default file/bin association)



const int NVariancePops = 5;
const int hyperOrder = 5;

//temporal and magnitude lengthscales + variances

const double sigmat = 5;
const double lm = 3;
const double lt = 240;

//prior and initialisation values
const double xmPrior = -3;
const double xmInitialised = 5;

const double xtPriorNonGap = 5;
const double xtPriorBorderCase =0;
const double xtPriorInsideGap = -5;

const double initialisationBounds = 0.1;

//normal approximation variances
const bool useMVarianceScaling = true;


Eigen::VectorXd initialisedVector(int n,std::string loadLocation);


/// OUTPUT STUFF

#define GLOBAL_LOGGING
//#define GLOBAL_DEBUGGING

const int GlobalLoggingLevel = 1;
const int GlobalDebuggingLevel = 8;


void PrintStatus(std::string location);




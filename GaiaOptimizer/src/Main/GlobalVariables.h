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

const int DataLoadCount = 5e5;	//set to a value > 0, this truncates any datafile readin to that many lines

const std::string TempDirName = "TempPositions";

//temporal, spatial and magnitude resolution

const int Nt = 1e5; // number of time bins, coarse, feel free to change
const int healpix_order = 1; // order of healpix map, can be any integer >= 0
const int needlet_order = 0; // maximum order of needlets used, can be any integ*needler >= -1
const int Nm = 10; // number of magnitude bins
const int magOffset = 0; //offset of loaded files from 0.csv (assuming default file/bin association)


//temporal and magnitude lengthscales + variances

const double sigmat = 5;
const double lm = 3;
const double lt = 24;

//prior and initialisation values
const double xmPrior = -3;
const double xmInitialised = 5;

const double xtPriorNonGap = 5;
const double xtPriorBorderCase = 0;
const double xtPriorInsideGap = -5;

const double initialisationBounds = 0.1;

//normal approximation variances

const double Population1_Fraction = 0.84415;
const double Population1_Baseline = 0.00215993;
const double Population1_Scaling = 0.00308153;
const double Population2_Baseline = 0.00877248;
const double Population2_Scaling = 0.0192088;


Eigen::VectorXd initialisedVector(int n,std::string loadLocation);


/// OUTPUT STUFF

#define GLOBAL_LOGGING
//#define GLOBAL_DEBUGGING

const int GlobalLoggingLevel = 1;
const int GlobalDebuggingLevel = 8;


void PrintStatus(std::string location);




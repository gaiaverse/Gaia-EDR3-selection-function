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
#include "../GenericFunctions/FileHandler.h"
#include <iomanip>
using Eigen::VectorXd;

const int RootID = 0; //<- declare that process 0 is always Root.

const int Nm = 10; // number of magnitude bins
const int Nt = 10000; // number of time bins, coarse, feel free to change

const int Nt_m = 0;

const int TotalScanningTime = 8967691; // number of time bins, must be 8967691, do not change!
const int healpix_order = 0; // order of healpix map, can be any integer >= 0
const int needlet_order = -1; // maximum order of needlets used, can be any integ*needler >= -1

const int N_SGD_Batches = 128;
const int DataLoadCount = 0;	//set to a value > 0, this truncates any datafile readin to that many lines

const double mum_prior = -2;
const double mum_init = 7;

const double mut_normal = 5;
const double mut_border = 0;
const double mut_gap = -5;

const double mut_mag = -5;
const double sigmat_mag = 3;
const double lt_mag = 5;

const double sigmat = 5;
const double lm = 3;
const double lt = 1;

const double density_alpha = 0.5*log(2.0);
const double density_cut = -3.0;
const double expm_density_cut = exp(-density_cut);

const int healpix_nside = pow(2,healpix_order);
const int Nl = 12*pow(healpix_nside,2);
const int Ns = pow(4,needlet_order+2) - 3;

const int totalRawParams = Nt + Nm*(Ns+Nt_m);
const int totalTransformedParams = Nt + Nm*(Nl+Nt_m);



const int NumberLargerThanMaxObservations = 1024;
const double VerySmallNumber = 1e-310;
const double VerySmallLog = -9999999999;


const double SingularityPreventer = 1e-18;
const int PipelineMinVisits = 5; 

const int SaveSteps = 2;
const bool SaveAllTemps = true;
const std::string TempDirName = "TempPositions";


const bool QuitOnLargeGrad = true;
#define FILEGAP << ", " << 

const double initialisationBounds = 0.1;
Eigen::VectorXd initialisedVector(int n,bool loadIn,std::string loadLocation);


/// OUTPUT STUFF

#define GLOBAL_LOGGING
//#define GLOBAL_DEBUGGING

const int GlobalLoggingLevel = 1;
const int GlobalDebuggingLevel = 8;


void PrintStatus(std::string location);


#define ERROR(exitCode, string){ std::cout << "\n\nCRITICAL ERROR: " << string << "\nTerminating Job\n"; exit(exitCode);}
 
#ifdef GLOBAL_LOGGING
	#define GlobalLog(level, ...){	if (level <= GlobalLoggingLevel){__VA_ARGS__}}
#else
	#define GlobalLog(level, ...){}
#endif	

#ifdef GLOBAL_DEBUGGING
	#define GlobalDebug(level, ...){	if (level <= GlobalDebuggingLevel){__VA_ARGS__}}
#else
	#define GlobalDebug(level, ...)
#endif	
	

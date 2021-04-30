#pragma once
#include <vector>
#include "libs/Eigen/Core"
#include <string.h>
#include <iostream>
#include <fstream>
#include <iomanip>
using Eigen::VectorXd;

const int Nm = 10;//35; // number of magnitude bins
const int Nt = 1; // number of time bins, coarse, feel free to change
const int TotalScanningTime = 8967691; // number of time bins, must be 8967691, do not change!
const int healpix_order = 0; // order of healpix map, can be any integer >= 0
const int needlet_order = -1; // maximum order of needlets used, can be any integ*needler >= -1


const double mut = 5;
const double sigmat = 1;
const double lm = 1;
const double lt = 1;

const int healpix_nside = pow(2,healpix_order);
const int Nl = 12*pow(healpix_nside,2);
const int Ns = pow(4,needlet_order+2) - 3;

const int totalRawParams = Nt + Nm*Ns;
const int totalTransformedParams = Nt + Nm*Nl;





const double SingularityPreventer = 1e-18;

const int SaveSteps = 150;
const int PipelineMinVisits = 5; 

const bool SaveAllTemps = false;
const std::string TempDirName = "TempPositions";



const bool QuitOnLargeGrad = true;
#define FILEGAP << ", " << 

const double initialisationBounds = 0.1;
Eigen::VectorXd initialisedVector(int n);


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
	

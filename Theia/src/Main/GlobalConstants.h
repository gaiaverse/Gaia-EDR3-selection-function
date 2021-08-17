#pragma once
#include "GlobalVariables.h"

/*
	Contains non-alterable constants and compile-time determined variables
	Things in this file should never, ever been changed without knowing what you are doing and why you are doing it
*/


//! Declares that process id=0 is always Root.
const int RootID = 0; 

//!	The maximum number of 10s time bins within the eDR3 catalogue
const int TotalScanningTime = 8967691; 

//! The minimum number of visitations that a star is required to have to enter into the astrometric pipeline
const int PipelineMinVisits = 5; 

//numbers derived from the spatial and temporal resolution, used for number of individual parameters + transforms

//! The number of spatial elements within the EfficiencyVector **raw** vector
const int Nl = 12*pow(2,2*healpix_order);

//! The number of spatial elements within the EfficiencyVector **transformed** vector
const int Ns = pow(4,needlet_order+2) - 3;

const int NHyper = (2 + hyperOrder) * NVariancePops;
const int hyperFractionOffset = (1+hyperOrder)*NVariancePops; 
const int rawNonHyperParams = Nt + Nm*(Ns);
const int transformedNonHyperParams =  Nt + Nm*Nl;

const int totalRawParams = rawNonHyperParams + NHyper;
const int totalTransformedParams = transformedNonHyperParams + NHyper;


const double secondsPerNtStep = 8967691.0 / Nt * 10;
const double lt = 3600*6 * lt_revs / secondsPerNtStep;


//some common re-used numbers
const int NumberLargerThanMaxObservations = 1024;
const double VerySmallNumber = 1e-310;
const double VerySmallLog = -9999999999;
const double SingularityPreventer = 1e-18;


//spatial probability paramaters
const double density_alpha = 0.5*log(2.0);
const double density_cut = -3.0;
const double expm_density_cut = exp(-density_cut);


//some useful macros
#define FILEGAP << ", " << 

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
	
static_assert(hyperOrder % 2 == 0,"hyperOrder must be an even integer");

#pragma once
#include "GlobalVariables.h"

/*
	Contains non-alterable constants and compile-time determined variables
	Things in this file should never, ever been changed without knowing what you are doing and why you are doing it
*/

/* 
 * 
 * Misc Parameters
 * 
*/
//! Declares that process id=0 is always Root.
const int RootID = 0; 

//! The minimum number of visitations that a star is required to have to enter into the astrometric pipeline
const int PipelineMinVisits = 5; 


/* 
 * 
 * Temporal Parameters 
 * 
*/


//!	The maximum number of 10s time bins within the eDR3 catalogue
const int TotalScanningTime = 8967691; 

//! The real-time length (in seconds) of each temporal bin
const double secondsPerNtStep = 8967691.0 / Nt * 10;

//! The temporal correlation length scale (measured in bins)
const double lt = 3600*6 * lt_revs / secondsPerNtStep;



/* 
 * 
 * Length + Index Parameters for Position Vector
 * 
*/

//! The number of spatial elements within the raw vector
const int Nl = 12*pow(2,2*healpix_order);

//! The number of spatial elements within the transformed vector
const int Ns = pow(4,needlet_order+2) - 3;

//! The number of hyperparameters within the vector
const int NHyper = (2 + hyperOrder) * NVariancePops;

//! The index of the first hyperparameter in the raw vector
const int rawNonHyperParams = Nt + Nm*(Ns);

//! The index of the first hyperparameter in the transformed vector
const int transformedNonHyperParams =  Nt + Nm*Nl;

//! The total length of the raw vector
const int totalRawParams = rawNonHyperParams + NHyper;

//!The total length of the transformed vector
const int totalTransformedParams = transformedNonHyperParams + NHyper;

//! The index of the population fraction hyperparameters, relative to the first hyperparameter index
const int hyperFractionOffset = (1+hyperOrder)*NVariancePops; 



/* 
 * 
 * Common Reused Numbers
 * 
*/#

//! The pre-set length of all {p_i} vectors - needs to be larger than max visitations to prevent overflow 
const int NumberLargerThanMaxObservations = 1024;

//! Numbers smaller than this are set to zero (or flagged as errors)
const double VerySmallNumber = 1e-300;
//! Very basic approx for log(0)
const double VerySmallLog = -9999999999;
//! Added to diagonal elements of matrices for stable inversion
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

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

//! The scaling prefactor used to convert x_ml(FOV_1) and x_ml(FOV_2) into a single probability associated with the looking posiiton of Gaia at a given time.
const double spatialAddingPrefactor = 0.5*log(2.0);

//!The truncation point for the elu() function
const double elu_transitionPoint = -3.0;

//! Scaling factor for elu() function below #elu_transitionPoint, equal to exp(elu_transitionPoint)
const double exp_elu_transitionPoint = exp(-elu_transitionPoint);


//! A useful macro used to exit the code immediately + output a given error message into the terminal.
#define ERROR(exitCode, string){ std::cout << "\n\nCRITICAL ERROR: " << string << "\nTerminating Job\n"; exit(exitCode);}
	
static_assert(hyperOrder % 2 == 0,"hyperOrder must be an even integer");

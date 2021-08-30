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





const int Nt =  896769; //!< Number of coarse time bins, can be 0 < Nt < TotalScanningTime 

const int healpix_order = 0; //!< Order of the healpix mapping, can be any integer >= 0. Higher numbers = more pixels on the map.

const int needlet_order = -1; //!< Maximum order of needlets used, can be any integer -1 <= needlet_order <= healpix_order. Higher numbers = finer detail within the map.

const int Nm = 10; //!< Number of magnitude bins

const int NVariancePops = 3; //!< Number of variance populations (hyperparameter)

const int hyperOrder = 4;	 //!< Order of the variance fitting, must be an even number


/*
 * 
 * Priors & Lengthscales
 * 
*/

const double xtPrior = 6; //!< Mean of prior on xt outside of gaps
const double studentNu = 0.5; //!<Student t nu parameter for zt prior


const double sigmat = 3; //!< Standard deviation of zt/xt prior (sort of)

const double lm = 3; //!< Magnitude coupling lengthscale (in mag-bins)

const double lt_revs = 1; //!<Temporal coupling lengthscale for zt (in Gaia revolution periods)

const double xmPrior = 3; //!< Mean of gaussian prior on xml 


const double gapPriorAlpha = 0.5; //!< Beta distribution alpha-value for the prior on the known gaps (in xt space)
const double gapPriorPeak = -6; //!< Beta distribution peak value of xt 
const double gapPriorBeta = gapPriorAlpha * exp(-gapPriorPeak); //!< Given peak value and alpha value, derive appropriate beta parameter

/*
 * 
 * Initialisation
 * 
*/

const double initialisationBounds = 0.3;//!< Unless otherwise specified, the efficiency vector randomly initialises values between +/- this value
const double xmInitialised = 3; //!< The initialisation value of the zeroth-order (whole-sky) spatial mode. 

/*
 * 
 * Data Properties
 * 
*/

//! The subdirectory within the output directory used to save non-converged positions
const std::string TempDirName="TempPositions";

const int DataLoadCount = 1e6;	//!< The maximum number of lines which can be loaded from a file. If the number is 0, or exceeds the number of lines within the file, the entire file is read.
const int magOffset = 0; //!<offset of loaded files from 0.csv (assuming default file/bin association)



//normal approximation variances
enum VarianceScaling {NScaling, MScaling, ActiveNScaling};
const VarianceScaling ScalingMode = ActiveNScaling;
const bool useHyperPrior = true;
const bool ignoreGapObs = false;
Eigen::VectorXd initialisedVector(int n,std::string loadLocation);

//! The 'significant value' - values below this are assumed to be zero and not included in #LogLikelihood::cholesky_w
const double cholesky_tol = 1e-4;

/// OUTPUT STUFF


void PrintStatus(std::string location);
std::vector<bool> AssembleGapList(bool buffered);
extern std::vector<bool> GapList;
extern std::vector<bool> BufferedGapList;


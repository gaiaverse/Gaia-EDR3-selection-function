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





const int Nt = 8967691; //!< Number of coarse time bins, can be 0 < Nt < TotalScanningTime 

const int healpix_order = 1; //!< Order of the healpix mapping, can be any integer >= 0. Higher numbers = more pixels on the map.

const int needlet_order = 0; //!< Maximum order of needlets used, can be any integer -1 <= needlet_order <= healpix_order. Higher numbers = finer detail within the map.

const int Nm = 1; //!< Number of magnitude bins

const int NVariancePops = 3; //!< Number of variance populations (hyperparameter)

const int hyperOrder = 4;	 //!< Order of the variance fitting, must be an even number


/*
 * 
 * Priors & Lengthscales
 * 
*/

const double xtPrior = 8; //!< Mean of prior on xt outside of gaps
const double studentNu = 0.8; //!<Student t nu parameter for zt prior


const double sigmat = 3; //!< Standard deviation of zt/xt prior (sort of)

const double lm = 3; //!< Magnitude coupling lengthscale (in mag-bins)

const double lt_revs = 2; //!<Temporal coupling lengthscale for zt (in Gaia revolution periods)

const double xmPrior = -3; //!< Mean of gaussian prior on xml 


const double gapPriorAlpha = 0.01; //!< Beta distribution alpha-value for the prior on the known gaps (in xt space)
const double gapPriorPeak = -10; //!< Beta distribution peak value of xt 
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

const int DataLoadCount = 0;	//!< The maximum number of lines which can be loaded from a file. If the number is 0, or exceeds the number of lines within the file, the entire file is read.
const int magOffset = 0; //!<offset of loaded files from 0.csv (assuming default file/bin association)


/*
 * 
 * Likelihood Stuff
 * 
*/

//! A list of length #Nt. GapList[i] is true if t[i] corresponds to a time during which our \verbatim embed:rst:inline :ref:`gap-list` \endverbatim indicates Gaia was deactivated.
extern std::vector<bool> GapList;


//! The three kinds of scaling we have incorporated into the Variance Model.
enum VarianceScaling 
{
		//! Scaling parameter is the (fixed) number of visitation of a star
		NScaling, 
		//! Scaling parameter is the mean predicted number of detections (i.e. the sum of probabilities p_t * p_ml)
		MScaling, 
		
		//! Scaling parameter is the mean predicted number of times the star was visited whilst Gaia was turned on (i.e. the sum of probabilities p_t)
		ActiveNScaling
};

//! The chosen kind of scaling
const VarianceScaling ScalingMode = ActiveNScaling;


//! The 'significant value' - values below this are assumed to be zero and not included in #LogLikelihood::cholesky_w
const double choleskyTolerance = 1e-4;

/// OUTPUT STUFF

void PrintStatus(std::string location);



#pragma once
#include "Star.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <mpi.h>
#include <vector>
#include <chrono>
//~ #include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"
#include "../libs/JSL/JSL.h"

/*! A simple structure for holding the name of a \verbatim embed:rst:inline :ref:`stellar data file <star-list>` \endverbatim and the number of stars within it.
 * 
*/
struct FileStarPairs
{
	//! The name of the proposed file
	std::string FileName;
	
	//! The number of stars within it, precalculated from the directory. Needed so that #DataLoadCount stars can be loaded easily from the file.
	int NStars;
};

/*!
  A way of keeping the name of the file and the associated magnitude bin together, forever, in perfect harmony.
 */ 
struct File
{
	std::string Name;
	int Bin;
};



/*!
 * If loaded by the #RootID process, calls the  \verbatim embed:rst:inline :doc:`assignment-creation` \endverbatim protocol (an external Python file), and uses MPI to block all other workers from continuing until the script completes. All workers then GetAssignments() and CalculateBatches(), before using that information to load the assigned stars into the \verbatim embed:rst:inline :doc:`required structures <minibatching>` \endverbatim. 
 * \param ProcessRank The MPI ID of the calling process
 * \param JobSize The number of MPI workers available (used to determine how to distribute files to each worker)
 * \param Data A reference to the container into which the data will be inserter 
 * \param TotalStars A reference to a counter which sums the total number of stars loaded into the worker, and hence into the parallel system
 * \param dataSource The location of the directory in which the stellar data is stored 
 * \param batches The number of minibatches the optimizer uses, and hence a determinant of the structure of the loaded data
 * \returns No explicit returns, but the ``Data`` object becomes populated with Star objects
 * 
*/ 
void LoadData(const int ProcessRank, const int JobSize, std::vector<std::vector<Star>> & Data, int & TotalStars,const std::string dataSource, int batches);

/*!
 * Compares the list of files assigned to the running process by the  \verbatim embed:rst:inline :doc:`assignment-creation` \endverbatim protocol and those within the \verbatim embed:rst:inline :ref:`stellar-directory` \endverbatim. If there is a match (which there **really** should be), it returns the assigned file objects as a list
 * \param id The MPI id of the running process, a number between 0 and NProcesses - 1
 * \param dataSource The directory in which the stellar data is stored
 * \returns A list of the File objects assigned to this process
*/
std::vector<File> GetAssignments(int id,std::string dataSource);

/*!
 * Generates a separation scheme to split the Nstars in each file approximately equally across the minibatches - where there is a rounding error or overflow, the final batch is used to store the excess. 
 * \param Files the output of GetAssignments()
 * \param batches the maximum number of minibatches
*/
void CalculateBatches(std::vector<File> & Files, int batches);

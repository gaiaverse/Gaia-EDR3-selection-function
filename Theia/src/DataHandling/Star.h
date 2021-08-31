#pragma once
#include <vector>
#include <string>
#include <iostream>
// The star class (once fully integrated) will be used to store line-by-line data from the Gaia data files
// There will be ~1.8bn of these objects, so minimising the memory they occupy will be a high priority

/*!
 * A class which represents a single line from one of the \verbatim embed:rst:inline :ref:`stellar data files <star-list>` \endverbatim There are up to 1.8bn of these objects held in memory at any one time, so optimising their memory usage is important.
*/
class Star
{
	public:
		//! The number of times that a star ended up in the Gaia Astrometric Pipeline (a 'detection'). Referred to as ``k`` in our probability work.
		unsigned short int nMeasure;
		
		//! The predicted number of times that Gaia would have seen the star, had it been operating at 100% efficiency for the entire time period, and perfectly followed the commanded scanning law. Referred to as ``n`` in our probability work.
		unsigned short int nVisit;
		
		//! The magnitude bin assigned to the star - a number between 0 and #Nm -1
		unsigned int gBin;
		
		//! A vector of detection times (in units of 10s timebins) that Gaia is predicted to have visited the star. A vector of length #nVisit, with members between 0 and #TotalScanningTime - 1. Annoyingly #TotalScanningTime exceeds the size of a ``short int``.
		std::vector<unsigned int> TimeSeries;

		
		/*!
		 * Constructor function
			\param data A tokenized string representing the line of the input file, assuming the prescription detailed in the \verbatim embed:rst:inline :ref:`data file specification <star-list>` \endverbatim
			\param bin The magnitude bin for #gBin
		*/
		Star(const std::vector<std::string> &data, int bin);
};

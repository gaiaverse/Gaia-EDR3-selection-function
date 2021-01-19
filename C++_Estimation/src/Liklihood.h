#pragma once
#include <vector>
#include "customClasses.h"
class Liklihood
{
	public:
		int ID;
		long double Value;
		std::vector<double> Gradient;
		const std::vector<Star> &Data;
		
		
		Liklihood(const std::vector<Star> & data, int nPoints, int id);
		
		void Calculate(std::vector<double> &position);
};


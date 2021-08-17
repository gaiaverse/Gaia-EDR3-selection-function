#pragma once
#include <vector.h>
#include "../Main/GlobalVariables.h"
#include "../Main/GlobalConstants.h"

class EfficiencyVector
{
	public:
		std::vector<double> Raw;
		std::vector<double> Transformed;
	
	private:
		Transform();
}

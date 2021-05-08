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
#include "../GenericFunctions/FileHandler.h"
#include "../GlobalVariables.h"
#include "../GenericFunctions/timeCodes.h"


struct FileStarPairs
{
	std::string FileName;
	int NStars;
};

void GetAssignments(int id,std::string dataSource);

void LoadData(const int ProcessRank, const int JobSize, std::vector<std::vector<Star>> & Data, int & TotalStars,const std::string dataSource);

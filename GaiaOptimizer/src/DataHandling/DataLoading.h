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

struct FileStarPairs
{
	std::string FileName;
	int NStars;
};
struct File
{
	std::string Name;
	int Bin;
};

std::vector<File> GetAssignments(int id,std::string dataSource);

void LoadData(const int ProcessRank, const int JobSize, std::vector<std::vector<Star>> & Data, int & TotalStars,const std::string dataSource, int batches);

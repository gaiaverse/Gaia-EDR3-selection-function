#pragma once

#include <fstream>
#include <string>
#include <vector>

#include "../Strings/split.h"


#define forLineIn(macroFileName, ...)\
{								\
	do 							\
	{							\
		std::ifstream macroFile(macroFileName);	\
		if (!macroFile.is_open())	\
		{							\
			std::cout << "\n\nERROR: Could not find the file '" << macroFileName << ".\n\nPlease provide a valid filepath.\n\n " << std::endl;	\
			exit(1);				\
		}							\
		std::string FILE_LINE;				\
		while (getline(macroFile,FILE_LINE))	\
		{							\
			__VA_ARGS__				\
		}							\
		macroFile.close();			\
	} while(0);						\
}									\

#define forLineVectorIn(macroFileName, token,...)\
{								\
	do 							\
	{							\
		std::ifstream macroFile(macroFileName);	\
		if (!macroFile.is_open())	\
		{							\
			std::cout << "\n\nERROR: Could not find the file '" << macroFileName << ".\n\nPlease provide a valid filepath.\n\n " << std::endl;	\
			exit(1);				\
		}							\
		std::string FILE_LINE;				\
		while (getline(macroFile,FILE_LINE))	\
		{							\
			std::vector<std::string> FILE_LINE_VECTOR = JSL::split(FILE_LINE,token);	\
			__VA_ARGS__;				\
		}							\
		macroFile.close();			\
	} while(0);						\
}		

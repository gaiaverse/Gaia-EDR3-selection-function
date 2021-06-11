#pragma once
#include <string>;
#include <stdexcept>
#include "../FileIO/FileIO.h";

namespace JSL
{
	template <class T>
	class Argument
	{
		public:
			T Value;
		
			Argument(){};
			Argument(std::string trigger)
			{
				TriggerString = trigger;
				StringCheck();
			}
			Argument(T defaultValue, std::string trigger)
			{
				Value = defaultValue;
				TriggerString = trigger;
				StringCheck();
			}
			Argument(T defaultValue, std::string trigger,const int argc,const char * argv[])
			{
				Value = defaultValue;
				TriggerString = trigger;
				
				for (int i = 1; i < argc-1; ++i)
				{
					Parse(argv[i],argv[i+1]);
				}
				StringCheck();
			}
			
			Argument(T defaultValue, std::string trigger, std::string configFile,const char configDelimiter)
			{
				Value = defaultValue;
				TriggerString = trigger;
				forLineVectorIn(configFile, configDelimiter,
					Parse(FILE_LINE_VECTOR[0].c_str(), FILE_LINE_VECTOR[1].c_str());
				);
				StringCheck();
			}
			
			void Parse(const char * arg,const char * value)
			{
				std::string sArg = arg;
				if (sArg == TriggerString)
				{
					AssignValue(value);
				}
			}
	
			operator T()
			{
				return Value;
			}

		private:
			std::string TriggerString;
			
			void AssignValue(const char * value){};
			
			void StringCheck()
			{
				const std::vector<std::string> ProtectedStrings = {"--help", "-help"};
				
				for (int i = 0; i < ProtectedStrings.size(); ++i)
				{
					if (TriggerString == ProtectedStrings[i])
					{
						throw std::invalid_argument(TriggerString + " is a protected triggername");
					}
				}
			}
	};
	
	template<>
	inline void Argument<int>::AssignValue(const char * value)
	{
		double testDouble = std::stod(value);
		int testInt = std::stoi(value);
		
		if ((double)testInt != testDouble)
		{
			throw std::invalid_argument("Argument passed to " + TriggerString + " was a double, expected an integer");
		}
		Value = std::stoi(value);
	}
	
	
	//specific parsing functions for different types
	template<>
	inline void Argument<double>::AssignValue(const char * value)
	{
		Value = std::stod(value);
	}
	
	template<>
	inline void Argument<std::string>::AssignValue(const char * value)
	{
		Value = (std::string)value;
	}
	
	template<>
	inline void Argument<bool>::AssignValue(const char * value)
	{
		int testValue = std::stoi(value);
		if (testValue != 0 && testValue != 1)
		{
			throw std::runtime_error("Argument passed to " + TriggerString + " was not a bool (0 or 1)");
		}
		Value = (bool)testValue;
	}
	

}

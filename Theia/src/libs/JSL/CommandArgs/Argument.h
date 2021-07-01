#pragma once
#include <string>;
#include <stdexcept>
#include "../FileIO/FileIO.h";

namespace JSL
{
	class ArgumentInterface
	{
		public:
		//~ ArgumentInterface(){};
		virtual void Parse( char * arg, char * value){};
		virtual void ListParse( int argc, char * argv[]){};
		virtual void Configure(std::string configFile, char configDelimiter){};
		protected:
		std::string TriggerString;
	};
	
	template <class T>
	class Argument : public ArgumentInterface
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
			Argument(T defaultValue, std::string trigger, int argc, char * argv[])
			{
				Value = defaultValue;
				TriggerString = trigger;
				
				ListParse(argc,argv);
				StringCheck();
			}
			
			Argument(T defaultValue, std::string trigger, std::string configFile, char configDelimiter)
			{
				Value = defaultValue;
				TriggerString = trigger;
				forLineVectorIn(configFile, configDelimiter,
				
					if (FILE_LINE_VECTOR[0] == TriggerString)
					{
						AssignValue(FILE_LINE_VECTOR[1]);
					}
				);
				StringCheck();
			}
			
			void Configure(std::string configFile, char configDelimiter)
			{
				forLineVectorIn(configFile, configDelimiter,
					if (FILE_LINE_VECTOR.size() > 1 && FILE_LINE_VECTOR[0] == TriggerString)
					{
						AssignValue(FILE_LINE_VECTOR[1]);
					}
				);
			}
			
			void ListParse( int argc, char * argv[])
			{
				for (int i = 1; i < argc-1; ++i)
				{
					Parse(argv[i],argv[i+1]);
				}
			}
			
			void Parse( char * arg, char * value)
			{
				std::string sArg = arg;
				if (sArg == "-" + TriggerString)
				{
					AssignValue(value);
				}
			}
	
			operator T()
			{
				return Value;
			}
			
		private:
			
			
			void AssignValue( char * value){};
			void AssignValue(std::string value)
			{
				char * v = value.data();
				AssignValue(v);
			}
			void StringCheck()
			{
				 std::vector<std::string> ProtectedStrings = {"help"};
				
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
	inline void Argument<int>::AssignValue( char * value)
	{
		double testDouble = std::stod(value);
		int testInt = std::stoi(value);
		
		if ((double)testInt != testDouble)
		{
			throw std::invalid_argument("Argument passed to " + TriggerString + " was a double, expected an integer");
		}
		
		//this method ensures that it can interpret in exponential notation as stod can do that!
		Value = (int)testDouble;
	}
	
	
	//specific parsing functions for different types
	template<>
	inline void Argument<double>::AssignValue( char * value)
	{
		Value = std::stod(value);
	}
	
	template<>
	inline void Argument<std::string>::AssignValue( char * value)
	{
		Value = (std::string)value;
	}
	
	template<>
	inline void Argument<bool>::AssignValue( char * value)
	{
		int testValue = std::stoi(value);
		if (testValue != 0 && testValue != 1)
		{
			throw std::runtime_error("Argument passed to " + TriggerString + " was not a bool (0 or 1)");
		}
		Value = (bool)testValue;
	}
	

}

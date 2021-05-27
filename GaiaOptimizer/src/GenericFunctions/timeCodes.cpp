#include "timeCodes.h"

void printTime()
{
	auto now = std::chrono::system_clock::now();
	std::time_t now_t = std::chrono::system_clock::to_time_t(now);
	std::cout << "Current time is: "<< std::ctime(&now_t);
}
std::string formatDuration(std::chrono::time_point<std::chrono::system_clock> start, std::chrono::time_point<std::chrono::system_clock> end)
{
	std::chrono::duration<double> diff = end - start;
	double seconds = diff.count();
	
	std::vector<std::string> divisions= {"Days", "Hours", "Minutes", "Seconds"};
	std::vector<int> duration = {86400, 3600, 60 , 1};
	//std::vector<int> chunks = std::vector(0,divisions.size());
	
	bool nothingAdded = true;
	std::string output = "";
	for (int i = 0; i < divisions.size(); ++i)
	{
		int v = seconds / duration[i];
		seconds -= v * duration[i];
		
		if (v > 0)
		{
			nothingAdded = false;
			output += std::to_string(v) + " " + divisions[i] + " ";
		}
	}
	if (nothingAdded)
	{
		output = "less than 1 second";
	}
	return output;
}

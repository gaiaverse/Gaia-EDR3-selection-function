#pragma once
#include <chrono>
#include <string>

#include <vector>
#include <iostream>

void printTime();
std::string formatDuration(std::chrono::time_point<std::chrono::system_clock> start, std::chrono::time_point<std::chrono::system_clock> end);

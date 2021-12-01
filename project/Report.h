#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <numeric>
class Report
{

public:
	
	// Method to recall the Report File
	Report();

	// Method to write the Report file 
	void WriteOnASCII(std::ofstream& fOut);

private:

	std::vector<std::string> fileName_ ;
	std::vector<std::string> experimentType_ ;
	std::vector<std::string> warningMessage_ ;
	std::vector<std::string> errorMessage_ ;
};


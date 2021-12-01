#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <numeric>

class Report
{
	
public:

	Report();

	void WriteOnASCII(std::ofstream& fOut);

private:

	std::vector<std::string> fileName_;
	std::vector<std::string> warningMessage_;
	std::vector<std::string> errorMessage_;
	std::vector<std::string> fatalErrorMessage_;
	//std::vector<std::string> fileName_;
};


#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

class Report
{
	
public:

	Report();

	void WriteReportOnASCII(std::ofstream& fOut);

	void WriterReport(const boost::filesystem::path file_name, const boost::filesystem::path output_folder);

	void UpdateErrors(std::string type, std::string error);

	void UpdateFileName(std::string name);

private:

	std::vector<std::string> fileName_;
	std::vector<std::string> warningMessage_;
	std::vector<std::string> genericMessage_;
	std::vector<std::string> errorMessage_;
};


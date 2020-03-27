#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

class DatabaseSpecies
{
public:

	DatabaseSpecies();

	void SetFromXML(const boost::filesystem::path file_name);

	void Summary();

	bool is_active() const { return is_active_; }

	const std::vector<std::string>& names() const { return names_; }

	const std::vector<std::string>& chem_names() const { return chem_names_; }

	const std::vector<std::string>& CAS() const { return CAS_; }

private:

	bool is_active_;

	unsigned int ns_;
	std::vector<std::string> names_;
	std::vector<std::string> chem_names_;
	std::vector<std::string> CAS_;
};


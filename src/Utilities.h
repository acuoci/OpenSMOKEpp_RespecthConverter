/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2020  Alberto Cuoci                                      |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#pragma once

#include <string>
#include <vector>
#include <iostream>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>

class Composition;
class DatabaseSpecies;
struct idtType;

void FatalErrorMessage(const std::string message);

void ReadConstantValueFromXML(boost::property_tree::ptree& subtree, const std::string name, std::vector<double>& value, std::string& units);

void ReadNonConstantValueFromXML(boost::property_tree::ptree& ptree, const std::string name, std::vector<double>& values, std::string& units);

void ReadNonConstantValueFromXML(	boost::property_tree::ptree& ptree, const std::string name, std::vector<Composition>& initial_compositions,
									const std::vector<std::string>& species_in_kinetic_mech, const bool case_sensitive, DatabaseSpecies& database_species);

void ReadProfileFromXML(boost::property_tree::ptree& ptree, const std::string name,
	const std::string name1, std::vector< std::vector<double> >& values1, std::vector<std::string>& units1,
	const std::string name2, std::vector< std::vector<double> >& values2, std::vector<std::string>& units2);

void WriteMixStatusOnASCII(const std::string name, std::ofstream& fOut, const double t, const std::string t_units, const double p, const std::string p_units, const Composition& composition);

void WriteOutputOptionsOnASCII(const std::string name, std::ofstream& fOut, const bool verbose_video, const int steps_video, const bool verbose_file, const int steps_file, const boost::filesystem::path& output_folder_simulation);

void WriteODEParametersOnASCII(const std::string name, std::ofstream& fOut, const double abs_tol, const double rel_tol);

void WriteIgnitionDelayTimesOnASCII(const std::string name, std::ofstream& fOut, const bool is_RCM, const idtType idt);

void WriteParametricAnalysisOnASCII(const std::string name, const std::string type, std::ofstream& fOut, const std::vector<double> values, const std::string units);

void WriteParametricAnalysisOnASCII(const std::string name, const std::string type, std::ofstream& fOut, const std::vector<double> values1, const std::string units1, const std::vector<double> values2, const std::string units2);

void WriteParametricAnalysisOnASCII(const std::string name, const std::string type, std::ofstream& fOut, const std::vector<boost::filesystem::path> file_names);

void WriteProfileOnCSV(const boost::filesystem::path file_name,
	const std::string variable1, const double value1, const std::string unit1,
	const std::string variable2, const double value2, const std::string unit2,
	const std::string name1, const std::vector<double> values1, const std::string units1,
	const std::string name2, const std::vector<double> values2, const std::string units2);

void ForceMonotonicProfiles(std::vector<double>& x, std::vector<double>& y);

void WriteReportFileOnASCII(const boost::filesystem::path file_name, const boost::filesystem::path output_folder, 
	std::vector<boost::filesystem::path> FilesList,
	std::vector<std::string> ErrorsList);

struct idtType
{
	std::string target_;
	std::string type_;
	double amount_;
	std::string units_;
};

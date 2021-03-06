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
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

class DatabaseSpecies;

class Composition
{
public:

	Composition() { };

	void ImportFromXMLTree(boost::property_tree::ptree& ptree, const std::vector<std::string> species_in_kinetic_mech, const bool case_sensitive, DatabaseSpecies& database_species);

	void Set(	const std::vector<std::string> names, const std::vector<double> values, const std::vector<std::string> units,
				const std::vector<std::string> names_chem, const std::vector<std::string> CAS,
				const std::vector<std::string> species_in_kinetic_mech, const bool case_sensitive, DatabaseSpecies& database_species);

	void WriteOnASCIIFile(std::ofstream& fOut) const;

private:

	void ImportFromXMLTree(boost::property_tree::ptree& ptree);

	void CheckForSpeciesNames(const std::vector<std::string> species_in_kinetic_mech, const bool case_sensitive);

	void CheckForSpeciesNames(DatabaseSpecies& database_species);

	void Convert2MoleFractionsAndCheckTheSum();

	void ErrorMessage(const std::string message);

	std::vector<std::string> names_key_;
	std::vector<std::string> names_chem_;
	std::vector<std::string> names_CAS_;
	std::vector<std::string> units_;
	std::vector<double> composition_;
};


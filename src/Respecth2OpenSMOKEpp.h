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

#include <iostream>
#include <numeric>
#include <boost/filesystem.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

#include "Bibliography.h"
#include "Composition.h"
#include "Conversions.h"
#include "Utilities.h"
#include "DatabaseSpecies.h"

class Respecth2OpenSMOKEpp
{
public:

	Respecth2OpenSMOKEpp(	const boost::filesystem::path file_name, 
							const boost::filesystem::path kinetics_folder, 
							const boost::filesystem::path output_folder,
							const std::vector<std::string> species_in_kinetic_mech,
							const bool case_sensitive,
							DatabaseSpecies& database_species);

	void ReadConstantValueFromXML();

	void ReadIdtTypeFromXML();

	void WriteOnASCIIFile(boost::filesystem::path file_name);

	void ErrorMessage(const std::string message);

protected:

	boost::filesystem::path file_name_xml_;
	boost::filesystem::path output_folder_;
	boost::filesystem::path output_folder_simulation_;

	boost::property_tree::ptree ptree_;

	boost::filesystem::path kinetics_folder_;
	std::vector<std::string> species_in_kinetic_mech_;
	bool case_sensitive_;

	std::vector<double> t_values_;
	std::string t_units_;
	bool constant_temperature_;

	std::vector<double> p_values_;
	std::string p_units_;
	bool constant_pressure_;

	std::vector<double> v_values_;
	std::string v_units_;
	bool constant_volume_;

	std::vector<double> phi_values_;
	std::string phi_units_;
	bool constant_equivalenceratio_;

	std::vector<double> m_values_;
	std::string m_units_;
	bool constant_massflowrate_;

	std::vector<double> sl_values_;
	std::string sl_units_;
	bool constant_laminarburningvelocity_;
	
	std::vector<double> tau_values_;
	std::string tau_units_;
	bool constant_residencetime_;

	std::vector<double> dpdt_values_;
	std::string dpdt_units_;
	bool constant_pressurerise_;

	std::vector<Composition> initial_compositions_;
	bool constant_composition_;

	idtType idt_;

protected:

	// Experiment type
	std::string experiment_type_;

	// File data
	std::string file_author_;
	std::string file_doi_;

	// File version
	unsigned int file_version_major_;
	unsigned int file_version_minor_;

	// Respect version
	unsigned int repecth_version_major_;
	unsigned int repecth_version_minor_;

	// Bibliography link
	Bibliography bibliography_;
	
	// Database of species names
	DatabaseSpecies& database_species_;

private:

	void WriteHeaderText(std::ofstream& fOut);
	void WriteMetaData(std::ofstream& fOut);

	virtual void WriteSimulationData(std::ofstream& fOut) = 0;

	virtual void WriteAdditionalFiles() = 0;
};


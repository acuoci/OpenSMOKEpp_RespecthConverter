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

#include "Respecth2OpenSMOKEpp.h"

Respecth2OpenSMOKEpp::Respecth2OpenSMOKEpp(	boost::filesystem::path file_name, 
											const boost::filesystem::path kinetics_folder,
											const boost::filesystem::path output_folder,
											const std::vector<std::string> species_in_kinetic_mech,
											const bool case_sensitive,
											DatabaseSpecies& database_species) :
database_species_(database_species) 
{
	kinetics_folder_ = kinetics_folder;

	species_in_kinetic_mech_ = species_in_kinetic_mech;

	case_sensitive_ = case_sensitive;

	boost::property_tree::read_xml(file_name.string(), ptree_);

	// File author (M)
	file_author_ = ptree_.get<std::string>("experiment.fileAuthor");

	// File DOI (O)
	file_doi_ = ptree_.get<std::string>("experiment.fileDOI", "");

	// File version (O)
	file_version_major_ = ptree_.get<unsigned int>("experiment.fileVersion.major", 0);
	file_version_minor_ = ptree_.get<unsigned int>("experiment.fileVersion.minor", 0);

	// Respect version (M)
	repecth_version_major_ = ptree_.get<unsigned int>("experiment.ReSpecThVersion.major");
	repecth_version_minor_ = ptree_.get<unsigned int>("experiment.ReSpecThVersion.minor");

	// Bibliography link
	bibliography_.ImportFromXMLTree(ptree_);

	// Experiment type (M)
	experiment_type_ = ptree_.get<std::string>("experiment.experimentType");

	// File name XML (complete path)
	file_name_xml_ = file_name;

	// Output folder for the specific experiment
	output_folder_ = output_folder / file_name_xml_.stem();
}

void Respecth2OpenSMOKEpp::WriteOnASCIIFile(boost::filesystem::path file_name)
{
	std::ofstream fOut(file_name.string(), std::ios::out);
	fOut.setf(std::ios::scientific);

	WriteHeaderText(fOut);
	WriteMetaData(fOut);
	WriteSimulationData(fOut);
	fOut.close();

	WriteAdditionalFiles();
}

void Respecth2OpenSMOKEpp::WriteHeaderText(std::ofstream& fOut)
{
	fOut << "//-----------------------------------------------------------------//" << std::endl;
	fOut << "//     ____                    ______ __  __  ____  _  ________    //" << std::endl;
	fOut << "//    / __ \\                  /  ___ |  \\/  |/ __ \\| |/ /  ____|   //" << std::endl;
	fOut << "//   | |  | |_ __   ___ _ __ |  (___ | \\  / | |  | | ' /| |__      //" << std::endl;
	fOut << "//   | |  | | '_ \\ / _ \\ '_ \\ \\___  \\| |\\/| | |  | |  < |  __|     //" << std::endl;
	fOut << "//   | |__| | |_) |  __/ | | |____)  | |  | | |__| | . \\| |____    //" << std::endl;
	fOut << "//    \\____/| .__/ \\___|_| |_|______/|_|  |_|\\____/|_|\\_\\______|   //" << std::endl;
	fOut << "//          | |                                                    //" << std::endl;
	fOut << "//          |_|                                                    //" << std::endl;
	fOut << "//                                                                 //" << std::endl;
	fOut << "//              http://www.opensmokepp.polimi.it/                  //" << std::endl;
	fOut << "//             http://creckmodeling.chem.polimi.it/                //" << std::endl;
	fOut << "//-----------------------------------------------------------------//" << std::endl;
	fOut << std::endl;
}

void Respecth2OpenSMOKEpp::WriteMetaData(std::ofstream& fOut)
{
	bibliography_.WriteOnASCII(fOut);
}

void Respecth2OpenSMOKEpp::ErrorMessage(const std::string message)
{
	std::cout << "Fatal Error Message in conversion" << message << std::endl;
	std::cout << " * Respecth file:   " << file_name_xml_.string() << std::endl;
	std::cout << " * Experiment type: " << experiment_type_ << std::endl;
	std::cout << " * Error message:   " << message << std::endl;
	std::cout << "Press enter to exit... ";
	getchar();
	exit(-1);
}

void Respecth2OpenSMOKEpp::ReadConstantValueFromXML()
{
	// Recognize the available data
	constant_temperature_ = false;
	constant_pressure_ = false;
	constant_composition_ = false;
	constant_residencetime_ = false;
	constant_equivalenceratio_ = false;
	constant_massflowrate_ = false;
	constant_pressurerise_ = false;
	constant_volume_ = false;

	// Check for constant variables
	{
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree_.get_child("experiment.commonProperties"))
		{
			boost::property_tree::ptree subtree = node.second;

			if (node.first == "property")
			{
				if (subtree.get<std::string>("<xmlattr>.name") == "temperature")
					constant_temperature_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "pressure")
					constant_pressure_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "initial composition")
					constant_composition_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "residence time")
					constant_residencetime_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "volume")
					constant_volume_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "equivalence ratio")
					constant_equivalenceratio_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "flow rate")
					constant_massflowrate_ = true;

				if (subtree.get<std::string>("<xmlattr>.name") == "pressure rise")
					constant_pressurerise_ = true;
			}
		}
	}

	// Read constant properties
	{
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree_.get_child("experiment.commonProperties"))
		{
			boost::property_tree::ptree subtree = node.second;

			if (node.first == "property")
			{
				// Temperature
				if (constant_temperature_ == true)
					::ReadConstantValueFromXML(subtree, "temperature", t_values_, t_units_);

				// Pressure
				if (constant_pressure_ == true)
					::ReadConstantValueFromXML(subtree, "pressure", p_values_, p_units_);

				// Residence time
				if (constant_residencetime_ == true)
					::ReadConstantValueFromXML(subtree, "residence time", tau_values_, tau_units_);

				// Volume
				if (constant_volume_ == true)
					::ReadConstantValueFromXML(subtree, "volume", v_values_, v_units_);

				// Equivalence ratio
				if (constant_equivalenceratio_ == true)
					::ReadConstantValueFromXML(subtree, "equivalence ratio", phi_values_, phi_units_);

				// Mass flow rate
				if (constant_massflowrate_ == true)
					::ReadConstantValueFromXML(subtree, "flow rate", m_values_, m_units_);

				// Pressure rise
				if (constant_pressurerise_ == true)
					::ReadConstantValueFromXML(subtree, "pressure rise", dpdt_values_, dpdt_units_);

				// Initial composition
				if (constant_composition_ == true)
				{
					if (subtree.get<std::string>("<xmlattr>.name") == "initial composition")
					{
						initial_compositions_.resize(1);
						initial_compositions_[0].ImportFromXMLTree(subtree, species_in_kinetic_mech_, case_sensitive_, database_species_);
					}
				}
			}
		}
	}
}

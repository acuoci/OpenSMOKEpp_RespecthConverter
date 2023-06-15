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

#include <boost/program_options.hpp>

//Utilities
#include "Utilities.h"

// Converters
#include "Respecth2OpenSMOKEpp_JetStirredReactor.h"
#include "Respecth2OpenSMOKEpp_LaminarBurningVelocity.h"
#include "Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation.h"
#include "Respecth2OpenSMOKEpp_ConcentrationTimeProfile.h"
#include "Respecth2OpenSMOKEpp_OutletConcentration.h"
#include "Respecth2OpenSMOKEpp_IgnitionDelay.h"


int ExitStatus = -1;
std::vector<std::string> ErrorList = {};
int IndexExperimentWithError;

int main(int argc, char** argv)
{

	std::string input_file_name_;
	std::string output_file_name_;

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_RespecthConverter");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the XML file in respecth format to be converted")
			("output", po::value<std::string>(), "name of the OpenSMOKE input file generated from the XML file");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return 0;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("output"))
				output_file_name_ = vm["output"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return -1;
		}
	}

	boost::filesystem::path path_input_file_name_ = input_file_name_;
	boost::filesystem::path path_output_file_name_ = output_file_name_;

	// if (!path_input_file_name_.is_absolute())
	//	path_input_file_name_ = boost::filesystem::current_path() / path_input_file_name_;

	if (!boost::filesystem::is_regular_file(path_input_file_name_))
	{
		std::cout << "The provided input file is not a valid" << std::endl;
		return -1;
	}

	DatabaseSpecies database_species;
	// database_species.SetFromXML("");

	// Convert files
	std::string apparatus_kind;
	std::string experiment_type;

	boost::property_tree::ptree ptree;
	boost::property_tree::read_xml(path_input_file_name_.string(), ptree);
	try
	{
		apparatus_kind = ptree.get<std::string>("experiment.apparatus.kind");
	}
	catch (const boost::property_tree::ptree_error& e)
	{
		std::cout << "Fatal error: " << e.what() << std::endl << std::endl;
	}
		
	try
	{
		experiment_type = ptree.get<std::string>("experiment.experimentType");
	}
	catch (const boost::property_tree::ptree_error& e)
	{
		std::cout << "Fatal error: " << e.what() << std::endl << std::endl;
	}

	if (experiment_type == "jet stirred reactor measurement")
	{
		Respecth2OpenSMOKEpp_JetStirredReactor reactor(path_input_file_name_, path_output_file_name_, database_species);
		reactor.WriteOnASCIIFile((path_output_file_name_.filename().string() + ".dic"));
	}
	else if (experiment_type == "laminar burning velocity measurement")
	{
		Respecth2OpenSMOKEpp_LaminarBurningVelocity reactor(path_input_file_name_, path_output_file_name_, database_species);
		reactor.WriteOnASCIIFile((path_output_file_name_.filename().string() + ".dic"));
	}
	else if (experiment_type == "burner stabilized flame speciation measurement")
	{
		Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation reactor(path_input_file_name_, path_output_file_name_, database_species);
		reactor.WriteOnASCIIFile((path_output_file_name_.filename().string() + ".dic"));
	}
	else if (experiment_type == "concentration time profile measurement")
	{
		Respecth2OpenSMOKEpp_ConcentrationTimeProfile reactor(path_input_file_name_, path_output_file_name_, database_species);
		reactor.WriteOnASCIIFile((path_output_file_name_.filename().string() + ".dic"));
	}
	else if (experiment_type == "outlet concentration measurement")
	{
		Respecth2OpenSMOKEpp_OutletConcentration reactor(path_input_file_name_, path_output_file_name_, database_species);
		reactor.WriteOnASCIIFile((path_output_file_name_.filename().string() + ".dic"));
	}	
	else if (experiment_type == "ignition delay measurement")
	{
		Respecth2OpenSMOKEpp_IgnitionDelay reactor(path_input_file_name_, path_output_file_name_, database_species);
		reactor.WriteOnASCIIFile((path_output_file_name_.filename().string() + ".dic"));
	}
	else
	{
		std::cout << "Unknown experiment type: " << experiment_type << std::endl;
		return -1;
	}

	return 0;
}
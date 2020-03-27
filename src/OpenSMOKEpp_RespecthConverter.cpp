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

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// Grammar
#include "Grammar_RespecthConverter.h"

// Converters
#include "Respecth2OpenSMOKEpp_JetStirredReactor.h"
#include "Respecth2OpenSMOKEpp_LaminarBurningVelocity.h"
#include "Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation.h"
#include "Respecth2OpenSMOKEpp_ConcentrationTimeProfile.h"
#include "Respecth2OpenSMOKEpp_OutletConcentration.h"
#include "Respecth2OpenSMOKEpp_IgnitionDelay.h"

int main(int argc, char** argv)
{
	boost::filesystem::path executable_file = OpenSMOKE::GetExecutableFileName(argv);
	boost::filesystem::path executable_folder = executable_file.parent_path();

	OpenSMOKE::OpenSMOKE_logo("OpenSMOKEpp_RespecthConverter", "Alberto Cuoci (alberto.cuoci@polimi.it)");

	std::string input_file_name_ = "input.dic";
	std::string main_dictionary_name_ = "RespecthConverter";
	unsigned int number_threads = 1;

	// Program options from command line
	{
		namespace po = boost::program_options;
		po::options_description description("Options for the OpenSMOKEpp_RespecthConverter");
		description.add_options()
			("help", "print help messages")
			("input", po::value<std::string>(), "name of the file containing the main dictionary (default \"input.dic\")")
			("dictionary", po::value<std::string>(), "name of the main dictionary to be used (default \"RespecthConverter\")");

		po::variables_map vm;
		try
		{
			po::store(po::parse_command_line(argc, argv, description), vm); // can throw 

			if (vm.count("help"))
			{
				std::cout << "Basic Command Line Parameters" << std::endl;
				std::cout << description << std::endl;
				return OPENSMOKE_SUCCESSFULL_EXIT;
			}

			if (vm.count("input"))
				input_file_name_ = vm["input"].as<std::string>();

			if (vm.count("dictionary"))
				main_dictionary_name_ = vm["dictionary"].as<std::string>();

			po::notify(vm); // throws on error, so do after help in case  there are any problems 
		}
		catch (po::error& e)
		{
			std::cerr << "Fatal error: " << e.what() << std::endl << std::endl;
			std::cerr << description << std::endl;
			return OPENSMOKE_FATAL_ERROR_EXIT;
		}
	}

	// Defines the grammar rules
	OpenSMOKE::Grammar_RespecthConverter grammar_respecthconverter;

	// Define the dictionaries
	OpenSMOKE::OpenSMOKE_DictionaryManager dictionaries;
	dictionaries.ReadDictionariesFromFile(input_file_name_);
	dictionaries(main_dictionary_name_).SetGrammar(grammar_respecthconverter);

	// Kinetics folder
	boost::filesystem::path path_kinetics_folder_remote;
	std::vector<std::string> species_in_kinetic_mech;
	if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolder") == true)
	{
		boost::filesystem::path path_kinetics_folder;

		dictionaries(main_dictionary_name_).ReadPath("@KineticsFolder", path_kinetics_folder);
		OpenSMOKE::CheckKineticsFolder(path_kinetics_folder);

		// Read names of species
		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml( (path_kinetics_folder / "kinetics.xml").string(), ptree);

		// Species in the kinetic mechanism
		const unsigned int ns = ptree.get<unsigned int>("opensmoke.NumberOfSpecies");
		species_in_kinetic_mech.resize(ns);

		std::stringstream stream;
		stream.str(ptree.get< std::string >("opensmoke.NamesOfSpecies"));
		for (unsigned int i = 0; i < ns; i++)
			stream >> species_in_kinetic_mech[i];

		path_kinetics_folder_remote = path_kinetics_folder;
	}

	if (dictionaries(main_dictionary_name_).CheckOption("@KineticsFolderRemote") == true)
		dictionaries(main_dictionary_name_).ReadPath("@KineticsFolderRemote", path_kinetics_folder_remote);

	boost::filesystem::path path_output_folder_remote;
	if (dictionaries(main_dictionary_name_).CheckOption("@OutputFolderRemote") == true)
		dictionaries(main_dictionary_name_).ReadPath("@OutputFolderRemote", path_output_folder_remote);

	bool case_sensitive = false;
	if (dictionaries(main_dictionary_name_).CheckOption("@CaseSensitiveSpecies") == true)
		dictionaries(main_dictionary_name_).ReadBool("@CaseSensitiveSpecies", case_sensitive);

	DatabaseSpecies database_species;
	if (dictionaries(main_dictionary_name_).CheckOption("@DatabaseSpecies") == true)
	{
		boost::filesystem::path path_database_species;
		dictionaries(main_dictionary_name_).ReadPath("@DatabaseSpecies", path_database_species);
		database_species.SetFromXML(path_database_species);
	}

	// Read list of xml files to be converted
	std::vector<boost::filesystem::path> list_xml_files;
	if (dictionaries(main_dictionary_name_).CheckOption("@InputFolder") == true)
	{
		boost::filesystem::path input_folder;
		dictionaries(main_dictionary_name_).ReadPath("@InputFolder", input_folder);
		if (input_folder.is_absolute() == false)
			input_folder = boost::filesystem::current_path() / input_folder;

		if (!boost::filesystem::exists(input_folder))
		{
			OpenSMOKE::FatalErrorMessage("The provided @InputFolder is not a directory, but a file");
		}
		else
		{
			if (boost::filesystem::is_regular_file(input_folder))
				OpenSMOKE::FatalErrorMessage("The provided @InputFolder is not a directory, but a file");
			else if (!boost::filesystem::is_directory(input_folder))
				OpenSMOKE::FatalErrorMessage("The provided @InputFolder exists, but is neither a regular file nor a directory");

			boost::filesystem::recursive_directory_iterator it(input_folder);
			boost::filesystem::recursive_directory_iterator endit;

			while (it != endit)
			{
				if (boost::filesystem::is_regular_file(*it) && it->path().extension() == ".xml")
					list_xml_files.push_back(it->path());
				++it;
			}
		}
	}

	// Convert files
	std::vector<std::string> apparatus_kind(list_xml_files.size());
	std::vector<std::string> experiment_type(list_xml_files.size());
	for (unsigned int j = 0; j < list_xml_files.size(); j++)
	{
		std::cout << list_xml_files[j].string() << std::endl;

		boost::property_tree::ptree ptree;
		boost::property_tree::read_xml(list_xml_files[j].string(), ptree);

		apparatus_kind[j] = ptree.get<std::string>("experiment.apparatus.kind");
		experiment_type[j] = ptree.get<std::string>("experiment.experimentType");

		std::cout << j+1 << "/" << list_xml_files.size() << " " << apparatus_kind[j] << " " << experiment_type[j] << std::endl;
	}

	for (unsigned int j = 0; j < list_xml_files.size(); j++)
	{
		std::cout << "Converting file: " << list_xml_files[j].filename().string() << std::endl;
		std::cout << apparatus_kind[j] << " " << experiment_type[j] << std::endl;

		if (experiment_type[j] == "jet stirred reactor measurement")
		{
			Respecth2OpenSMOKEpp_JetStirredReactor reactor(list_xml_files[j], path_kinetics_folder_remote, path_output_folder_remote, species_in_kinetic_mech, case_sensitive, database_species);
			reactor.WriteOnASCIIFile( (list_xml_files[j].filename().string() + ".dic" ) );
		}

		else if (experiment_type[j] == "laminar burning velocity measurement")
		{
			Respecth2OpenSMOKEpp_LaminarBurningVelocity reactor(list_xml_files[j], path_kinetics_folder_remote, path_output_folder_remote, species_in_kinetic_mech, case_sensitive, database_species);
			reactor.WriteOnASCIIFile((list_xml_files[j].filename().string() + ".dic"));
		}

		else if (experiment_type[j] == "burner stabilized flame speciation measurement")
		{
			Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation reactor(list_xml_files[j], path_kinetics_folder_remote, path_output_folder_remote, species_in_kinetic_mech, case_sensitive, database_species);
			reactor.WriteOnASCIIFile((list_xml_files[j].filename().string() + ".dic"));
		}

		else if (experiment_type[j] == "concentration time profile measurement")
		{
			Respecth2OpenSMOKEpp_ConcentrationTimeProfile reactor(list_xml_files[j], path_kinetics_folder_remote, path_output_folder_remote, species_in_kinetic_mech, case_sensitive, database_species);
			reactor.WriteOnASCIIFile((list_xml_files[j].filename().string() + ".dic"));
		}

		else if (experiment_type[j] == "outlet concentration measurement")
		{
			Respecth2OpenSMOKEpp_OutletConcentration reactor(list_xml_files[j], path_kinetics_folder_remote, path_output_folder_remote, species_in_kinetic_mech, case_sensitive, database_species);
			reactor.WriteOnASCIIFile((list_xml_files[j].filename().string() + ".dic"));
		}

		else if (experiment_type[j] == "ignition delay measurement")
		{
			Respecth2OpenSMOKEpp_IgnitionDelay reactor(list_xml_files[j], path_kinetics_folder_remote, path_output_folder_remote, species_in_kinetic_mech, case_sensitive, database_species);
			reactor.WriteOnASCIIFile((list_xml_files[j].filename().string() + ".dic"));
		}
	}
}
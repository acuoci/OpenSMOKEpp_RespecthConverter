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

#include "Composition.h"
#include "DatabaseSpecies.h"
#include <algorithm>
#include <iterator>
#include <cmath>

void Composition::ImportFromXMLTree(boost::property_tree::ptree& ptree)
{
	BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child(""))
	{
		boost::property_tree::ptree subtree = node.second;

		if (node.first == "component")
		{
			// Name of species
			names_key_.push_back(subtree.get<std::string>("speciesLink.<xmlattr>.preferredKey"));

			// Composition (together with units)
			units_.push_back(subtree.get<std::string>("amount.<xmlattr>.units"));
			composition_.push_back(subtree.get<double>("amount"));

			// Optional names
			names_chem_.push_back(subtree.get<std::string>("speciesLink.<xmlattr>.chemName", ""));
			names_CAS_.push_back(subtree.get<std::string>("speciesLink.<xmlattr>.CAS", ""));
		}
	}
}

void Composition::ImportFromXMLTree(boost::property_tree::ptree& ptree, const std::vector<std::string> species_in_kinetic_mech, const bool case_sensitive, DatabaseSpecies& database_species)
{
	// Import species
	ImportFromXMLTree(ptree);
	
	// Check for species in the database
	CheckForSpeciesNames(database_species);

	// Check species names
	CheckForSpeciesNames(species_in_kinetic_mech, case_sensitive);
	
	// Convert to mole fractions and check the sum
	Convert2MoleFractionsAndCheckTheSum();
}

void Composition::Set(	const std::vector<std::string> names, const std::vector<double> values, const std::vector<std::string> units,
						const std::vector<std::string> names_chem, const std::vector<std::string> CAS,
						const std::vector<std::string> species_in_kinetic_mech, const bool case_sensitive, DatabaseSpecies& database_species)
{
	names_key_ = names;
	composition_ = values;
	units_ = units;
	names_chem_ = names_chem;
	names_CAS_ = CAS;

	// Check for species in the database
	CheckForSpeciesNames(database_species);

	// Check species names
	CheckForSpeciesNames(species_in_kinetic_mech, case_sensitive);

	// Convert to mole fractions and check the sum
	Convert2MoleFractionsAndCheckTheSum();
}

void Composition::CheckForSpeciesNames(const std::vector<std::string> species_in_kinetic_mech, const bool case_sensitive)
{
	// Check for the existence of species (case sensitive)
	if (case_sensitive == true)
	{
		for (unsigned int i = 0; i < names_key_.size(); i++)
		{
			std::vector<std::string>::const_iterator it = std::find(species_in_kinetic_mech.begin(), species_in_kinetic_mech.end(), names_key_[i]);
			if (it == species_in_kinetic_mech.end())
				ErrorMessage("Case sensitive check: Species " + names_key_[i] + " is not available in the kinetic mechanism.");
		}
	}
	// Check for the existence of species (case insensitive)
	else
	{
		std::vector<std::string> species_in_kinetic_mech_upper(species_in_kinetic_mech.size());
		for (unsigned int i = 0; i < species_in_kinetic_mech.size(); i++)
			species_in_kinetic_mech_upper[i] = boost::to_upper_copy(species_in_kinetic_mech[i]);

		for (unsigned int i = 0; i < names_key_.size(); i++)
		{
			const std::string names_key_upper = boost::to_upper_copy(names_key_[i]);
			std::vector<std::string>::const_iterator it = std::find(species_in_kinetic_mech_upper.begin(), species_in_kinetic_mech_upper.end(), names_key_upper);
			if (it != species_in_kinetic_mech_upper.end())
			{
				std::vector<std::string>::const_iterator first = species_in_kinetic_mech_upper.begin();
				names_key_[i] = species_in_kinetic_mech[std::distance(first, it)];
			}
			else
				ErrorMessage("Case unsensitive check: Species " + names_key_[i] + " is not available in the kinetic mechanism.");
		}
	}
}

void Composition::CheckForSpeciesNames(DatabaseSpecies& database_species)
{
	if (database_species.is_active())
	{
		for (unsigned int i = 0; i < names_key_.size(); i++)
		{
			if (names_CAS_[i] != "")
			{
				std::vector<std::string>::const_iterator it = std::find(database_species.CAS().begin(), database_species.CAS().end(), names_CAS_[i]);
				if (it != database_species.CAS().end())
				{
					std::cout << "Before: " << names_key_[i] << " " << database_species.CAS()[std::distance(database_species.CAS().begin(), it)] << std::endl;
					names_key_[i] = database_species.names()[std::distance(database_species.CAS().begin(), it)];
					//std::cout << "After: " << names_key_[i] << std::endl;
				}
				else
				{
					if (names_chem_[i] != "")
					{
						std::vector<std::string>::const_iterator it = std::find(database_species.chem_names().begin(), database_species.chem_names().end(), names_chem_[i]);
						if (it != database_species.chem_names().end())
						{
							//std::cout << "Before: " << names_key_[i] << " " << database_species.chem_names()[std::distance(database_species.chem_names().begin(), it)] << std::endl;
							names_key_[i] = database_species.names()[std::distance(database_species.chem_names().begin(), it)];
							//std::cout << "After: " << names_key_[i] << std::endl;
						}
					}
				}
			}
		}
	}
}

void Composition::Convert2MoleFractionsAndCheckTheSum()
{	
	
	bool composition_in_concentration = false;

	for (unsigned int i = 0; i < composition_.size(); i++)
	{
		if (units_[i] == "mole fraction")
		{
			// do nothing
		}
		else if (units_[i] == "percent")
		{
			composition_[i] /= 100.;
			units_[i] = "mole fraction";
		}
		else if (units_[i] == "ppm")
		{
			composition_[i] *= 1.e6;
			units_[i] = "mole fraction";
		}
		else if (units_[i] == "ppb")
		{
			composition_[i] *= 1.e9;
			units_[i] = "mole fraction";
		}
		else if (units_[i] == "mol/cm3")
		{
			composition_in_concentration = true;
		}
		else
		{
			ErrorMessage("Unknown units for composition: " + units_[i] + ". Available units: mole fraction | percent | ppm | ppb | mol/cm3");
		}
	}

	// In case of composition given in terms of concentration
	if (composition_in_concentration == true)
	{
		// Check if concentration is provided for all the species (and in the same units)
		for (unsigned int i = 0; i < composition_.size(); i++)
		{
			if (units_[i] != "mol/cm3")
				ErrorMessage("If composition is given in terms of concentration, this must be done for all the species.");
			units_[i] = "mole fraction";
		}

		const double ctot = std::accumulate(composition_.begin(), composition_.end(), 0.);
		for (unsigned int i = 0; i < composition_.size(); i++)
			composition_[i] /= ctot;
	}

	// Sum
	const double sum_threshold = 1.0001e-4;
	const double sum = std::accumulate(composition_.begin(), composition_.end(), 0.);
	if (std::fabs(sum - 1.) > sum_threshold) {
		ErrorMessage("Sum is not equal to 1: " + std::to_string(sum));
	}

	// Normalization
	for (unsigned int i = 0; i < composition_.size(); i++)
		composition_[i] /= sum;
}

void Composition::WriteOnASCIIFile(std::ofstream& fOut) const
{
	if (units_[0] == "mole fraction")
		fOut << "        @MoleFractions    ";

	for (unsigned int i = 0; i < composition_.size(); i++)
		fOut << names_key_[i] << " " << composition_[i] << " ";
	fOut << ";" << std::endl;
}

void Composition::ErrorMessage(const std::string message)
{
	std::cout << "Fatal Error Message: " << message << std::endl;
	std::cout << "Press enter to exit... ";
	getchar();
	exit(-1);
	
}

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

#include "Utilities.h"

#include "Conversions.h"
#include "Composition.h"
#include "DatabaseSpecies.h"

void FatalErrorMessage(const std::string message)
{
	std::cout << "Fatal Error Message: " << message << std::endl;
	std::cout << "Press enter to exit... ";
	getchar();
	exit(-1);
}

void ReadConstantValueFromXML(boost::property_tree::ptree& subtree, const std::string name, std::vector<double>& values, std::string& units)
{
	if (subtree.get<std::string>("<xmlattr>.name") == name)
	{
		double value = subtree.get<double>("value");
		units = subtree.get<std::string>("<xmlattr>.units");
		CheckAndConvertUnits(name, value, units);
		values.push_back(value);
	}
}

void ReadNonConstantValueFromXML(boost::property_tree::ptree& ptree, const std::string name, std::vector<double>& values, std::string& units)
{
	std::string id = "n.a.";

	units = "n.a.";
	BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("experiment.dataGroup"))
	{
		boost::property_tree::ptree subtree = node.second;

		if (node.first == "property")
		{
			if (subtree.get<std::string>("<xmlattr>.name") == name)
			{
				id = subtree.get<std::string>("<xmlattr>.id");
				units = subtree.get<std::string>("<xmlattr>.units");
			}
		}
	}

	// Read the profile, only if it exists
	if (id != "n.a.")
	{
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("experiment.dataGroup"))
		{
			boost::property_tree::ptree subtree = node.second;

			if (node.first == "dataPoint")
				values.push_back(subtree.get<double>(id));
		}
	}

	CheckAndConvertUnits(name, values, units);
}

void ReadProfileFromXML(boost::property_tree::ptree& ptree, const std::string name, 
						const std::string name1, std::vector< std::vector<double> >& values1, std::vector<std::string>& units1,
						const std::string name2, std::vector< std::vector<double> >& values2, std::vector<std::string>& units2)
{
	std::vector<std::string> id1;
	std::vector<std::string> id2;
	
	BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("experiment"))
	{
		boost::property_tree::ptree subtree = node.second;

		if (node.first == "dataGroup")
		{
			const std::string label = subtree.get<std::string>("<xmlattr>.label","");

			if (label == name)
			{
				BOOST_FOREACH(boost::property_tree::ptree::value_type const& subnode, subtree)
				{
					boost::property_tree::ptree subsubtree = subnode.second;

					if (subnode.first == "property")
					{
						if (subsubtree.get<std::string>("<xmlattr>.name") == name1)
						{
							id1.push_back( subsubtree.get<std::string>("<xmlattr>.id") );
							units1.push_back( subsubtree.get<std::string>("<xmlattr>.units") );
						}
						if (subsubtree.get<std::string>("<xmlattr>.name") == name2)
						{
							id2.push_back( subsubtree.get<std::string>("<xmlattr>.id") );
							units2.push_back( subsubtree.get<std::string>("<xmlattr>.units") );
						}
					}
				}
			}
		}
	}

	if (id1.size() != 0)
	{
		values1.resize(id1.size());
		values2.resize(id2.size());

		unsigned int i = 0;
		BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("experiment"))
		{
			boost::property_tree::ptree subtree = node.second;

			if (node.first == "dataGroup")
			{
				const std::string label = subtree.get<std::string>("<xmlattr>.label", "");

				if (label == name)
				{
					BOOST_FOREACH(boost::property_tree::ptree::value_type const& subnode, subtree)
					{
						boost::property_tree::ptree subsubtree = subnode.second;

						if (subnode.first == "dataPoint")
						{
							values1[i].push_back(subsubtree.get<double>(id1[i]));
							values2[i].push_back(subsubtree.get<double>(id2[i]));
						}
					}

					i++;
				}
			}
		}

		for (unsigned int i = 0; i < id1.size(); i++)
			CheckAndConvertUnits(name1, values1[i], units1[i]);

		for (unsigned int i = 0; i < id2.size(); i++)
			CheckAndConvertUnits(name2, values2[i], units2[i]);

	}
}

void ReadNonConstantValueFromXML(boost::property_tree::ptree& ptree, const std::string name, std::vector<Composition>& initial_compositions,
								const std::vector<std::string>& species_in_kinetic_mech, const bool case_sensitive, DatabaseSpecies& database_species)
{
	std::vector<std::string> id;
	std::vector<std::string> composition_names;
	std::vector<std::string> composition_units;
	std::vector<std::string> composition_names_chem;
	std::vector<std::string> composition_CAS;

	BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("experiment.dataGroup"))
	{
		boost::property_tree::ptree subtree = node.second;

		if (node.first == "property")
		{
			if (subtree.get<std::string>("<xmlattr>.name") == "composition")
			{
				id.push_back(subtree.get<std::string>("<xmlattr>.id"));

				composition_units.push_back(subtree.get<std::string>("<xmlattr>.units"));
				composition_names.push_back(subtree.get<std::string>("speciesLink.<xmlattr>.preferredKey"));
				composition_names_chem.push_back(subtree.get<std::string>("speciesLink.<xmlattr>.chemName", ""));
				composition_CAS.push_back(subtree.get<std::string>("speciesLink.<xmlattr>.CAS", ""));
			}
		}
	}

	BOOST_FOREACH(boost::property_tree::ptree::value_type const& node, ptree.get_child("experiment.dataGroup"))
	{
		std::vector<double> composition_values;
		boost::property_tree::ptree subtree = node.second;

		if (node.first == "dataPoint")
		{
			for (unsigned int j = 0; j < id.size(); j++)
				composition_values.push_back(subtree.get<double>(id[j]));

			Composition tmp; tmp.Set(composition_names, composition_values, composition_units,
				composition_names_chem, composition_CAS,
				species_in_kinetic_mech, case_sensitive, database_species);

			initial_compositions.push_back(tmp);
		}
	}
}

void WriteMixStatusOnASCII(const std::string name, std::ofstream& fOut, const double t, const std::string t_units, const double p, const std::string p_units, const Composition& composition)
{
	std::cout << "   - mix status" << std::endl;

	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @Temperature      " << t << " " << t_units << " ;" << std::endl;
	fOut << "        @Pressure         " << p << " " << p_units << " ;" << std::endl;
	composition.WriteOnASCIIFile(fOut);
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteOutputOptionsOnASCII(const std::string name, std::ofstream& fOut, const int steps_file, const boost::filesystem::path& output_folder)
{
	std::cout << "   - output options" << std::endl;

	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @StepsVideo       1000;" << std::endl;
	fOut << "        @StepsFile        " << steps_file << ";" << std::endl;
	fOut << "        @VerboseVideo     false;" << std::endl;
	fOut << "        @VerboseASCIIFile true;" << std::endl;
	fOut << "        @OutputFolder     " << output_folder.string() << ";" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteODEParametersOnASCII(const std::string name, std::ofstream& fOut, const double abs_tol, const double rel_tol)
{
	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @OdeSolver         OpenSMOKE;" << std::endl;
	fOut << "        @AbsoluteTolerance " << abs_tol << ";" << std::endl;
	fOut << "        @RelativeTolerance " << rel_tol << ";" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteParametricAnalysisOnASCII(const std::string name, const std::string type, std::ofstream& fOut, const std::vector<double> values, const std::string units)
{
	std::cout << "   - parametric analysis" << std::endl;

	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @Type          " << type << ";" << std::endl;
	fOut << "        @ListOfValues  ";
	for (unsigned int i = 0; i < values.size(); i++)
		fOut << values[i] << " ";
	fOut << units << " ;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteParametricAnalysisOnASCII(const std::string name, const std::string type, std::ofstream& fOut, const std::vector<boost::filesystem::path> file_names)
{
	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @Type           " << type << ";" << std::endl;
	fOut << "        @ListOfProfiles " << std::endl;
	for (unsigned int i = 0; i < file_names.size(); i++)
	fOut << "                        " << file_names[i].string() << std::endl;
	fOut << "                         ;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteParametricAnalysisOnASCII(const std::string name, const std::string type, std::ofstream& fOut, const std::vector<double> values1, const std::string units1, const std::vector<double> values2, const std::string units2)
{
	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @Type           " << type << ";" << std::endl;
	fOut << "        @ListOfValues   ";
	for (unsigned int i = 0; i < values1.size(); i++)
		fOut << values1[i] << " ";
	fOut << units1 << " ;" << std::endl;
	fOut << "        @ListOfValues2  ";
	for (unsigned int i = 0; i < values2.size(); i++)
		fOut << values2[i] << " ";
	fOut << units2 << " ;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteIgnitionDelayTimesOnASCII(const std::string name, std::ofstream& fOut)
{
	fOut << "Dictionary " << name << std::endl;
	fOut << "{" << std::endl;
		fOut << "        @Temperature                       true;" << std::endl;
		fOut << "        @Pressure                          true;" << std::endl;
		fOut << "        @RapidCompressionMachine           true;" << std::endl;
		fOut << "        @FilterWidth                       0.1 ms;" << std::endl;
		fOut << "        @RegularizationTimeInterval        2.0 ms;" << std::endl;
		fOut << "        @TemperatureDerivativeThreshold    1.0 K/ms;" << std::endl;
		fOut << "        @Verbose                           true;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

void WriteProfileOnCVS(const boost::filesystem::path file_name,
	const std::string variable1, const double value1, const std::string unit1,
	const std::string variable2, const double value2, const std::string unit2,
	const std::string name1, const std::vector<double> values1, const std::string units1,
	const std::string name2, const std::vector<double> values2, const std::string units2)
{
	std::ofstream fOut(file_name.string(), std::ios::out);
	fOut.setf(std::ios::scientific);
	
	fOut << variable1 << ";" << value1 << " " << unit1 << std::endl;
	fOut << variable2 << ";" << value2 << " " << unit2 << std::endl;
	fOut << name1 << ";" << units1 << std::endl;
	fOut << name2 << ";" << units2 << std::endl;
	fOut << "profile;" << std::endl;
	for (unsigned int i = 0; i < values1.size(); i++)
		fOut << values1[i] << ";" << values2[i] << std::endl;

	fOut.close();
}

void ForceMonotonicProfiles(std::vector<double>& x, std::vector<double>& y)
{
	std::vector<int> indices;
	for (unsigned int i = 1; i < x.size(); i++)
		if (x[i] <= x[i - 1])
			indices.push_back(i - 1);

	std::sort(indices.begin(), indices.end());
	for (std::vector<int>::reverse_iterator i= indices.rbegin(); i != indices.rend(); ++i)
	{
		x.erase(x.begin() + *i);
		y.erase(y.begin() + *i);
	}
}
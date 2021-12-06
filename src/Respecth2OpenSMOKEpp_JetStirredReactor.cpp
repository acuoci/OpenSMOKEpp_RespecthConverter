
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

#include "Respecth2OpenSMOKEpp_JetStirredReactor.h"
#include "Utilities.h"

Respecth2OpenSMOKEpp_JetStirredReactor::Respecth2OpenSMOKEpp_JetStirredReactor
(const boost::filesystem::path file_name,
	const boost::filesystem::path kinetics_folder,
	const boost::filesystem::path output_folder,
	const std::vector<std::string> species_in_kinetic_mech,
	const bool case_sensitive,
	DatabaseSpecies& database_species) :
	Respecth2OpenSMOKEpp(file_name, kinetics_folder, output_folder, species_in_kinetic_mech, case_sensitive, database_species)
{
	// Recognize the apparatus kind
	const std::string apparatus_kind = ptree_.get<std::string>("experiment.apparatus.kind", "unspecified");
	if (apparatus_kind == "stirred reactor")	apparatus_kind_ = ApparatusKind::STIRRED;
	else ErrorMessage("Unknown kind: " + apparatus_kind + ". Available: stirred reactor");
	
	// Read constant values
	std::cout << " * Reading commonProperties section..." << std::endl;
	ReadConstantValueFromXML();

	// Check constant values
	std::cout << " * Checking input data from commonProperties section..." << std::endl;
	if (constant_temperature_ == false && constant_pressure_ == true && constant_composition_ == true &&
		constant_residencetime_ == true && constant_volume_ == true)
		type_ = Type::VARIABLE_T;
	else if (constant_temperature_ == true && constant_pressure_ == false && constant_composition_ == true &&
		constant_residencetime_ == true && constant_volume_ == true)
		type_ = Type::VARIABLE_P;
	else if (constant_temperature_ == true && constant_pressure_ == true && constant_composition_ == true &&
		constant_residencetime_ == true && constant_volume_ == false)
		type_ = Type::VARIABLE_V;
	else if (constant_temperature_ == true && constant_pressure_ == true && constant_composition_ == true &&
		constant_residencetime_ == false && constant_volume_ == true)
		type_ = Type::VARIABLE_TAU;
	else if (constant_temperature_ == true && constant_pressure_ == true && constant_composition_ == false &&
		constant_residencetime_ == true && constant_volume_ == true)
		type_ = Type::VARIABLE_COMPOSITION;
	else if (constant_temperature_ == false && constant_pressure_ == true && constant_composition_ == true &&
		constant_residencetime_ == false && constant_volume_ == true)
		type_ = Type::VARIABLE_TEMPERATURE_TAU;
	else
	{
		ErrorMessage(" Possible combinations of constant variables: (P,X,V,tau) | (T,X,V,tau) | (T,P,X,tau) | (T,P,X,V) | (T,P,V,tau) | (P,X,V)");
	}
	// Read non constant variables
	std::cout << " * Reading dataGroup section..." << std::endl;
	if (constant_temperature_ == false)
	{
		std::cout << "    - temperature..." << std::endl;
		::ReadNonConstantValueFromXML(ptree_, "temperature", t_values_, t_units_);
	}
	if (constant_pressure_ == false)
	{
		std::cout << "    - pressure..." << std::endl;
		::ReadNonConstantValueFromXML(ptree_, "pressure", p_values_, p_units_);
	}
	if (constant_volume_ == false)
	{
		std::cout << "    - volume..." << std::endl;
		::ReadNonConstantValueFromXML(ptree_, "volume", v_values_, v_units_);
	}
	if (constant_residencetime_ == false)
	{
		std::cout << "    - time..." << std::endl;
		::ReadNonConstantValueFromXML(ptree_, "time", tau_values_, tau_units_);
	}
	if (constant_composition_ == false)
	{
		std::cout << "    - composition..." << std::endl;
		::ReadNonConstantValueFromXML(ptree_, "composition", initial_compositions_, species_in_kinetic_mech_, case_sensitive_, database_species_);
	}

	// TODO: Management of variable composition and residence time
	if (type_ == Type::VARIABLE_COMPOSITION) {
		ErrorMessage("Combination of constant variable (P,T,V,tau) not yet implemented!");
	}
	// TODO: Double parametric at variable temperature and residence time
	if (type_ == Type::VARIABLE_COMPOSITION) {
		ErrorMessage("Combination of constant variable (P,X,V) not yet implemented!");
	}
}

void Respecth2OpenSMOKEpp_JetStirredReactor::WriteSimulationData(std::ofstream& fOut)
{
	std::cout << "   - simulation data" << std::endl;

	fOut << "Dictionary PerfectlyStirredReactor" << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @KineticsFolder      " << kinetics_folder_.string() << ";" << std::endl;
	fOut << "        @Type                Isothermal-ConstantPressure;" << std::endl;
	fOut << "        @InletStatus         inlet-status;" << std::endl;
	fOut << "        @ResidenceTime       " << tau_values_[0] << " " << tau_units_ << " ;" << std::endl;
	fOut << "        @Volume              " << v_values_[0] << " " << v_units_ << " ;" << std::endl;
	fOut << "        @EndTime             " << "5" << " " << "s" << " ;" << std::endl;
	fOut << "        @Options             output-options;" << std::endl;
	fOut << "        @ParametricAnalysis  parametric-analysis;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;

	WriteMixStatusOnASCII("inlet-status", fOut, t_values_[0], t_units_, p_values_[0], p_units_, initial_compositions_[0]);

	if (type_ == Type::VARIABLE_T)
		WriteParametricAnalysisOnASCII("parametric-analysis", "temperature", fOut, t_values_, t_units_);
	else if (type_ == Type::VARIABLE_P)
		WriteParametricAnalysisOnASCII("parametric-analysis", "pressure", fOut, p_values_, p_units_);
	else if (type_ == Type::VARIABLE_V)
		WriteParametricAnalysisOnASCII("parametric-analysis", "volume", fOut, v_values_, v_units_);
	else if (type_ == Type::VARIABLE_TAU)
		WriteParametricAnalysisOnASCII("parametric-analysis", "time", fOut, tau_values_, tau_units_);
	
	WriteOutputOptionsOnASCII("output-options", fOut, true, 1000, true, 5000, output_folder_simulation_);
} 

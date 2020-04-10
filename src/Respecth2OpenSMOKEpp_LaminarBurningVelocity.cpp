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

#include "Respecth2OpenSMOKEpp_LaminarBurningVelocity.h"
#include "Utilities.h"

Respecth2OpenSMOKEpp_LaminarBurningVelocity::Respecth2OpenSMOKEpp_LaminarBurningVelocity
(	const boost::filesystem::path file_name,
	const boost::filesystem::path kinetics_folder,
	const boost::filesystem::path output_folder,
	const std::vector<std::string> species_in_kinetic_mech,
	const bool case_sensitive,
	DatabaseSpecies& database_species) :
	Respecth2OpenSMOKEpp(file_name, kinetics_folder, output_folder, species_in_kinetic_mech, case_sensitive, database_species)
{
	// Recognize the apparatus kind
	const std::string apparatus_kind = ptree_.get<std::string>("experiment.apparatus.kind", "unspecified");
	if (apparatus_kind == "flame")	apparatus_kind_ = ApparatusKind::FLAME;
	else ErrorMessage("Unknown kind: " + apparatus_kind + ". Available: flame");

	// Read constant values
	std::cout << " * Reading commonProperties section..." << std::endl;
	ReadConstantValueFromXML();

	// Check constant values
	std::cout << " * Checking input data from commonProperties section..." << std::endl;
	if (constant_temperature_ == true && constant_composition_ == true && constant_pressure_ == false)
		type_ = Type::VARIABLE_P;
	else if (constant_temperature_ == false && constant_composition_ == true && constant_pressure_ == true)
		type_ = Type::VARIABLE_T;
	else if (constant_temperature_ == true && constant_composition_ == false && constant_pressure_ == true)
		type_ = Type::VARIABLE_COMPOSITION;
	else if (constant_temperature_ == false && constant_composition_ == false && constant_pressure_ == true)
		type_ = Type::VARIABLE_T_COMPOSITION;
	else if (constant_temperature_ == true && constant_composition_ == false && constant_pressure_ == false)
		type_ = Type::VARIABLE_P_COMPOSITION;
	else
		ErrorMessage("Possible combinations of constant variables: (P,X) | (T,X) | (T,P) | (P) | (T)");

	// Find pressures
	if (constant_pressure_ == false)
		ReadNonConstantValueFromXML(ptree_, "pressure", p_values_, p_units_);
	
	// Find temperatures
	if (constant_temperature_ == false)
		ReadNonConstantValueFromXML(ptree_, "temperature", t_values_, t_units_);

	// Find list of composition
	if (constant_composition_ == false)
		ReadNonConstantValueFromXML(ptree_, "composition", initial_compositions_, species_in_kinetic_mech_, case_sensitive_, database_species_);

	// Number of simulations
	{
		const unsigned int ns = static_cast<unsigned int>( std::max(p_values_.size(), std::max(t_values_.size(), initial_compositions_.size())) );

		list_inlet_dicts_.resize(ns);
		for (unsigned int i = 0; i < ns; i++)
			list_inlet_dicts_[i] = "inlet-stream-" + std::to_string(i + 1);
	}
}

void Respecth2OpenSMOKEpp_LaminarBurningVelocity::WriteSimulationData(std::ofstream& fOut)
{
	std::cout << "   - simulation data" << std::endl;

	fOut << "Dictionary PremixedLaminarFlame1D" << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @KineticsFolder      " << kinetics_folder_.string() << ";" << std::endl;
	fOut << "        @Type                FlameSpeed;" << std::endl;

	fOut << "        @InletStream         ";
	for (unsigned int i = 0; i < list_inlet_dicts_.size(); i++)
		fOut << list_inlet_dicts_[i] << " ";
	fOut << ";" << std::endl;

	fOut << "        @InletVelocity       50 cm/s;" << std::endl;
	fOut << "        @Grid                grid;" << std::endl;
	fOut << "        @Output              " << output_folder_.string() << ";" << std::endl;
	fOut << "        @UseDaeSolver        true;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;

	if (type_ == Type::VARIABLE_T)
		for (unsigned int i = 0; i < list_inlet_dicts_.size(); i++)
			WriteMixStatusOnASCII(list_inlet_dicts_[i], fOut, t_values_[i], t_units_, p_values_[0], p_units_, initial_compositions_[0]);

	if (type_ == Type::VARIABLE_P)
		for (unsigned int i = 0; i < list_inlet_dicts_.size(); i++)
			WriteMixStatusOnASCII(list_inlet_dicts_[i], fOut, t_values_[0], t_units_, p_values_[i], p_units_, initial_compositions_[0]);

	if (type_ == Type::VARIABLE_COMPOSITION)
		for (unsigned int i = 0; i < list_inlet_dicts_.size(); i++)
			WriteMixStatusOnASCII(list_inlet_dicts_[i], fOut, t_values_[0], t_units_, p_values_[0], p_units_, initial_compositions_[i]);

	if (type_ == Type::VARIABLE_T_COMPOSITION)
		for (unsigned int i = 0; i < list_inlet_dicts_.size(); i++)
			WriteMixStatusOnASCII(list_inlet_dicts_[i], fOut, t_values_[i], t_units_, p_values_[0], p_units_, initial_compositions_[i]);

	if (type_ == Type::VARIABLE_P_COMPOSITION)
		for (unsigned int i = 0; i < list_inlet_dicts_.size(); i++)
			WriteMixStatusOnASCII(list_inlet_dicts_[i], fOut, t_values_[0], t_units_, p_values_[i], p_units_, initial_compositions_[i]);

	fOut << "Dictionary grid" << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @Length                5 cm;" << std::endl;
	fOut << "        @InitialPoints         12;" << std::endl;
	fOut << "        @Type                  database;" << std::endl;
	fOut << "        @MaxPoints             400;" << std::endl;
	fOut << "        @MaxAdaptivePoints     15;" << std::endl;
	fOut << "        @GradientCoefficient   0.05;" << std::endl;
	fOut << "        @CurvatureCoefficient  0.5;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;
}

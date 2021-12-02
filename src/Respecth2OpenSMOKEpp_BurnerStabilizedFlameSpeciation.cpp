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

#include "Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation.h"
#include "Utilities.h"

Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation::Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation
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
	if (apparatus_kind == "flame")	apparatus_kind_ = ApparatusKind::FLAME;
	else 
		ErrorMessage("Unknown kind: " + apparatus_kind + ". Available: flame");
	// Recognize the operation mode
	const std::string apparatus_mode = ptree_.get<std::string>("experiment.apparatus.kind.mode", "burner-stabilized");
	if (apparatus_mode == "burner-stabilized")	apparatus_mode_ = ApparatusMode::BURNER_STABILIZED;
	else 
		ErrorMessage("Unknown mode: " + apparatus_mode + ". Available: burner-stabilized");

	// Read constant values
	std::cout << " * Reading commonProperties section..." << std::endl;
	ReadConstantValueFromXML();

	// Check constant values
	std::cout << " * Checking input data from commonProperties section..." << std::endl;
	if (constant_temperature_ == true && constant_pressure_ == true && constant_composition_ == true && constant_massflowrate_ == true && constant_laminarburningvelocity_ == false)
		type_ = Type::ASSIGNED_M;
	else if (constant_temperature_ == true && constant_pressure_ == true && constant_composition_ == true && constant_massflowrate_ == false && constant_laminarburningvelocity_ == true)
		type_ = Type::ASSIGNED_SL;
	else
		ErrorMessage("(T,P,X,m) or (T,P,X,sl) must be defined as constant variables");

	// Read space-temperature profile
	std::cout << " * Reading dataGroup section (distance)..." << std::endl;
	ReadNonConstantValueFromXML(ptree_, "distance", x_profile_values_, x_profile_units_);
	std::cout << " * Reading dataGroup section (temperature)..." << std::endl;
	ReadNonConstantValueFromXML(ptree_, "temperature", t_profile_values_, t_profile_units_);

	// Recognize if the temperature profile is fixed or not
	fixed_temperature_profile_ = true;
	if (t_profile_values_.size() == 0)
		fixed_temperature_profile_ = false;

	// Check for first point
	if (fixed_temperature_profile_ == true)
	{
		if (x_profile_values_[0] != 0.)
		{
			x_profile_values_.insert(x_profile_values_.begin(), 0.);
			t_profile_values_.insert(t_profile_values_.begin(), t_values_[0]);
		}
	}
}

void Respecth2OpenSMOKEpp_BurnerStabilizedFlameSpeciation::WriteSimulationData(std::ofstream& fOut)
{
	std::cout << "   - simulation data" << std::endl;

	fOut << "Dictionary PremixedLaminarFlame1D" << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @KineticsFolder          " << kinetics_folder_.string() << ";" << std::endl;
	fOut << "        @Type                    BurnerStabilized;" << std::endl;
	fOut << "        @InletStream             inlet-stream;" << std::endl;
	
	if (type_ == Type::ASSIGNED_M)
		fOut << "        @InletMassFlux           " << m_values_[0] << " " << m_units_ << ";" << std::endl;
	else if (type_ == Type::ASSIGNED_SL)
		fOut << "        @InletVelocity           " << sl_values_[0] << " " << sl_units_ << ";" << std::endl;

	fOut << "        @Grid                    grid;" << std::endl;
	fOut << "        @Output                  " << output_folder_simulation_.string() << ";" << std::endl;
	fOut << "        @UseDaeSolver            true;" << std::endl;
	if (fixed_temperature_profile_ == true)
		fOut << "        @FixedTemperatureProfile T-Profile;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;

	WriteMixStatusOnASCII("inlet-stream", fOut, t_values_[0], t_units_, p_values_[0], p_units_, initial_compositions_[0]);

	fOut << "Dictionary grid" << std::endl;
	fOut << "{" << std::endl;
	if (fixed_temperature_profile_ == true)
		fOut << "        @Length                " << x_profile_values_.back() << " " << x_profile_units_ << " ;" << std::endl;
	else
		fOut << "        @Length                10 cm;" << std::endl;
	fOut << "        @InitialPoints			12;" << std::endl;
	fOut << "        @Type					database;" << std::endl;
	fOut << "        @MaxPoints				400;" << std::endl;
	fOut << "        @MaxAdaptivePoints		15;" << std::endl;
	fOut << "        @GradientCoefficient	0.05;" << std::endl;
	fOut << "        @CurvatureCoefficient	0.5;" << std::endl;
	fOut << "}" << std::endl;
	fOut << std::endl;

	if (fixed_temperature_profile_ == true)
	{
		fOut << "Dictionary T-Profile" << std::endl;
		fOut << "{" << std::endl;
		fOut << "        @XVariable length;" << std::endl;
		fOut << "        @YVariable temperature;" << std::endl;
		fOut << "        @XUnits    " << x_profile_units_ << " ;" << std::endl;
		fOut << "        @YUnits    " << t_profile_units_ << " ;" << std::endl;
		fOut << "        @Profile" << std::endl;
		for (unsigned int j = 0; j < t_profile_values_.size(); j++)
			fOut << "        " << x_profile_values_[j] << " " << t_profile_values_[j] << std::endl;
		fOut << "        ;" << std::endl;
		fOut << "}" << std::endl;
	}

	fOut << std::endl;
}
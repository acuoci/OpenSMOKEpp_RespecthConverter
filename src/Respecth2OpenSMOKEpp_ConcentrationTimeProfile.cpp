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

#include "Respecth2OpenSMOKEpp_ConcentrationTimeProfile.h"
#include "Utilities.h"

Respecth2OpenSMOKEpp_ConcentrationTimeProfile::Respecth2OpenSMOKEpp_ConcentrationTimeProfile
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
	if (apparatus_kind == "flow reactor")		apparatus_kind_ = ApparatusKind::FLOW_REACTOR;
	else if (apparatus_kind == "shock tube")	apparatus_kind_ = ApparatusKind::SHOCK_TUBE;
	else if (apparatus_kind == "batch")			apparatus_kind_ = ApparatusKind::BATCH;
	else ErrorMessage("Unknown kind: " + apparatus_kind + ". Available: flow reactor | shock tube | batch");

	// Read constant values
	std::cout << " * Reading commonProperties section..." << std::endl;
	ReadConstantValueFromXML();

	// Check constant values
	std::cout << " * Checking input data from commonProperties section..." << std::endl;
	if (constant_temperature_ == false || constant_pressure_ == false || constant_composition_ == false)
		ErrorMessage(" T,P and X must be defined as constant variables");

	// Read time profile
	ReadNonConstantValueFromXML(ptree_, "time", time_profile_values_, time_profile_units_);
}

void Respecth2OpenSMOKEpp_ConcentrationTimeProfile::WriteSimulationData(std::ofstream& fOut)
{
	std::cout << "   - simulation data" << std::endl;

	if (apparatus_kind_ == ApparatusKind::FLOW_REACTOR)
	{
		fOut << "Dictionary PlugFlowReactor" << std::endl;
		fOut << "{" << std::endl;
		fOut << "        @KineticsFolder          " << kinetics_folder_.string() << ";" << std::endl;
		fOut << "        @Type                    Isothermal;" << std::endl;
		fOut << "        @InletStatus             mix-status;" << std::endl;
		fOut << "        @ResidenceTime           " << time_profile_values_.back() << " " << time_profile_units_ << ";" << std::endl;
		fOut << "        @ConstantPressure        true;" << std::endl;
		fOut << "        @Velocity                10 cm/s;" << std::endl;
		fOut << "        @Options                 output-options;" << std::endl;
		fOut << "}" << std::endl;
		fOut << std::endl;
	}
	else if (apparatus_kind_ == ApparatusKind::SHOCK_TUBE)
	{
		fOut << "Dictionary ShockTubeReactor" << std::endl;
		fOut << "{" << std::endl;
		fOut << "        @KineticsFolder          " << kinetics_folder_.string() << ";" << std::endl;
		fOut << "        @Type                    ReflectedShock;" << std::endl;
		fOut << "        @ReflectedShockStatus    mix-status;" << std::endl;
		fOut << "        @EndTime                 " << time_profile_values_.back() << " " << time_profile_units_ << ";" << std::endl;
		fOut << "        @Options                 output-options;" << std::endl;
		fOut << "}" << std::endl;
		fOut << std::endl;
	}
	else if (apparatus_kind_ == ApparatusKind::BATCH)
	{
		fOut << "Dictionary BatchReactor" << std::endl;
		fOut << "{" << std::endl;
		fOut << "        @KineticsFolder          " << kinetics_folder_.string() << std::endl;
		fOut << "        @Type                    Isothermal-ConstantPressure;" << std::endl;
		fOut << "        @IninitialtStatus        mix-status;" << std::endl;
		fOut << "        @EndTime                 " << time_profile_values_.back() << " " << time_profile_units_ << ";" << std::endl;
		fOut << "        @Options                 output-options;" << std::endl;
		fOut << "}" << std::endl;
		fOut << std::endl;
	}

	WriteMixStatusOnASCII("mix-status", fOut, t_values_[0], t_units_, p_values_[0], p_units_, initial_compositions_[0]);

	WriteOutputOptionsOnASCII("output-options", fOut, 1, output_folder_);
}

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

#include "Respecth2OpenSMOKEpp_IgnitionDelay.h"
#include "Utilities.h"

Respecth2OpenSMOKEpp_IgnitionDelay::Respecth2OpenSMOKEpp_IgnitionDelay
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
	if (apparatus_kind == "flow reactor")					apparatus_kind_ = ApparatusKind::FLOW_REACTOR;
	else if (apparatus_kind == "shock tube")				apparatus_kind_ = ApparatusKind::SHOCK_TUBE;
	else if (apparatus_kind == "rapid compression machine")	apparatus_kind_ = ApparatusKind::RCM;
	else ErrorMessage("Unknown kind: " + apparatus_kind + ". Available: flow reactor | shock tube | rapid compression machine");

	// Recognize the ignition type
	std::cout << " * Reading ignition delay time type section..." << std::endl;
	ReadIdtTypeFromXML();

	// Read constant values
	std::cout << " * Reading commonProperties section..." << std::endl;
	ReadConstantValueFromXML();

	// Check constant values 
	std::cout << " * Checking input data from commonProperties section..." << std::endl;
	if (constant_temperature_ == false && constant_pressure_ == true && constant_composition_ == true)
		type_ = Type::VARIABLE_T;
	if (constant_temperature_ == true && constant_pressure_ == false && constant_composition_ == true)
		type_ = Type::VARIABLE_P;
	if (constant_temperature_ == false && constant_pressure_ == false && constant_composition_ == true)
		type_ = Type::VARIABLE_TP;
	if (constant_composition_ == false)
		ErrorMessage("Only constant composition is allowed");
	if (constant_temperature_ == true && constant_pressure_ == true)
		ErrorMessage("Pressure and Temperature cannot be constant at the same time");

	// Read temperatures
	if (constant_temperature_ == false)
	{
		std::cout << " * Reading dataGroup section (temperature)..." << std::endl;
		ReadNonConstantValueFromXML(ptree_, "temperature", t_values_, t_units_);
	}

	// Read pressures
	if (constant_pressure_ == false)
	{
		std::cout << " * Reading dataGroup section (pressure)..." << std::endl;
		ReadNonConstantValueFromXML(ptree_, "pressure", p_values_, p_units_);
	}

	// Read ignition delay times
	std::cout << " * Reading dataGroup section (ignition delay)..." << std::endl;
	ReadNonConstantValueFromXML(ptree_, "ignition delay", tau_values_, tau_units_);

	// Check for possible v-t history
	ReadProfileFromXML(ptree_, "V-t history", "volume", v_history_values_, v_history_units_, "time", tau_history_values_, tau_history_units_ );

	// Check for monoticity of profiles (if any)
	for (unsigned int i = 0; i < v_history_values_.size(); i++)
		ForceMonotonicProfiles(tau_history_values_[i], v_history_values_[i]);

	// Select a suitable maximum time for integration
	tau_max_ = *std::max_element(std::begin(tau_values_), std::end(tau_values_)) * 2;
}

void Respecth2OpenSMOKEpp_IgnitionDelay::WriteSimulationData(std::ofstream& fOut)
{
	std::cout << "   - simulation data" << std::endl;

	fOut << "Dictionary BatchReactor" << std::endl;
	fOut << "{" << std::endl;
	fOut << "        @KineticsFolder          " << kinetics_folder_.string() << ";" << std::endl;
	if (v_history_values_.size() != 0 ) // per ora non consideriamo PrssureCoefficient andrebbe messo così|| dpdt_values_.size() != 0
		fOut << "        @Type                    NonIsothermal-UserDefinedVolume;" << std::endl;
	else
		fOut << "        @Type                    NonIsothermal-ConstantVolume;" << std::endl;
	fOut << "        @InitialStatus           mix-status;" << std::endl;
	fOut << "        @EndTime                 " << tau_max_ << " " << tau_units_ << ";" << std::endl;
	fOut << "        @Volume                  1 cm3;" << std::endl;
	fOut << "        @OdeParameters           ode-parameters;" << std::endl;
	fOut << "        @Options                 output-options;" << std::endl;
	fOut << "        @ParametricAnalysis      parametric-analysis;" << std::endl;
	fOut << "        @IgnitionDelayTimes      ignition-delay-times;" << std::endl;
	
	if (dpdt_values_.size() != 0)
		fOut << "        //@PressureCoefficient     " << dpdt_values_[0] << dpdt_units_ << ";" << std::endl;

	fOut << "}" << std::endl;
	fOut << std::endl;

	WriteMixStatusOnASCII("mix-status", fOut, t_values_[0], t_units_, p_values_[0], p_units_, initial_compositions_[0]);

	WriteODEParametersOnASCII("ode-parameters", fOut, 1e-14, 1e-7);

	if (apparatus_kind_ == ApparatusKind::RCM)
		WriteIgnitionDelayTimesOnASCII("ignition-delay-times", fOut, true, idt_);
	else
		WriteIgnitionDelayTimesOnASCII("ignition-delay-times", fOut, false, idt_);

	if (v_history_units_.size() == 0)
	{
		if (type_ == Type::VARIABLE_TP)
			WriteParametricAnalysisOnASCII("parametric-analysis", "temperature-pressure", fOut, t_values_, t_units_, p_values_, p_units_);

		if (type_ == Type::VARIABLE_T)
			WriteParametricAnalysisOnASCII("parametric-analysis", "temperature", fOut, t_values_, t_units_);

		if (type_ == Type::VARIABLE_P)
			WriteParametricAnalysisOnASCII("parametric-analysis", "pressure", fOut, p_values_, p_units_);
	}
	else
	{
		std::vector<boost::filesystem::path> list_of_files(v_history_units_.size());
		for (unsigned int i = 0; i < list_of_files.size(); i++)
		{
			list_of_files[i] = file_name_xml_.stem();
			list_of_files[i] += ".";  list_of_files[i] += std::to_string(i + 1); list_of_files[i] += ".csv";
		}
		WriteParametricAnalysisOnASCII("parametric-analysis", "temperature-pressure", fOut, list_of_files);
	}

	WriteOutputOptionsOnASCII("output-options", fOut, true, 1000, true, 5, output_folder_simulation_);
}

void Respecth2OpenSMOKEpp_IgnitionDelay::WriteAdditionalFiles()
{
	if (v_history_units_.size() != 0)
	{
		std::cout << "   - additional files" << std::endl;

		for (unsigned int i = 0; i < v_history_units_.size(); i++)
		{
			boost::filesystem::path file_name_csv = output_folder_ / file_name_xml_.stem();
			file_name_csv += ".";  file_name_csv += std::to_string(i + 1); file_name_csv += ".csv";
			
			if (type_ == Type::VARIABLE_TP)
				WriteProfileOnCSV(file_name_csv,
					"temperature", t_values_[i], t_units_, "pressure", p_values_[i], p_units_,
					"time", tau_history_values_[i], tau_history_units_[i], "volume", v_history_values_[i], v_history_units_[i]);
			else if(type_ == Type::VARIABLE_T)
				WriteProfileOnCSV(file_name_csv,
					"temperature", t_values_[i], t_units_, "pressure", p_values_[0], p_units_,
					"time", tau_history_values_[i], tau_history_units_[i], "volume", v_history_values_[i], v_history_units_[i]);
			else if (type_ == Type::VARIABLE_P)
				WriteProfileOnCSV(file_name_csv,
					"temperature", t_values_[0], t_units_, "pressure", p_values_[i], p_units_,
					"time", tau_history_values_[i], tau_history_units_[i], "volume", v_history_values_[i], v_history_units_[i]);
		}
	}	
}

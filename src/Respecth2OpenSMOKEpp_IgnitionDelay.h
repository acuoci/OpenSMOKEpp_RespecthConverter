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

#pragma once

#include "Respecth2OpenSMOKEpp.h"

class Respecth2OpenSMOKEpp_IgnitionDelay : public Respecth2OpenSMOKEpp
{
public:

	Respecth2OpenSMOKEpp_IgnitionDelay(	const boost::filesystem::path file_name,
										const boost::filesystem::path kinetics_folder,
										const boost::filesystem::path output_folder,
										const std::vector<std::string> species_in_kinetic_mech,
										const bool case_sensitive,
										DatabaseSpecies& database_species);

private:

	enum class ApparatusKind { FLOW_REACTOR, SHOCK_TUBE, RCM }	apparatus_kind_;
	enum class Type { VARIABLE_T, VARIABLE_P, VARIABLE_TP }		type_;

	double tau_max_;

	std::vector<std::vector<double>> v_history_values_;
	std::vector<std::string> v_history_units_;

	std::vector<std::vector<double>> tau_history_values_;
	std::vector<std::string> tau_history_units_;

	virtual void WriteSimulationData(std::ofstream& fOut);

	virtual void WriteAdditionalFiles();

};
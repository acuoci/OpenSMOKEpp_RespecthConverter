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

#include "Conversions.h"

void CheckAndConvertUnits(const std::string name, double& value, std::string& units)
{
	if (name == "equivalence ratio")
	{
		// do nothing
	}
	else if (name == "temperature")
	{
		if (units != "K")
			ConversionErrorMessage("Unknown temperature units: " + units + ". Available units: K");
	}
	else if (name == "pressure")
	{
		if (units != "atm" && units != "bar" && units != "torr" && units != "Torr" && units != "Pa" && units != "kPa" && units != "MPa" && units != "mbar")
			ConversionErrorMessage("Unknown pressure units: " + units + ". Available units: atm | bar | mbar | torr | Torr | Pa | kPa | MPa");

		if (units == "torr" || units == "Torr")
		{
			value /= 760.;
			units = "atm";
		}
		else if (units == "kPa")
		{
			value *= 1.e3;
			units = "Pa";
		}
		else if (units == "MPa")
		{
			value *= 1.e6;
			units = "Pa";
		}
		else if (units == "mbar")
		{
			value /= 1000.;
			units = "bar";
		}
	}
	else if (name == "residence time" || name == "ignition delay" || name == "time")
	{
		if (units != "s" && units != "ms" && units != "us" && units != "ns" && units != "min")
			ConversionErrorMessage("Unknown time units: " + units + ". Available units: s | ms | us | ns | min");

		if (units == "ns")
		{
			value /= 1.e6;
			units = "ms";
		}
		else if (units == "us")
		{
			value /= 1.e3;
			units = "ms";
		}
	}
	else if (name == "volume")
	{
		if (units != "m3" && units != "cm3" && units != "mm3" && units != "dm3" && units != "L")
			ConversionErrorMessage("Unknown volume units: " + units + ". Available units: m3 | dm3 | cm3 | mm3 | L");

		if (units == "L")
		{
			units = "dm3";
		}
	}
	else if (name == "flow rate")
	{
		if (units != "g cm-2 s-1" && units != "kg m-2 s-1" && units != "us" && units != "ns" && units != "min")
			ConversionErrorMessage("Unknown flow rate units: " + units + ". Available units: g cm-2 s-1 | kg m-2 s-1");

		if (units == "g cm-2 s-1")
		{
			units = "g/cm2/s";
		}
		else if (units == "kg m-2 s-1")
		{
			units = "kg/m2/s";
		}
	}
	else if (name == "pressure rise")
	{
		if (units != "ms-1" && units != "s-1")
			ConversionErrorMessage("Unknown pressure rise units: " + units + ". Available units: ms-1 | s-1");

		if (units == "ms-1")
		{
			units = "1/ms";
		}
		else if (units == "s-1")
		{
			units = "1/s";
		}
	}
	else if (name == "distance")
	{
		if (units != "m" && units != "dm" && units != "cm" && units != "mm")
			ConversionErrorMessage("Unknown distance units: " + units + ". Available units: m | dm | cm | mm");
	}
	else if (name == "laminar burning velocity")
	{
		if (units != "m/s" && units != "cm/s" && units != "mm/s")
			ConversionErrorMessage("Unknown laminar burning velocity units: " + units + ". Available units: m/s | cm/s | mm/s");
	}
	else
	{
		ConversionErrorMessage("Unknown variable: " + name);
	}
}

void CheckAndConvertUnits(const std::string name, std::vector<double>& values, std::string& units)
{
	for (unsigned int i = 0; i < values.size(); i++)
		CheckAndConvertUnits(name, values[i], units);
}

void ConversionErrorMessage(const std::string message)
{
	std::cout << "Error in conversion of units: " << message << std::endl;
	std::cout << "Press enter to exit... ";
	getchar();
	exit(-1);
}

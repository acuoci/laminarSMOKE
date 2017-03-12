/*----------------------------------------------------------------------*\
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
|   Copyright(C) 2014, 2013, 2012  Mattia Bissoli and Alberto Cuoci       |
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

#ifndef OpenSMOKE_CONVERSION_H
#define	OpenSMOKE_CONVERSION_H

#include "math/PhysicalConstants.h"

namespace Conversions
{
	double Length(const double value, const std::string units);
	double U_Length(const double value, const std::string units);
	double Area(const double value, const std::string units);
	double Volume(const double value, const std::string units);
	double Specific_Volume(const double value, const std::string units);
	double Pressure(const double value, const std::string units);
	double Time(const double value, const std::string units);
	double Energy(const double value, const std::string units);
	double Entropy(const double value, const std::string units);
	double Specific_Energy(const double value, const std::string units);
	double Specific_Entropy(const double value, const std::string units);
	double Mass(const double value, const std::string units);
	double Temperature(const double value, const std::string units);
	double Frequency(const double value, const std::string units);
	double Velocity(const double value, const std::string units);
	double Mass_Flow_Rate(const double value, const std::string units);
	double Mole_Flow_Rate(const double value, const std::string units);
	double Volumetric_Flow_Rate(const double value, const std::string units);
	double Heat_Flux(const double value, const std::string units);
	double Heat_Exchange_Coefficient(const double value, const std::string units);
	double Dynamic_Viscosity(const double value, const std::string units);
	double Density(const double value, const std::string units);
	double Angle(const double value, const std::string units);
	double Angular_Velocity(const double value, const std::string units);
	double Area_Velocity(const double value, const std::string units);
	double Specific_Energy_Molar(const double value, const std::string units);

	const double Pa_from_atm = 101325.0;
	const double Pa_from_bar = 100000.0;
	const double Pa_from_mbar = 100.0;
	const double Pa_from_torr = 133.32237;
	const double Pa_from_kPa = 1000.;
	const double Pa_from_psi = 6894.75728;

	const double m_from_mm = 0.0010;
	const double m_from_cm = 0.0100;
	const double m_from_in = 0.0254;
	const double m_from_ft = 0.3048;
	const double m_from_km = 0.3048;

	const double m2_from_mm2 = 1.e-6;
	const double m2_from_cm2 = 1.e-4;
	const double m2_from_in2 = 6.4516e-4;
	const double m2_from_ft2 = 0.09290304;

	const double m3_from_mm3 = 1.e-9;
	const double m3_from_cm3 = 1.e-6;
	const double m3_from_l = 1.e-3;
	const double m3_from_in3 = 1.6387064e-5;
	const double m3_from_ft3 = 0.028316846592;
	const double m3_from_gallon_UK = 0.00454609;
	const double m3_from_oz_UK = 0.0000284130625;
	const double m3_from_oz_USA = 0.000029573529563;
	const double m3_from_gallon_dry_USA = 0.0044048838;
	const double m3_from_gallon_liq_USA = 0.003785411784;

	const double kg_from_g = 1.e-3;
	const double kg_from_lb = 0.454;
	const double kg_from_oz = 0.028349;

	const double s_from_hr = 3600.;
	const double s_from_min = 60.00;
	const double s_from_ms = 1e-3;

	const double J_from_kJ = 1.e3;
	const double J_from_erg = 1.e-7;
	const double J_from_BTU = 1055.0559;
	const double J_from_kWh = 3600000.;
	const double J_from_cal = PhysicalConstants::R_J_mol / PhysicalConstants::R_cal_mol;
	const double J_from_kcal = 1.e3 * J_from_cal;
	const double J_from_eV = 1.6021765314e-19;

	const double kmol_from_mol = 1.e-3;

	const double W_from_kW = 1.e3;

	const double rad_to_degrees = 360. / (2. * PhysicalConstants::pi);
	const double degrees_to_rad = (2. * PhysicalConstants::pi) / 360.;
}

#include "math/Conversions.hpp"

#endif	/* OpenSMOKE_CONVERSION_H */

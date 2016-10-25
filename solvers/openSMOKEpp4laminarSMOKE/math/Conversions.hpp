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

namespace Conversions
{

double Length(const double value, const std::string units)
{
    if (units == "cm") return value * m_from_cm;
    else if (units == "m") return value;
    else if (units == "mm") return value * m_from_mm;
    else if (units == "in") return value * m_from_in;
    else if (units == "ft") return value * m_from_ft;

    else return OpenSMOKE::FatalErrorMessage("Wrong length units: " + units);
}

double U_Length(const double value, const std::string units)
{
    if (units == "1/cm") return value / m_from_cm;
    else if (units == "1/m") return value;
    else if (units == "1/mm") return value / m_from_mm;
    else if (units == "1/in") return value / m_from_in;
    else if (units == "1/ft") return value / m_from_ft;

    else return OpenSMOKE::FatalErrorMessage("Wrong u_length units: " + units);
}

double Area(const double value, const std::string units)
{
    if (units == "cm2") return value * m2_from_cm2;
    else if (units == "m2") return value;
    else if (units == "mm2") return value * m2_from_mm2;
    else if (units == "in2") return value * m2_from_in2;
    else if (units == "ft2") return value * m2_from_ft2;

    else return OpenSMOKE::FatalErrorMessage("Wrong area units: " + units);
}

double Volume(const double value, const std::string units)
{
    if (units == "mm3") return value * m3_from_mm3;
    else if (units == "m3") return value;
    else if (units == "cm3") return value * m3_from_cm3;
    else if (units == "l") return value * m3_from_l;
    else if (units == "in3") return value * m3_from_in3;
    else if (units == "ft3") return value * m3_from_ft3;
    else if (units == "gallon_UK") return value * m3_from_gallon_UK;
    else if (units == "oz_UK") return value * m3_from_oz_UK;
    else if (units == "oz_USA") return value * m3_from_oz_USA;
    else if (units == "gallon_dry_USA") return value * m3_from_gallon_dry_USA;
    else if (units == "gallon_liq_USA") return value * m3_from_gallon_liq_USA;

    else return OpenSMOKE::FatalErrorMessage("Wrong volume units: " + units);
}

double Specific_Volume(const double value, const std::string units)
{

    if (units == "m3/kg") return value;
    else if (units == "cm3/g") return value / kg_from_g * m3_from_cm3;
    else if (units == "ft3/lb") return value / kg_from_lb * m3_from_ft3;
    else if (units == "ft3/oz") return value / kg_from_oz * m3_from_ft3;

    else return OpenSMOKE::FatalErrorMessage("Wrong specific volume units: " + units);
}

double Pressure(const double value, const std::string units)
{
    if (units == "atm") return value * Pa_from_atm;
    else if (units == "Pa") return value;
    else if (units == "bar") return value * Pa_from_bar;
    else if (units == "mbar") return value * Pa_from_mbar;
    else if (units == "torr") return value * Pa_from_torr;
    else if (units == "kPa") return value * Pa_from_kPa;
    else if (units == "psi") return value * Pa_from_psi;

    else return OpenSMOKE::FatalErrorMessage("Wrong pressure units: " + units);
}

double Mass(const double value, const std::string units)
{
    if (units == "g") return value * kg_from_g;
    else if (units == "kg") return value;
    else if (units == "lb") return value * kg_from_lb;
    else if (units == "oz") return value * kg_from_oz;

    else return OpenSMOKE::FatalErrorMessage("Wrong mass units: " + units);
}

double Time(const double value, const std::string units)
{
    if (units == "s") return value;
    else if (units == "min") return value * s_from_min;
    else if (units == "ms") return value * s_from_ms;
    else if (units == "hr") return value * s_from_hr;

    else return OpenSMOKE::FatalErrorMessage("Wrong time units: " + units);
}

double Energy(const double value, const std::string units)
{
    if (units == "J") return value;
    else if (units == "kJ") return value * J_from_kJ;
    else if (units == "cal") return value * J_from_cal;
    else if (units == "kcal") return value * J_from_kcal;
    else if (units == "kWh") return value * J_from_kWh;
    else if (units == "BTU") return value * J_from_BTU;
    else if (units == "erg") return value * J_from_erg;
    else if (units == "eV") return value * J_from_eV;

    else return OpenSMOKE::FatalErrorMessage("Wrong energy units: " + units);
}

double Entropy(const double value, const std::string units)
{
    if (units == "J/K") return value;
    else if (units == "kJ/K") return value * J_from_kJ;
    else if (units == "cal/K") return value * J_from_cal;
    else if (units == "kcal/K") return value * J_from_kcal;
    else if (units == "kWh/K") return value * J_from_kWh;
    else if (units == "BTU/K") return value * J_from_BTU;
    else if (units == "erg/K") return value * J_from_erg;
    else if (units == "eV/K") return value * J_from_eV;

    else return OpenSMOKE::FatalErrorMessage("Wrong entropy units: " + units);
}

double Specific_Energy(const double value, const std::string units)
{
    if (units == "J/kg") return value;
    else if (units == "kJ/kg") return value * J_from_kJ;
    else if (units == "cal/kg") return value * J_from_cal;
    else if (units == "kcal/kg") return value * J_from_kcal;
    else if (units == "kWh/kg") return value * J_from_kWh;
    else if (units == "BTU/kg") return value * J_from_BTU;
    else if (units == "erg/kg") return value * J_from_erg;
    else if (units == "eV/kg") return value * J_from_eV;
    else if (units == "J/g") return value / kg_from_g;
    else if (units == "kJ/g") return value * J_from_kJ / kg_from_g;
    else if (units == "cal/g") return value * J_from_cal / kg_from_g;
    else if (units == "kcal/g") return value * J_from_kcal / kg_from_g;
    else if (units == "kWh/g") return value * J_from_kWh / kg_from_g;
    else if (units == "BTU/g") return value * J_from_BTU / kg_from_g;
    else if (units == "erg/g") return value * J_from_erg / kg_from_g;
    else if (units == "eV/g") return value * J_from_eV / kg_from_g;

    else return OpenSMOKE::FatalErrorMessage("Wrong specific energy units: " + units);
}

double Specific_Energy_Molar(const double value, const std::string units)
{
    if (units == "J/kmol") return value;
    else if (units == "kJ/kmol") return value * J_from_kJ;
    else if (units == "cal/kmol") return value * J_from_cal;
    else if (units == "kcal/kmol") return value * J_from_kcal;
    else if (units == "kWh/kmol") return value * J_from_kWh;
    else if (units == "BTU/kmol") return value * J_from_BTU;
    else if (units == "erg/kmol") return value * J_from_erg;
    else if (units == "eV/kmol") return value * J_from_eV;
    else if (units == "J/mol") return value / kmol_from_mol;
    else if (units == "kJ/mol") return value * J_from_kJ / kmol_from_mol;
    else if (units == "cal/mol") return value * J_from_cal / kmol_from_mol;
    else if (units == "kcal/mol") return value * J_from_kcal / kmol_from_mol;
    else if (units == "kWh/mol") return value * J_from_kWh / kmol_from_mol;
    else if (units == "BTU/mol") return value * J_from_BTU / kmol_from_mol;
    else if (units == "erg/mol") return value * J_from_erg / kmol_from_mol;
    else if (units == "eV/mol") return value * J_from_eV / kmol_from_mol;

    else return OpenSMOKE::FatalErrorMessage("Wrong specific energy units: " + units);
}

double Specific_Entropy(const double value, const std::string units)
{
    if (units == "J/kg/K") return value;
    else if (units == "kJ/kg/K") return value * J_from_kJ;
    else if (units == "cal/kg/K") return value * J_from_cal;
    else if (units == "kcal/kg/K") return value * J_from_kcal;
    else if (units == "kWh/kg/K") return value * J_from_kWh;
    else if (units == "BTU/kg/K") return value * J_from_BTU;
    else if (units == "erg/kg/K") return value * J_from_erg;
    else if (units == "eV/kg/K") return value * J_from_eV;
    else if (units == "J/g/K") return value / kg_from_g;
    else if (units == "kJ/g/K") return value * J_from_kJ / kg_from_g;
    else if (units == "cal/g/K") return value * J_from_cal / kg_from_g;
    else if (units == "kcal/g/K") return value * J_from_kcal / kg_from_g;
    else if (units == "kWh/g/K") return value * J_from_kWh / kg_from_g;
    else if (units == "BTU/g/K") return value * J_from_BTU / kg_from_g;
    else if (units == "erg/g/K") return value * J_from_erg / kg_from_g;
    else if (units == "eV/g/K") return value * J_from_eV / kg_from_g;

    else return OpenSMOKE::FatalErrorMessage("Wrong specific entropy units: " + units);
}

double Temperature(const double value, const std::string units)
{
    if (units == "K") return value;
    else if (units == "C") return value + 273.15;
    else if (units == "F") return (value + 459.67)*5. / 9.;
    else if (units == "R") return value / 1.80;

    else return OpenSMOKE::FatalErrorMessage("Wrong temperature units: " + units);
}

double Frequency(const double value, const std::string units)
{
    if (units == "Hz") return value;
    else if (units == "1/s") return value;
    else if (units == "1/min") return value / s_from_min;
    else if (units == "1/ms") return value / s_from_ms;
    else if (units == "1/hr") return value / s_from_hr;

    else return OpenSMOKE::FatalErrorMessage("Wrong frequency units: " + units);
}

double Velocity(const double value, const std::string units)
{
    if (units == "m/s") return value;
    else if (units == "cm/s") return value * m_from_cm;
    else if (units == "mm/s") return value * m_from_mm;
    else if (units == "km/hr") return value * m_from_km / s_from_hr;
    else if (units == "cm/min") return value * m_from_cm / s_from_min;
    else if (units == "m/min") return value / s_from_min;
    else if (units == "in/s") return value * m_from_in;
    else if (units == "in/min") return value * m_from_in / s_from_min;
    else if (units == "ft/s") return value * m_from_ft;
    else if (units == "ft/min") return value * m_from_ft / s_from_min;

    else return OpenSMOKE::FatalErrorMessage("Wrong velocity units: " + units);
}

double Area_Velocity(const double value, const std::string units)
{
    if (units == "m2/s") return value;
    else if (units == "cm2/s") return value * m_from_cm * m_from_cm;
    else if (units == "mm2/s") return value * m_from_mm * m_from_mm;
    else if (units == "km2/hr") return value * m_from_km * m_from_km / s_from_hr;
    else if (units == "cm2/min") return value * m_from_cm * m_from_cm / s_from_min;
    else if (units == "m2/min") return value / s_from_min;
    else if (units == "in2/s") return value * m_from_in * m_from_in;
    else if (units == "in2/min") return value * m_from_in * m_from_in / s_from_min;
    else if (units == "ft2/s") return value * m_from_ft * m_from_ft;
    else if (units == "ft2/min") return value * m_from_ft * m_from_ft / s_from_min;

    else return OpenSMOKE::FatalErrorMessage("Wrong velocity units: " + units);
}

double Mass_Flow_Rate(const double value, const std::string units)
{
    if (units == "kg/s") return value;
    else if (units == "kg/min") return value / s_from_min;
    else if (units == "kg/hr") return value / s_from_hr;
    else if (units == "g/s") return value * kg_from_g;
    else if (units == "g/min") return value * kg_from_g / s_from_min;
    else if (units == "g/hr") return value * kg_from_g / s_from_hr;
    else if (units == "lb/s") return value * kg_from_lb;
    else if (units == "lb/min") return value * kg_from_lb / s_from_min;
    else if (units == "lb/hr") return value * kg_from_lb / s_from_hr;

    else return OpenSMOKE::FatalErrorMessage("Wrong mass flow rate units: " + units);
}

double Mole_Flow_Rate(const double value, const std::string units)
{
    if (units == "kmol/s") return value;
    else if (units == "kmol/min") return value / s_from_min;
    else if (units == "kmol/hr") return value / s_from_hr;
    else if (units == "mol/s") return value * kmol_from_mol;
    else if (units == "mol/min") return value * kmol_from_mol / s_from_min;
    else if (units == "mol/hr") return value * kmol_from_mol / s_from_hr;

    else return OpenSMOKE::FatalErrorMessage("Wrong mole flow rate units: " + units);
}

double Volumetric_Flow_Rate(const double value, const std::string units)
{
    if (units == "m3/s") return value;
    else if (units == "m3/min") return value / s_from_min;
    else if (units == "m3/hr") return value / s_from_hr;
    else if (units == "l/s") return value * m3_from_l;
    else if (units == "l/min") return value * m3_from_l / s_from_min;
    else if (units == "l/hr") return value * m3_from_l / s_from_hr;
    else if (units == "cm3/s") return value * m3_from_cm3;
    else if (units == "cm3/min") return value * m3_from_cm3 / s_from_min;
    else if (units == "cm3/hr") return value * m3_from_cm3 / s_from_hr;
    else if (units == "mm3/s") return value * m3_from_mm3;
    else if (units == "mm3/min") return value * m3_from_mm3 / s_from_min;
    else if (units == "mm3/hr") return value * m3_from_mm3 / s_from_hr;
    else if (units == "ft3/s") return value * m3_from_ft3;
    else if (units == "ft3/min") return value * m3_from_ft3 / s_from_min;
    else if (units == "ft3/hr") return value * m3_from_ft3 / s_from_hr;

    else return OpenSMOKE::FatalErrorMessage("Wrong volumetric flow rate units: " + units);
}

double Heat_Flux(const double value, const std::string units)
{
    if (units == "W/m2") return value;
    else if (units == "J/m2/s") return value;
    else if (units == "J/cm2/s") return value / m2_from_cm2;
    else if (units == "cal/m2/s") return value * J_from_cal;
    else if (units == "cal/cm2/s") return value * J_from_cal / m2_from_cm2;
    else if (units == "kW/m2") return value * W_from_kW;
    else if (units == "kJ/m2/s") return value * J_from_kJ;
    else if (units == "kcal/m2/s") return value * J_from_kcal;
    else if (units == "kWh/m2/s") return value * J_from_kWh;

    else return OpenSMOKE::FatalErrorMessage("Wrong heat flux units: " + units);
}

double Heat_Exchange_Coefficient(const double value, const std::string units)
{
    if (units == "W/m2/K") return value;
    else if (units == "J/m2/K/s") return value;
    else if (units == "J/cm2/K/s") return value / m2_from_cm2;
    else if (units == "cal/m2/K/s") return value * J_from_cal;
    else if (units == "cal/cm2/K/s") return value * J_from_cal / m2_from_cm2;
    else if (units == "kW/m2/K") return value * W_from_kW;
    else if (units == "kJ/m2/K/s") return value * J_from_kJ;
    else if (units == "kcal/m2/K/s") return value * J_from_kcal;
    else if (units == "kWh/m2/K/s") return value * J_from_kWh;

    else return OpenSMOKE::FatalErrorMessage("Wrong heat exchange coefficient units: " + units);
}

double Dynamic_Viscosity(const double value, const std::string units)
{
    if (units == "kg/m/s") return value;
    else if (units == "Pa.s") return value;
    else if (units == "g/cm/s") return value * kg_from_g / m_from_cm;
    else if (units == "kg/cm/s") return value / m_from_cm;

    else return OpenSMOKE::FatalErrorMessage("Wrong dynamic viscosity units: " + units);
}

double Density(const double value, const std::string units)
{
    if (units == "kg/m3") return value;
    else if (units == "g/cm3") return value * kg_from_g / m3_from_cm3;
    else if (units == "lb/ft3") return value * kg_from_lb / m3_from_ft3;
    else if (units == "oz/ft3") return value * kg_from_oz / m3_from_ft3;

    else return OpenSMOKE::FatalErrorMessage("Wrong density units: " + units);
}

double Angle(const double value, const std::string units)
{
    if (units == "rad") return value;
    else if (units == "deg") return value * 2. * PhysicalConstants::pi / 360.;
    else return OpenSMOKE::FatalErrorMessage("Wrong angle units: " + units);
}

double Angular_Velocity(const double value, const std::string units)
{
    if (units == "rad/s") return value;
    else if (units == "rad/min") return value / 60.;
    else if (units == "rad/hr") return value / 3600.;
    else if (units == "rad/ms") return value * 1000.;

    else if (units == "deg/s") return value * 2. * PhysicalConstants::pi / 360.;
    else if (units == "deg/min") return value * 2. * PhysicalConstants::pi / 360. / 60.;
    else if (units == "deg/hr") return value * 2. * PhysicalConstants::pi / 360. / 3600.;
    else if (units == "deg/ms") return value * 2. * PhysicalConstants::pi / 360. * 1000.;

    else if (units == "rpm") return value * 2. * PhysicalConstants::pi / 60.;

    else return OpenSMOKE::FatalErrorMessage("Wrong angular velocity units: " + units);
}

}
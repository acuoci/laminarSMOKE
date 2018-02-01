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
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
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

#include "math/PhysicalConstants.h"
#include "boost/date_time/posix_time/posix_time.hpp"

namespace OpenSMOKE
{
	template<typename Species>
	PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport() 
	{
	}

	template<typename Species>
	PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport(const PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport& orig) {
	}

	template<typename Species>
	PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::~PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport() 
	{
		species_.clear();
		names_.clear();
	}


	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::Allocate()
	{
		NC = boost::lexical_cast<unsigned int>(species_.size());

		ChangeDimensions(NC,&MW,true);
		ChangeDimensions(NC,&uMW,true);

		// Thermodynamic properties
		ChangeDimensions(NC,&Cp,true);
		ChangeDimensions(NC,&Cv,true);

		// Molecular wights [kg/kmol]
		for(unsigned int i=1;i<=NC;i++)
		{
			MW[i] = species_[i-1].MolecularWeight();
			uMW[i] = 1./MW[i];
		}
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::Setup()
	{
		Allocate();

		return true;
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::SpeciesCp(const double T)
	{
		for(unsigned int i=1;i<=NC;i++)
			Cp[i] = species_[i-1].cp(T);
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::SpeciesCv(void)
	{
		for(unsigned int k=1;k<=NC;k++)
			Cv[k] = Cp[k] - PhysicalConstants::R_J_kmol*uMW[k];
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteThermodynamicDataOnASCIIFileOldStyle(std::ofstream &fOutput) const
	{
		std::cout << " * Writing the interpreted thermodynamic and transport properties in ASCII format for old versions of OpenSMOKE..." << std::endl;

		std::string preprocessing_name   = "CHEMKIN-Preprocessor";
		std::string preprocessing_author = "CRECK-Modeling";
		std::string preprocessing_place  = "PoliMI";
		std::string preprocessing_date   = "Mar2013";

		double TMIN = 200.;
		double TMAX = 5000.;

		// Header
		fOutput << "V101116" << std::endl;
		fOutput << preprocessing_name << std::endl;
		fOutput << preprocessing_author << std::endl;
		fOutput << preprocessing_place << std::endl;
		fOutput << preprocessing_date << std::endl;

		// Minimum and Maximum temperatures
		fOutput << TMIN << std::endl;
		fOutput << TMAX << std::endl;

		// Species
		fOutput << NC << std::endl;
		for (unsigned int i=1;i<=NC;i++)
			fOutput << names_[i] << std::endl;

		// NASA Coefficients
		for (unsigned int i=1;i<=NC;i++)
			species_[i-1].WriteOnASCIIFileOldStyle(fOutput);

		// Molecular weights
		fOutput << NC << std::endl;
		for (unsigned int i=1;i<=NC;i++)
			fOutput << species_[i-1].MolecularWeight() << " ";
		fOutput << std::endl;

		// T low
		fOutput << NC << std::endl;
		for (unsigned int i=1;i<=NC;i++)
			fOutput << species_[i-1].Tlow() << " ";
		fOutput << std::endl;

		// T medium
		fOutput << NC << std::endl;
		for (unsigned int i=1;i<=NC;i++)
			fOutput << species_[i-1].Tmedium() << " ";
		fOutput << std::endl;

		// T high
		fOutput << NC << std::endl;
		for (unsigned int i=1;i<=NC;i++)
			fOutput << species_[i-1].Thigh() << " ";
		fOutput << std::endl;

		return true;
	}


	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteThermodynamicsDataOnASCIIFile(std::ofstream &fOutput) const
	{
		// NASA Coefficients
		for (unsigned int i=1;i<=NC;i++)
			species_[i-1].WriteOnASCIIFile(fOutput);

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteThermodynamicsDataOnXMLFile(std::stringstream& xml_string) const
	{	
		xml_string << "<Thermodynamics type=\"NASA\">" << std::endl;
		xml_string << "<NASA-coefficients>" << std::endl;
		
		for (unsigned int i=1;i<=NC;i++)
			species_[i-1].WriteOnASCIIFile(xml_string);

		xml_string << "</NASA-coefficients>" << std::endl;
		xml_string << "</Thermodynamics>" << std::endl;

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteSpeciesOnXMLFile(std::stringstream& xml_string) const
	{	
		xml_string << "<NumberOfSpecies>" << std::endl;
		xml_string << NC << std::endl;
		xml_string << "</NumberOfSpecies>" << std::endl;
		
		xml_string << "<NamesOfSpecies>" << std::endl;
		for (unsigned int i=1;i<=NC;i++)
		{
			xml_string << names_[i] << " ";
			if (i%20==0) xml_string << std::endl;
		}
		xml_string << std::endl;
		xml_string << "</NamesOfSpecies>" << std::endl;
		
		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::CheckThermodynamicConsistency(std::ostream& flog)
	{
		unsigned int corrections = 0;
		unsigned int inconsistencies = 0;

		for(unsigned int i=0;i<species_.size();i++)
		{
			int tag = species_[i].CheckThermodynamicConsistency(flog);
			if ( tag == -1)	inconsistencies++;
			else if ( tag == 0)	corrections++;
		}

		if (inconsistencies > 0)
		{
			std::cout << std::endl;
			std::cout << " ! WARNING MESSAGE: Inconsistencies were found in the thermodynamic properties for " << inconsistencies << " species." << std::endl;
			std::cout << "                    Please check the log file for additional details." << std::endl;
			std::cout << std::endl;
		}
		else if (corrections > 0)
		{
			std::cout << std::endl;
			std::cout << " @ Info message: The thermodynamic coefficients were corrected for " << corrections << " species (relative errors < 1.e-6)" << std::endl;
			std::cout << "                 You can ignore this message." << std::endl;
		}

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::ReformulationOfThermodynamics(const std::string file_name, const std::string original_file_name)
	{
		boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
		std::stringstream date_today;
		date_today << static_cast<int>(now.date().month()) << "/" << now.date().day() << "/" << now.date().year();

		std::ofstream flog(file_name.c_str(), std::ios::out);
		CheckIfFileIsOpen(flog, file_name);

		flog << "! This thermodynamic database was obtained by fitting the thermodynamic properties" << std::endl;
		flog << "! extracted from the following file: " << original_file_name << std::endl;
		flog << "! The thermodynamic properties are fitted in order to preserve not only the " << std::endl;
		flog << "! continuity of each function at the intermediate temperature, but also the  continuity" << std::endl;
		flog << "! of the derivatives, from the 1st to the 3rd order" << std::endl;
		flog << "! The intermediate temperatures are chosen in order to minimize the fitting error." << std::endl;
		flog << "! Last update: " << date_today.str() << std::endl;
		flog << std::endl;
		
		flog << "THERMO ALL" << std::endl;
		flog << "270.   1000.   3500. " << std::endl;
		for(unsigned int i=0;i<species_.size();i++)
			species_[i].ReformulationOfThermodynamics(flog, 1, 1000., 3500.);
		flog << "END" << std::endl;
		flog.close();

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::ReformulationOfThermodynamicsFixedIntermediateTemperature(const std::string file_name, const std::string original_file_name)
	{
		boost::posix_time::ptime now = boost::posix_time::second_clock::local_time();
		std::stringstream date_today;
		date_today << static_cast<int>(now.date().month()) << "/" << now.date().day() << "/" << now.date().year();

		std::ofstream flog(file_name.c_str(), std::ios::out);
		CheckIfFileIsOpen(flog, file_name);

		flog << "! This thermodynamic database was obtained by fitting the thermodynamic properties" << std::endl;
		flog << "! extracted from the following file: " << original_file_name << std::endl;
		flog << "! The thermodynamic properties are fitted in order to preserve not only the " << std::endl;
		flog << "! continuity of each function at the intermediate temperature, but also the  continuity" << std::endl;
		flog << "! of the derivatives, from the 1st to the 3rd order" << std::endl;
		flog << "! The intermediate temperatures are the same for all the species." << std::endl;
		flog << "! Last update: " << date_today.str() << std::endl;
		flog << std::endl;

		flog << "THERMO ALL" << std::endl;
		flog << "270.   1000.   3500. " << std::endl;
		for (unsigned int i = 0; i<species_.size(); i++)
			species_[i].ReformulationOfThermodynamicsFixedIntermediateTemperature(flog, 1, 1000., 3500.);
		flog << "END" << std::endl;
		flog.close();

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::StatusOfThermodynamics(const std::string file_name)
	{
		std::ofstream fout(file_name.c_str(), std::ios::out);
		CheckIfFileIsOpen(fout, file_name);

		fout << std::setfill('=') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ') << "   Status of specific heat"  << std::endl;
		fout << std::setfill('-') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ');

		fout << std::setw(18) << std::left << "Name";
		fout << std::setw(9) << std::left  << "Tlow[K]";
		fout << std::setw(9) << std::left  << "Tmed[K]";
		fout << std::setw(9) << std::left  << "Thigh[K]";

		fout << std::setw(18) << std::left  << "Cp/R(-)";
		fout << std::setw(18) << std::left  << "Cp/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "dCp/R(-)";
		fout << std::setw(18) << std::left  << "dCp/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d2Cp/R(-)";
		fout << std::setw(18) << std::left  << "d2Cp/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d3Cp/R(-)";
		fout << std::setw(18) << std::left  << "d3Cp/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d4Cp/R(-)";
		fout << std::setw(18) << std::left  << "d4Cp/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";
		fout << std::endl;

		fout << std::setfill('=') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ');

		for(unsigned int i=0;i<species_.size();i++)
			int tag = species_[i].ReportStatusSpecificHeatCoefficient(fout);

		fout << std::endl << std::endl;


		fout << std::setfill('=') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ') << "   Status of enthalpy"  << std::endl;
		fout << std::setfill('-') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ');

		fout << std::setw(18) << std::left << "Name";
		fout << std::setw(9) << std::left  << "Tlow[K]";
		fout << std::setw(9) << std::left  << "Tmed[K]";
		fout << std::setw(9) << std::left  << "Thigh[K]";

		fout << std::setw(18) << std::left  << "H/RT(-)";
		fout << std::setw(18) << std::left  << "H/RT(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "dH/RT(-)";
		fout << std::setw(18) << std::left  << "dH/RT(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d2H/RT(-)";
		fout << std::setw(18) << std::left  << "d2H/RT(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d3H/RT(-)";
		fout << std::setw(18) << std::left  << "d3H/RT(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d4H/RT(-)";
		fout << std::setw(18) << std::left  << "d4H/RT(+)";
		fout << std::setw(18) << std::left  << "error(%)";
		fout << std::endl;

		fout << std::setfill('=') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ');

		for(unsigned int i=0;i<species_.size();i++)
			int tag = species_[i].ReportStatusEnthalpy(fout);

		fout << std::endl << std::endl;

		fout << std::setfill('=') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ') << "   Status of entropy"  << std::endl;
		fout << std::setfill('-') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ');

		fout << std::setw(18) << std::left << "Name";
		fout << std::setw(9) << std::left  << "Tlow[K]";
		fout << std::setw(9) << std::left  << "Tmed[K]";
		fout << std::setw(9) << std::left  << "Thigh[K]";

		fout << std::setw(18) << std::left  << "S/R(-)";
		fout << std::setw(18) << std::left  << "S/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "dS/R(-)";
		fout << std::setw(18) << std::left  << "dS/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d2S/R(-)";
		fout << std::setw(18) << std::left  << "d2S/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d3S/R(-)";
		fout << std::setw(18) << std::left  << "d3S/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";

		fout << std::setw(18) << std::left  << "d4S/R(-)";
		fout << std::setw(18) << std::left  << "d4S/R(+)";
		fout << std::setw(18) << std::left  << "error(%)";
		fout << std::endl;
		fout << std::setfill('=') << std::setw(320) << "" << std::endl;
		fout << std::setfill(' ');

		for(unsigned int i=0;i<species_.size();i++)
			int tag = species_[i].ReportStatusEntropy(fout);

		fout << std::endl;
		fout << "================================================================================================================================" << std::endl;
		fout << "   Anomalous behaviour of thermodynamic properties"  << std::endl;
		fout << "================================================================================================================================" << std::endl;

		for(unsigned int i=0;i<species_.size();i++)
			int tag = species_[i].CheckForAnomalies(fout);

		fout.close();

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteThermodynamicsForMATLAB(const std::string file_name)
	{
		std::ofstream fout(file_name.c_str(), std::ios::out);
		CheckIfFileIsOpen(fout, file_name);

		fout << species_.size() << std::endl;
		for(unsigned int i=0;i<species_.size();i++)
			fout << species_[i].name() << " ";
	
		fout.setf(std::ios::scientific);
		fout.precision(10);

		for(unsigned int i=0;i<species_.size();i++)
		{
			OpenSMOKEVectorDouble thermodynamics_parameters;
			thermodynamics_parameters = species_[i].thermodynamics_parameters();
			for(unsigned int k=1;k<=17;k++)
				fout << std::setw(18) << std::left << thermodynamics_parameters[k] << " ";
			fout << std::endl;
		}

		fout.close();

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteThermodynamicCoefficientsOnASCIIFile(const std::string file_name) const
	{
		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);

		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                     SPECIFIC HEAT COEFFICIENTS                                                                     " << std::endl;
		fOutput << std::endl;
		fOutput << "             Cp = A + B*T + C*T^2 + D*T^3 + E*T^4  [J/kg/K]                                                         " << std::endl;
		fOutput << std::endl;
		fOutput << "    Species                      Cp(298K)         Cp(1000K)          LT-HT            A(LT)            B(LT)            C(LT)            D(LT)            E(LT)             A(HT)           B(HT)            C(HT)            D(HT)           E(HT)      " << std::endl;
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
 
		for(unsigned int i=0;i<NC;i++)
		{
			const double Cp298   = species_[i].cp(298.);
			const double Cp1000  = species_[i].cp(1000.);
			
			fOutput << std::right << std::setw(5) << i+1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << names_[i+1];

			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << Cp298;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << Cp1000;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].Tmedium();

			for(unsigned int j=11;j<=15;j++)
				fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[j]*species_[i].RGAS_over_MW();

			for(unsigned int j=4;j<=8;j++)
				fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[j]*species_[i].RGAS_over_MW();
		
			fOutput << std::endl;
		}
		fOutput << std::endl;
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl << std::endl;


		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                           ENTHALPY COEFFICIENTS                                                                     " << std::endl;
		fOutput << std::endl;
		fOutput << "             H/(RT) = A + B*T + C*T^2 + D*T^3 + E*T^4 + F/T  [-]                                                         " << std::endl;
		fOutput << std::endl;
		fOutput << "    Species                   H(298K)[J/kmol]  H(1000K)[J/kmol]      LT-HT            A(LT)            B(LT)            C(LT)            D(LT)            E(LT)           F(LT)             A(HT)            B(HT)            C(HT)           D(HT)             E(HT)            F(HT)      " << std::endl;
		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
 
		for(unsigned int i=0;i<NC;i++)
		{
			const double H298   = species_[i].enthalpy(298.)  * MW[i+1];	// [J/kmol]
			const double H1000  = species_[i].enthalpy(1000.) * MW[i+1];	// [J/kmol]
			

			fOutput << std::right << std::setw(5) << i+1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << names_[i+1];

			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << H298;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << H1000;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].Tmedium();

			// LT
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[11];
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[12]/2.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[13]/3.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[14]/4.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[15]/5.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[16];

			// HT
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[4];
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[5]/2.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[6]/3.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[7]/4.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[8]/5.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[9];
		
			fOutput << std::endl;
		}
		fOutput << std::endl;
		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl << std::endl;

		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                             ENTROPY COEFFICIENTS                                                                     " << std::endl;
		fOutput << std::endl;
		fOutput << "             S/R = A*lnT + B*T + C*T^2 + D*T^3 + E*T^4 + F  [-]                                                         " << std::endl;
		fOutput << std::endl;
		fOutput << "    Species                 S(298K)[J/kmol/K] S(1000K)[J/kmol/K]     LT-HT            A(LT)            B(LT)            C(LT)            D(LT)            E(LT)           F(LT)             A(HT)            B(HT)            C(HT)           D(HT)             E(HT)            F(HT)      " << std::endl;
		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
 
		for(unsigned int i=0;i<NC;i++)
		{
			const double S298   = species_[i].entropy(298.)  * MW[i+1];		// [J/kmol/K]
			const double S1000  = species_[i].entropy(1000.) * MW[i+1];		// [J/kmol/K]

			fOutput << std::right << std::setw(5) << i+1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << names_[i+1];

			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << S298;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << S1000;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].Tmedium();

			// LT
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[11];
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[12];
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[13]/2.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[14]/3.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[15]/4.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[17];

			// HT
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[4];
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[5];
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[6]/2.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[7]/3.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[8]/4.;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << species_[i].thermodynamics_parameters()[10];
		
			fOutput << std::endl;
		}
		fOutput << std::endl;
		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl << std::endl;

		fOutput.close();

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteThermodynamicTablesOnASCIIFile(const std::string file_name) const
	{
		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);

		// Print table
		fOutput << "-----------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "Name                              MW       Cp@298K        H@298K        G@298K        S@298K      Cp@1000K       H@1000K       G@1000K       S@1000K " << std::endl;
	    fOutput << "                           [kg/kmol]   [cal/mol/K]    [kcal/mol]    [kcal/mol]   [cal/mol/K]   [cal/mol/K]    [kcal/mol]    [kcal/mol]   [cal/mol/K] " << std::endl;
		fOutput << "-----------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i<NC; i++)
			species_[i].ThermodynamicsStatus(fOutput);
		fOutput << "-----------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;

		// Print individual species
		std::vector<double> T(12);
		T[0]=300.;	T[1]=600.;
		for(unsigned int i=2;i<12;i++)
			T[i]=T[i-1]+200.;
		for(unsigned int i=0;i<NC;i++)
			species_[i].ThermodynamicsStatus(fOutput, T);

		fOutput.close();

		return true;
	}


	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteTransportDataOnASCIIFileOldStyle(std::ofstream &fOutput) const
	{
		std::cout << " ! Fatal error: Old-style idealgas file requires transport data!" << std::endl;
		return false;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteTransportDataOnASCIIFile(std::ofstream &fOutput) const
	{
		fOutput << "no-transport-data" << std::endl;
		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteTransportDataOnXMLFile(std::stringstream &xml_string) const
	{
		xml_string << "<Transport type=\"none\">" << std::endl;
		xml_string << "</Transport>" << std::endl;

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>::WriteTransportTableOnASCIIFile(std::ostream& fOutput) const
	{
		return true;
	}
}

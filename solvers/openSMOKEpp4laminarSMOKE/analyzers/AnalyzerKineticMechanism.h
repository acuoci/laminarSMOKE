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

#ifndef OpenSMOKE_AnalyzerKineticMechanism_H
#define	OpenSMOKE_AnalyzerKineticMechanism_H

namespace OpenSMOKE
{

	//!  A class to analyze a pre-processed kinetic mechanism 
	/*!
			This class provides the tools to analyze a kinetic mechanism which was previously
			pre-processed. The class requires a xml file containing the kinetic mechanism to be 
			analyzed. The current version performs only the two following operations, but the idea
			is to extend the class with additional, useful tools:

			i. writes reaction tables on ascii files, containing forward and backward kinetic constants,
			   equilibrium constants and thermodynamic data for each reaction as a function of the temperature
			ii. analyzes the reverse reactions (if any) and tries to write the reverse kinetic constants
			    using the Arrhenius' law, by performing a linear regression analysis on a specific 
				temperature interval
	*/

	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	class AnalyzerKineticMechanism
	{
	public:

		/**
		* Default constructor
		* @param kinetics_preprocessor preprocessor containing all the data about the kinetic scheme 
		* @param kinetics_map kinetics_map the kinetic map used toevaluate the reaction and formation rates
		*/
		AnalyzerKineticMechanism(Kinetics_PreProcessor& kinetics_preprocessor, Kinetics_Map& kinetics_map);

		/**
		* This function writes reaction tables on ascii files, containing forward and backward kinetic constants,
		  equilibrium constants and thermodynamic data for each reaction as a function of the temperature
		* @param file_name the name of the file where the reaction tables will be written
		* @param list_of_temperatures list of temperatures at which the detailed kinetic data are written (in K)
		*/
		bool WriteReactionTablesOnASCIIFile(const std::string& file_name, const std::vector<double> list_of_temperatures) const;

		/**
		* This function analyzes the reverse reactions (if any) and tries to write the reverse kinetic constants
		   using the Arrhenius' law, by performing a linear regression analysis on a specific 
		   temperature interval
		* @param file_name the name of the file where the fitted kinetic constants will be written
		*/
		bool WriteFittedInverseKineticConstantsOnASCIIFile(const std::string& file_name) const;

		/**
		* This function analyzes the Chebyshev reactions (if any) and tries to write the kinetic constant
		using the Arrhenius' law, by performing a linear regression analysis on the given T and P ranges
		* @param file_name the name of the file where the fitted kinetic constants will be written
		*/
		bool WriteFittedChebyshevOnASCIIFile(const std::string& file_name) const;

		/**
		* This function writes on a file the results of the collision rate analysis carried
		  out on all the bimolecular reactions to recognize unfeasible reactions
		* @param transport map containing the description of transport properties
		* @param file_name the name of the file where the analysis of collision rates is summarized
		* @param list_of_temperatures list of temperatures (in K) at which the analysis is carried out
		*/
		bool WriteCollisionRatesOnASCIIFile(TransportPropertiesMap_CHEMKIN& transport, const std::string& file_name, const std::vector<double> list_of_temperatures) const;

		/**
		* TODO
		* @param file_name the name of the file where the results are written
		*/
		bool SparsityPatternAnalysis(const std::string& file_name) const;
		

	private:

		Kinetics_PreProcessor& kinetics_preprocessor_;	//!< reference to the kinetic preprocessor 
		Kinetics_Map& kinetics_map_;					//!< reference to the kinetic map 
	};

}

#include "AnalyzerKineticMechanism.hpp"

#endif	/* OpenSMOKE_AnalyzerKineticMechanism_H */


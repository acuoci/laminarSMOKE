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

#ifndef OpenSMOKE_AbstractionReactions_H
#define	OpenSMOKE_AbstractionReactions_H

#include "AbstractionExploded.h"

namespace OpenSMOKE
{
	//!  Abstract class for managing abstraction reactions from a CHEMKIN input file
	/*!
		 This class provides the tools for managing abstraction reactions from a CHEMKIN input file
	*/

	class AbstractionReactions 
	{
	public:

		/**
		* Default constructor
		*/
		AbstractionReactions(InputFileCHEMKIN& kinetics);

		/**
		* Checks if the different sections (ABSTRACTORS and CORRECTIONS) are correctly written in the CHMKIN file
		*/
		bool Checking(InputFileCHEMKIN& kinetics);

		/**
		* Reads the ABSTRACTORS and CORRECTIONS tables
		*/
		bool ExtractTables();

		/**
		* Prints a summary on a stream
		*/
		void Summary(std::ostream& out);

		/**
		* Parses a single reaction to check if it is consistent with the abstraction formalism
		*/
		int ParseReaction(ReactionPolicy_CHEMKIN& reaction, const std::vector<std::string>& map_of_species);

		/**
		* Checks if the reactions are already written in the kinetic mechansim
		*/
		unsigned int RemoveExistingReactions(ReactionPolicy_CHEMKIN& reaction, const std::vector<std::string>& map_of_species);

		/**
		* Returns a vector to each class of abstraction reactions
		*/
		const std::vector<AbstractionExploded>& exploded() const { return exploded_;  }

		/**
		* Parses and cleans the reaction lines
		*/
		int ReplaceAbstractionReaction(std::string& line);

		/**
		* Checks if all the species included in the ABSTRACTORS tables are available in the main list
		*/
		bool CheckListOfSpecies(const std::vector<std::string> names_of_species);

		/**
		* Returns true if the abstraction reactions are active
		*/
		bool is_active() const { return is_active_;  }

		/**
		* Explodes all the classes of abstraction reactions
		*/
		void ExplodeReactions();
    
	private:

		/**
		* Reads the ABSTRACTORS table
		*/
		bool ExtractTableAbstractors();

		/**
		* Reads the CORRECTIONS table
		*/
		bool ExtractTableCorrections();

		/**
		* Reads the ABSTRACTIONS table
		*/
		bool ExtractTableAbstractions();

	private:

		unsigned int first_line_of_abstractors_;				//!< starting line of ABSTRACTORS table
		unsigned int first_line_of_corrections_;				//!< starting line of CORRECTIONS table
		unsigned int first_line_of_abstractions_;				//!< starting line of ABSTRACTIONS table
		std::vector<unsigned int> abstractions_end_lines_;		//!< final line of ABSTRACTORS table


		std::vector<std::string> lines_abstractors_;			//!< lines corresponding to the ABSTRACTORS table
		std::vector<std::string> lines_corrections_;			//!< lines corresponding to the CORRECTIONS table
		std::vector<std::string> lines_abstractions_;			//!< lines corresponding to the ABSTRACTIONS table

		std::vector<std::string> corrections_types_;			//!< list of correction types
		std::vector<double> corrections_A_;						//!< list of correction coefficients for A (read from CHEMKIN input file)
		std::vector<double> corrections_E_;						//!< list of correction coefficients for E (read from CHEMKIN input file) 
		std::vector<double> corrections_Er_;					//!< list of correction coefficients for Er (calculated internally) 

		std::vector<std::string> abstractors_;					//!< list of abstractors
		std::vector<std::string> abstracted_;					//!< list of abstracted species
		std::vector<double> A_;									//!< list of uncorrected frequency factors
		std::vector<double> Beta_;								//!< list of temperature exponents
		std::vector<double> E_;									//!< list of uncorrected activation energies

		double Tref_;											//!< reference temperature
		double alpha_;											//!< fitting parameter
		double Eref_;											//!< reference activation energy

		std::vector<AbstractionExploded> exploded_;				//!< list of abstraction reaction classes

		bool is_active_;										//!< true if the abstraction formalism is on
	};

}

#include "AbstractionReactions.hpp"

#endif	/* OpenSMOKE_AbstractionReactions_H */


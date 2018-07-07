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
|	License                                                           |
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

#ifndef OpenSMOKE_AbstractionExploded_H
#define	OpenSMOKE_AbstractionExploded_H

namespace OpenSMOKE
{
	//!  Class for managing a single class of abstraction reactions 
	/*!
	This class provides the tools for managing a single class of abstraction reactions: R+R'H=>R'+RH
	*/
	class AbstractionExploded
	{

	public:

		/**
		* Default constructor
		*/
		AbstractionExploded(	const std::vector<std::string>& species_R, const std::string& species_RpH,
								const std::string& species_Rp, const std::vector<std::string>& species_RH,
								const double n_H, const unsigned int type_H);

		/**
		* Explodes the reactions for the given class
		*/
		void ExplodeReactions(	const std::vector<std::string>& abstractors,
								const std::vector<double>& abstractors_A,
								const std::vector<double>& abstractors_Beta,
								const std::vector<double>& abstractors_E,
								const std::vector<double>& corrections_A,
								const std::vector<double>& corrections_E,
								const std::vector<double>& corrections_Er,
								const double Eref, const double alpha, const double Tref);

		/**
		* Checks if one or more exploded abstraction reactions are already present in the kinetic mechanism
		*/
		unsigned int RemoveExistingReactions(const std::vector<std::string>& reactant_species, const std::vector<std::string>& product_species);

		/**
		* Prints the exploded list of reactions on a output stream
		*/
		void PrintExplodedReactions(std::ostream& out) const;

		/**
		* Returns the number of hydrogen for the given abstraction class
		*/
		double n_H() const { return n_H_; }

		/**
		* Returns the type of hydrogen for the given abstraction class (from 0 to N)
		*/
		double type_H() const { return type_H_; }

		/**
		* Returns the RH species for the given class of abstraction reactions
		*/
		const std::vector<std::string>& species_RH() const { return species_RH_; }

		/**
		* Returns the R'H species for the given class of abstraction reactions
		*/
		const std::string species_RpH() const { return species_RpH_; }

	public:

		std::vector<double> A_;						//!< the corrected frequency factors [mol,l,s]
		std::vector<double> Beta_;					//!< the exponent of temperature
		std::vector<double> E_;						//!< the corrected activation energies [cal/mol]
		std::vector<bool> is_dummy_;				//!< true if the reaction is dummy (i.e. products equal reactants)
		std::vector<bool> is_existing_;				//!< true if the reaction is already present in the kinetic mechanism

		std::vector<std::string> species_R_;		//!< list of R species: R+R'H=>R'+RH
		std::string species_RpH_;					//!< R'H species: R+R'H=>R'+RH
		std::string species_Rp_;					//!< R' species: R+R'H=>R'+RH
		std::vector<std::string> species_RH_;		//!< RH species: R+R'H=>R'+RH

		unsigned int n_;							//!< number of exploded reactions
		double n_H_;								//!< number of hydrogens
		unsigned int type_H_;						//!< type of hydrogen

		std::vector<std::string> label_;			//!< explicit reactions
	};
}

#include "AbstractionExploded.hpp"

#endif	/* OpenSMOKE_AbstractionExploded_H */


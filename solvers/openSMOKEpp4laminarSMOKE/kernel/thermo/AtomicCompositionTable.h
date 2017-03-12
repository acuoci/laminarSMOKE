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

#ifndef OpenSMOKE_AtomicCompositionTable_H
#define	OpenSMOKE_AtomicCompositionTable_H

#include <Eigen/Dense>
#include "kernel/thermo/AtomicComposition.h"

namespace OpenSMOKE
{
	//!  A class to manage the atomic composition of all the species in a kinetic scheme
	/*!
			This class basically stores the atmonic composition in a matrix, containing for each 
			species the number of each element in the kinetic scheme. All the elements really
			present among the species in the kinetic scheme are accounted for (i.e. if in input file
			the user declared an atomic element which is not present in any species, this atomic
			element is automatically excluded from the list.
	*/

	class AtomicCompositionTable 
	{

	public:
	
		/**
		*@brief Default Constructor
		*/
		AtomicCompositionTable() {};

		/**
		*@brief Default copy constructor
		*/
		AtomicCompositionTable(const AtomicCompositionTable& orig);

		/**
		*@brief Default copy constructor
		*/
		void CopyFrom(const AtomicCompositionTable& orig);

		/**
		*@brief Default destructor
		*/
		virtual ~AtomicCompositionTable() {};

		/**
		*@brief Given the map of species in the kinetic scheme, the matrix and the vectors
				containing the information about the atomic coefficients are built
		*/
		template<typename SpeciesMap>
		void CreateAtomicElementLists(const SpeciesMap& species);

		/**
		*@brief Given a reaction, this function performs a check about the stoichiometry. 
				If the test fails an error message is returned. A warning message is 
				returned if the test fails for a small relative error
		 @param fOut the stream where to print the error or warning message
		 @param reaction the reaction to be checked
		 @param epsilon maximum allowed relative error
		*/
		template<typename Reaction>
		unsigned int CheckStoichiometry(std::ostream& fOut, const Reaction& reaction, const double epsilon);

		/**
		*@brief Given a reaction for which the elemental balances have small errors in the closure
	            this function performs a least-squares analysis to correct the stoichiometric coefficients
				of products in order to ensure the perfect closure of atomic balances.
		@param reaction the reaction to be checked
		*/
		template<typename Reaction>
		unsigned int CorrectStoichiometry(Reaction& reaction);

		/**
		*@brief Return a list of the names of the elements 
		 @return The list of the elements' names (zero-index)
		*/
		const std::vector<std::string>& element_names_list() const { return element_names_list_;}

		/**
		*@brief Return a the weights of the elements
		 @return The weights of the elements (zero-index)
		*/
		const Eigen::VectorXd& element_weights_list() const { return element_weights_list_; }

		/**
		*@brief Return a the matrix containing the atomic composition
		 @return The matrix containing the atomic composition (zero-index): number of species x number of elements
		*/
		const Eigen::MatrixXd& element_coefficients_list() const { return element_coefficients_list_; }


	protected:

		std::vector<std::string> element_names_list_;		//!< list containing the names of the elements
		Eigen::VectorXd element_weights_list_;				//!< vector containing the weights of the elements [kg/kmol]
		Eigen::MatrixXd element_coefficients_list_;			//!< matrix containing the atomic composition (ns x ne)
	};
}

#include "AtomicCompositionTable.hpp"

#endif	/* OpenSMOKE_AtomicCompositionTable_H */


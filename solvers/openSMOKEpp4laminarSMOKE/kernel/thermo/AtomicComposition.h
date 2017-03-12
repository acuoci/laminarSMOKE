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

#ifndef OpenSMOKE_AtomicComposition_H
#define OpenSMOKE_AtomicComposition_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	//!  A class describing the atomic composition of a single species
	/*!
		 This class provides the information about the atomic composition of a single
		 species (list and number of elements). Specific functions to calculate the
		 atomic weight of the species, and to write the atomic composition on a file
		 are also provided.
	*/

	class AtomicComposition : public OpenSMOKEClass
	{
	
	public:

		/**
		*@brief Default Constructor
		*/
		AtomicComposition();

		/**
		*@brief Copy Constructor
		*/
		AtomicComposition(const AtomicComposition& orig);

		/**
		*@brief Constructor: the atomic composition is read from an external input file
		*@param fInput The input stream from which coefficients have to be read
		*@param fileFormat The format of the input stream 
		*/
		AtomicComposition(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat);
        
		/**
		*@brief Constructor: from a list of elements names and the corresponding coefficients (zero-index)
		*@param names The list of names 
		*@param coefficients The list of coefficients
		*/	
		AtomicComposition( std::vector<std::string> const& names, std::vector<double> const& coefficients );    

		/**
		*@brief Setup function
		*@param names The list of names 
		*@param coefficients The list of coefficients
		*/	
		void operator() ( std::vector<std::string>& names, std::vector<double>& coefficients );  

		/**
		*@brief Assignment operator
		*@param rhs AtomicComposition to be assigned 
		*/	
		AtomicComposition& operator=(const AtomicComposition& rhs);

		/**
		*@brief Calculates the corresponding molecular weight
		*@return The molecular weight in kg/kmol
		*/	
		double mw() const;

		/**
		*@brief Returns the atomic coefficient corresponding to the requested element (with checks)
		*@return The requested coefficient
		*/	
		double GetCoefficient(const std::string name);

		/**
		*@brief Write the elemental composition
		*@param fOut The output stream
		*/	
		void Status(std::ostream& fOut);

		/**
		*@brief Saves the atomic composition on a file
		*@param fOutput The output stream
		*@param fileFormat The format (ascii or binary) of output file
		*/	
		void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat);
		
		/**
		*@brief Loads the atomic composition from a file
		*@param fInput The input stream
		*@param fileFormat The format (ascii or binary) of input file
		*/	
		void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat);

		/**
		*@brief Checks the atomic composition
		*/	
		bool CheckAtomicComposition() const;

		/**
		*@brief Returns the names of the atomic elements
		*/	
		const std::vector<std::string>& element_names() const;

		/**
		*@brief Returns the coefficients of the atomic elements
		*/	
		const std::vector<double>& element_coefficients() const;

	protected:

		std::vector<double>			coefficients_;		//!< element coefficients
		std::vector<std::string>	names_;             //!< element names

	};
}

#include "AtomicComposition.hpp"

#endif // OpenSMOKE_AtomicComposition_H
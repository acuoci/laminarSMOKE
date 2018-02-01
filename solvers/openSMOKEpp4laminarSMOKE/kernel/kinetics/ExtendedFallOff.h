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


#ifndef OpenSMOKE_ExtendedFallOff_H
#define	OpenSMOKE_ExtendedFallOff_H

#include "math/PhysicalConstants.h"
#include <vector>
#include <string.h>

namespace OpenSMOKE
{
	enum ExtendedFalloffType { EXTENDED_FALLOFF_LINDEMANN, EXTENDED_FALLOFF_TROE, EXTENDED_FALLOFF_SRI };

	//!  A class to manage a Pressure Dependence through logarithmic interpolation rate expression
	/*!
		 TODO
	*/

	class ExtendedFallOff
	{
		
		public:

			/**
			* Default constructor
			*/
			ExtendedFallOff();

			/**
			* Default copy constructor
			*/
			ExtendedFallOff(const ExtendedFallOff& orig);

			/**
			* Default destructor
			*/
			//virtual ~ExtendedFallOff();

			/**
			*@brief Prepares all the data to be used for evaluating the reaction rate
			*/
			void Setup(const std::vector< std::vector<std::string> >& coefficients);

			/**
			*@brief Reads from a file the data about the reaction
			*/
			void ReadFromASCIIFile(std::istream& fInput);

			/**
			*@brief Writes a short summary on file
			*/
			void WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const;

			/**
			*@brief Writes the reaction in CHEMKIN on file
			*/
			void WriteCHEMKINOnASCIIFile(std::ostream& fOutput) const;

			/**
			*@brief Evaluates the kinetic constant
			*/
			double KineticConstant(const double T, const double P, const double cTot, const double* c);

			void WriteAdditionalParameters(std::ofstream& fOut, const int i);
			int number_of_species() const { return n_; }
			std::string species(const int i) const { return species_[i]; }
			double A0(const int i) const { return A0_[i]; }
			double Beta0(const int i) const { return Beta0_[i]; }
			double E0(const int i) const { return E0_over_R_[i]*PhysicalConstants::R_J_kmol; }

		private:

			/**
			*@brief Prepares all the data to be used for evaluating the reaction rate
			*/
			void Setup(unsigned int index, std::vector<double> coefficients, const double conversion_A, const double conversion_E);

			/**
			*@brief Writes a short summary on file
			*/
			void WriteShortSummaryOnASCIIFile(const unsigned int index, std::ostream& fOutput) const;

        private:
            
			/**
			*@brief Returns an error message
			*/
			void ErrorMessage(const std::string message);

			/**
			*@brief Returns a warning message
			*/
			void WarningMessage(const std::string message);

			/**
			*@Check if redundant data are provided
			*/
			void CheckForRedundancy(const std::vector< std::vector<std::string> >& coefficients);

			/**
			*@Check if the low-pressure parameters for the mixture are provided
			*/
			void CheckForMixtureLowPressureParameters(const std::vector< std::vector<std::string> >& coefficients);

			/**
			*@Check if the low-pressure parameters are provided for every species
			*/
			void CheckForLowPressureKineticParameters(const std::vector< std::vector<std::string> >& coefficients);

			/**
			*@Check if the third-body efficiencies of specified species are equal to 1
			*/
			void CheckForUnitaryThirdBodyEfficiencies();

	private:

			int n_;									//!< number of species (mixture included)
			std::vector<std::string> species_;		//!< list of bath species
			std::vector<int> species_indices_;		//!< indices of bath species (0-based) (-1 means mixture)
			
			double conversion_A0_;							//!< conversion factor from the CHEMKIN units to the SI units
			double conversion_AInf_;						//!< conversion factor from the CHEMKIN units to the SI units

			double AInf_;									//!< high pressure frequency factor [kmol,m3,s]
			double BetaInf_;								//!< high pressure temperature exponent
			double EInf_over_R_;							//!< high pressure activation temperature [K]
			std::vector<double> A0_;						//!< low pressure frequency factors [kmol,m3,s]
			std::vector<double> Beta0_;						//!< low pressure temperature exponents
			std::vector<double> E0_over_R_;					//!< low pressure activation temperatures [K]
			std::vector<double> k_;							//!< individual kinetic constants [kmol,m3,s]
			std::vector<ExtendedFalloffType > type_;		//!< type of falloff reaction [kmol,m3,s]
			std::vector<std::vector<double> > teta_;		//!< parameters of falloff reactions (Troe, SRI)
			std::vector<double> third_body_efficiencies_;	//!< third body efficiencies
			std::vector<int> third_body_indices_;			//!< indices of thrird body
	};

}

#include "ExtendedFallOff.hpp"

#endif	/* OpenSMOKE_ExtendedFallOff_H */


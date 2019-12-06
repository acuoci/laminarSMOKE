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


#ifndef OpenSMOKE_ChebyshevPolynomialRateExpression_H
#define	OpenSMOKE_ChebyshevPolynomialRateExpression_H

#include "math/PhysicalConstants.h"
#include <vector>
#include <string.h>
#include <Eigen/Dense>

namespace OpenSMOKE
{
	//!  A class to manage a Chebyshev Polynomial Rate Expression
	/*!
		 This class provides the interface and the tools to manage a Chebyshev Polynomial Rate Expression
		 (CHEB in CHEMKIN standard)
		 (see Section 3.6.4 of CHEMKIN manual)
	*/

	class ChebyshevPolynomialRateExpression
	{
	public:

		/**
		* Default constructor
		*/
		ChebyshevPolynomialRateExpression();

		/**
		* Default copy constructor
		*/
		ChebyshevPolynomialRateExpression(const ChebyshevPolynomialRateExpression& orig);

		/**
		* Default destructor
		*/
		//virtual ~ChebyshevPolynomialRateExpression();

		/**
		*@brief Prepares all the data to be used for evaluating the reaction rate
		*/
		void Setup(std::vector<double> coefficients_, std::vector<double> pressures_, std::vector<double> temperatures_);

		/**
		*@brief Sets the policy to be used for managing temperatures and pressures outside the limits
		*/
		void SetViolationAllowed(const bool flag);

		/**
		*@brief Evaluates the kinetic constant
		*/
		double KineticConstant(const double T, const double P);

		/**
		*@brief Evaluates the kinetic constant
		*/
		void FitttingArrheniusLaw(std::ostream& fOutput);

		/**
		*@brief Writes on a file the data about the Chebyshev Polynomial 
		*/
        void WriteStatus(std::ostream& fOutput) const;
		
		/**
		*@brief Writes on a file the data about the Chebyshev Polynomial (compatibility with previous version of OpenSMOKE)
		*/
		void WriteOnASCIIFileOldStyle(std::ostream& fOutput) const;

		/**
		*@brief Reads from a file the data about the Chebyshev Polynomial
		*/		
		void ReadFromASCIIFile(std::istream& fInput);

		/**
		*@brief Writes a short summary on file
		*/	
		void WriteShortSummaryOnASCIIFile(std::ostream& fOutput) const;
                
        private:
            
            	/**
		*@brief Evaluates the Chebyshev Polynomial of the first kind of degree n-1
		*/
		double Phi(const int n, const double x);

		/**
		*@brief Returns an error message
		*/
		void ErrorMessage(const std::string message);

		/**
		*@brief Returns a warning message
		*/
		void WarningMessage(const std::string message);

	private:

		unsigned int N;					//!< number of basis functions along the temperature axis
		unsigned int M;					//!< number of basis functions along the pressure axis
		Eigen::MatrixXd a;				//!< Chebyshev coefficients
		Eigen::VectorXd phi_n;			//!< temperature Chebyshev polynomials
		Eigen::VectorXd phi_m;			//!< pressure Chebyshev polynomials

		double conversion;				//!< conversion factor for the reaction rate

		double Tmin;					//!< minimum temperature [K]
		double Tmax;					//!< maximum temperature [K]
		double Pmin;					//!< minimum pressure [Pa]
		double Pmax;					//!< maximum temperature [Pa]

		double log10_Pmin;				//!< log of minimum pressure (only for effciency reasons)
		double log10_Pmax;				//!< log of maximum pressure (only for effciency reasons)

		bool is_violation_allowed_;		//!< if true, the polynomials are used also outside the limits
	};

}

#include "ChebyshevPolynomialRateExpression.hpp"

#endif	/* OpenSMOKE_ChebyshevPolynomialRateExpression_H */


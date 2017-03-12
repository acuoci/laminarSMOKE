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

#ifndef OpenSMOKE_TransportPolicy_CHEMKIN_H
#define	OpenSMOKE_TransportPolicy_CHEMKIN_H

namespace OpenSMOKE
{
	//!  A class to manage the transport properties of a species in CHEMKIN format
	/*!
		 This class reads and analyze the transport properties in CHEMKIN format
	*/

	class TransportPolicy_CHEMKIN 
	{
		public:

			/**
			*@brief Default constructor
			*/
			TransportPolicy_CHEMKIN();

			/**
			*@brief Default copy constructor
			*/
			TransportPolicy_CHEMKIN(const TransportPolicy_CHEMKIN& orig);

			/**
			*@brief Default destructor
			*/
			virtual ~TransportPolicy_CHEMKIN();

			/**
			*@brief Evaluates the dynamic viscosity of the species in [kg/m/s]
			*@param T Temperature in [K]
			*@param P Pressure in [Pa]
			*/
			double Viscosity(const double T, const double P);

			/**
			*@brief Evaluates the mass diffusivities of the species in [m2/s]
			*@param T Temperature in [K]
			*@param P Pressure in [Pa]
			*@param c Concentration in [kmol/m3]
			*/
			double Diffusivity(const double T, const double P, const OpenSMOKE::OpenSMOKEVectorDouble &c);

			/**
			*@brief Reads the transport properties from a string
			*/
			bool ReadTransportFromStrings(const std::vector<std::string> &lines);

			/**
			*@brief Copy the transport properties to the rhs object
			*/
			void CopyTransportProperties(TransportPolicy_CHEMKIN& rhs) const;

			/**
			*@brief Returns the name name of the species
			*/
			std::string name_transport() const { return name_transport_; }
    
		public:

			/**
			*@brief Returns the shape factor of the species
			        An index indicating whether the molecule has a monatomic, linear or nonlinear geometrical
					configuration. If the index is 0, the molecule is a single atom. If the 
                    index is 1 the molecule is linear, and if it is 2, the molecule is nonlinear.
			*/
			int shape_factor() const;

			/**
			*@brief Returns the Lennard-Jones potential well depth in [K]
			*/
			double epsylon_over_kb() const;

			/**
			*@brief Returns the Lennard-Jones collision diameter in [angstroms]
			*/
			double sigma() const;

			/**
			*@brief Returns the dipole moment in [Debye]
			*/
			double mu() const;

			/**
			*@brief Returns the polarizability in [angstroms^3]
			*/
			double alfa() const;

			/**
			*@brief Returns the rotational relaxation collision number at 298 K
			*/
			double zRot298() const;

			/**
			*@brief Returns data about the transport properties
			*/
			void TransportStatus(std::ostream& fOut) const;
    
		protected:

			std::string name_transport_;							//!< name of the species
			std::vector<int> transport_int_parameters_;				//!< vector of transport properties (integer values)
			std::vector<double> transport_double_parameters_;		//!< vector of transport properties (double values)
	};

}

#include "TransportPolicy_CHEMKIN.hpp"

#endif	/* OpenSMOKE_TransportPolicy_CHEMKIN_H */


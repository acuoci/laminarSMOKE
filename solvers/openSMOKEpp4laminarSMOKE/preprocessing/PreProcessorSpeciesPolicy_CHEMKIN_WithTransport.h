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

#ifndef OpenSMOKE_PreProcessorSpeciesPolicy_CHEMKIN_WithTransport_H
#define	OpenSMOKE_PreProcessorSpeciesPolicy_CHEMKIN_WithTransport_H

#include <Eigen/Dense>
#include "math/OpenSMOKEVector.h"
#include "PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport.h"

namespace OpenSMOKE
{

	//!  A class to preprocess thermodynamic and transport properties provided in CHEMKIN format
	/*!
		 This class preprocess the thermodynamic and transport properties provided
		 in CHEMKIN format. For the transport properties the fitting procedure commonly 
		 adopted by CHEMKIN is performed.
		 The equations are the same reported in the Theory Manual of CHEMKIN Software
	*/

	template<typename Species>
	class PreProcessorSpeciesPolicy_CHEMKIN_WithTransport : public PreProcessorSpeciesPolicy_CHEMKIN_WithoutTransport<Species>
	{
	public:

		/**
		* Default constructor
		*/
		PreProcessorSpeciesPolicy_CHEMKIN_WithTransport();

		/**
		* Copy constructor
		*/
		PreProcessorSpeciesPolicy_CHEMKIN_WithTransport(const PreProcessorSpeciesPolicy_CHEMKIN_WithTransport& orig);
    
		/**
		* Default destructor
		*/
		virtual ~PreProcessorSpeciesPolicy_CHEMKIN_WithTransport();

		/** 
		* This function allocates the memory and prepares the data to be used for
		  the internal preprocessing
		*/
		bool Setup();

		/**
		* Fitting of transport properties
		*/
		bool Fitting();

		/**
		* This function write the transport properties in a readeable format
		  This function can be called only after the transport properties have 
		  been preprocessed
		*/
		bool WriteTransportDataOnASCIIFileOldStyle(std::ofstream &fOutput) const;

		/**
		* This function writes the transport properties in a readeable format
		  to be read by the TransportProperties Maps (obsolete, TOREMOVE)
		*/
		bool WriteTransportDataOnASCIIFile(std::ostream &fOutput) const;

		/**
		* This function writes the transport properties in a readeable XML format
		  to be read by the TransportProperties Maps
		*/
		bool WriteTransportDataOnXMLFile(std::stringstream& xml_string) const;

		/**
		* This function writes the fitting coefficients calculated from the transport
		  properties in a readeable format. This function can be called only after 
		  the transport properties have been preprocessed
		*/
		bool WriteFittingCoefficientsOnASCIIFile(const std::string file_name) const;

		/**
		* This function writes the transport properties (lennard-jones parameters) in a table.
		*/
		bool WriteTransportTableOnASCIIFile(std::ostream& fOutput) const;

		/**
		* This function writes the species bundling result in a readeable XML format
		  to be read by the TransportProperties Maps
		*/
		bool WriteSpeciesBundlingOnXMLFile(std::stringstream& xml_string, const double epsilon) const;


	private:

		void Allocate();
		void InitializeMassDiffusivity();
		void InitializeThermalConductivity();
		void InitializeViscosity();

		void SpeciesViscosities(const double T);								//!< calculates the viscosities of species
		void SpeciesThermalConductivities(const double T);						//!< calculates the thermal conductivites of species
		void SpeciesBinaryDiffusivities(const double T, const double P_bar);	//!< calculates the binary diffusivity coefficients of species
		void SpeciesThermalDiffusionRatios(const double T);						//!< calculates the thermal diffusivities for the selected species

		void WriteViscosityFittingCoefficients(std::ostream &fOutput) const;
		void WriteThermalConductivityFittingCoefficients(std::ostream &fOutput) const;
		void WriteBinaryDiffusivityFittingCoefficients(std::ostream &fOutput) const;
		void WriteThermalDiffusionRatiosFittingCoefficients(std::ostream &fOutput) const;

		void SpeciesBundling(std::ostream &fOutput, const double epsilon) const;


	private:

		/** @name Evaluation of viscosities of single species
		 *  This group provides the functions to calculate the viscosities of single species
		 */
		///@{
		/** calculates the reduced temperature for each species */
		void ReducedTemperatureSingleSpecies(const double T);
		/** calculates the collisional integral */
		void SingleSpeciesViscosity(const double T);
		///@}


		/** @name Evaluation of thermal conductivites of single species
		 *  This group provides the functions to calculate the thermal conductivities of single species
		 */
		///@{
		/** calculates the FT */
		void ComputeFT(const double T);
		/** calculates the Z rotational */
		void ComputeZrotational();
		/** calculates the vibrational Cv */
		void ComputeCvVibrational();
		/** calculates the densities */
		void DensitySingleSpecies(const double T);
		/** calculates the thermal conductivies */
		void ThermalConductivitySingleSpecies(void);
		/** calculates the collisional integral */
		void Omega11kk(void);
		/** calculates the self diffusion coefficient */
		void SelfDiffusionCoefficient(const double T);
		///@}


		/** @name Evaluation of mass diffusivities of single species
		 *  This group provides the functions to calculate the mass diffusivities of single species
		 */
		///@{
		/** calculates the collisional integral */
		void BinaryDiffusionCoefficients(const double T, const double P_bar);
		/** calculates the Lennard Jones Potential And The Collision Diameter */
		void LennardJonesPotentialAndCollisionDiameter(const int j, const int k, double& epsjk_over_kb, double& sigmajk);
		/** calculates the uncorrected Lennard Jones Potential */
		double UncorrectedLennardJonesPotential(const int j, const int k) const;
		///@}

	private:

		OpenSMOKE::OpenSMOKEVectorInt shape_factor;				
		OpenSMOKE::OpenSMOKEVectorDouble epsylon_over_kb;
		OpenSMOKE::OpenSMOKEVectorDouble sigma;
		OpenSMOKE::OpenSMOKEVectorDouble mu;
		OpenSMOKE::OpenSMOKEVectorDouble alfa;
		OpenSMOKE::OpenSMOKEVectorDouble zRot298;

		OpenSMOKE::OpenSMOKEVectorDouble kb_over_epsylon;
		OpenSMOKE::OpenSMOKEVectorDouble muStar;
		OpenSMOKE::OpenSMOKEVectorDouble epsylon;

		// Mass Diffusivity
		OpenSMOKE::OpenSMOKEMatrixDouble deltajkStar;
		OpenSMOKE::OpenSMOKEMatrixDouble coeff_Djk;
		OpenSMOKE::OpenSMOKEMatrixDouble Djk;

		// Viscosity
		OpenSMOKE::OpenSMOKEVectorDouble deltakStar;
		OpenSMOKE::OpenSMOKEVectorDouble omega22k;
		OpenSMOKE::OpenSMOKEVectorDouble TkStar;
		OpenSMOKE::OpenSMOKEVectorDouble eta;
		OpenSMOKE::OpenSMOKEVectorDouble coeff_eta;

		// Thermal conductivity
		OpenSMOKE::OpenSMOKEVectorDouble fVib;
		OpenSMOKE::OpenSMOKEVectorDouble fRot;
		OpenSMOKE::OpenSMOKEVectorDouble fTrans;
		OpenSMOKE::OpenSMOKEVectorDouble cVtrans;
		OpenSMOKE::OpenSMOKEVectorDouble cVrot;
		OpenSMOKE::OpenSMOKEVectorDouble cVvib;
		OpenSMOKE::OpenSMOKEVectorDouble f298;
		OpenSMOKE::OpenSMOKEVectorDouble fT;
		OpenSMOKE::OpenSMOKEVectorDouble zRot;
		OpenSMOKE::OpenSMOKEVectorDouble rho;
		OpenSMOKE::OpenSMOKEVectorDouble Dkk;
		OpenSMOKE::OpenSMOKEVectorDouble coeff_Dkk;
		OpenSMOKE::OpenSMOKEVectorDouble omega11kk;
		OpenSMOKE::OpenSMOKEVectorDouble lambda;

		// Thermal diffusion ratios
		OpenSMOKE::OpenSMOKEVectorInt		iThermalDiffusionRatios;	//!< list of species for which the thermal diffusivity has to be calculated
		OpenSMOKE::OpenSMOKEMatrixDouble	TetaStar;                   //!< thermal diffusivities

		// Fitting coefficients
		Eigen::MatrixXd fittingEta;							//!< matrix containng the fitting coefficients (viscosity) for each species
		Eigen::MatrixXd fittingLambda;						//!< matrix containng the fitting coefficients (thermal conductivity) for each species
	//	Eigen::MatrixXd fittingBinaryDiffusivities;			//!< matrix containng the fitting coefficients (mass diffusivity) for each species
		Eigen::MatrixXd* fittingBinaryDiffusivities;		//!< matrix containng the fitting coefficients (mass diffusivity) for each species
		Eigen::MatrixXd fittingTetaBinary;					//!< matrix containng the fitting coefficients (thermal diffusivity) for each species

	private:

		static const double BOLTZMANN;
		static const double BOLTZMANN3;
		static const double MU_MIN;
		static const double ZROTA;
		static const double ZROTB;
		static const double ZROTC;
		static const double CONST_2SUPI;
		static const double CONST_5SU3R;
	};

}

#include "PreProcessorSpeciesPolicy_CHEMKIN_WithTransport.hpp"

#endif	/* OpenSMOKE_PreProcessorSpeciesPolicy_CHEMKIN_WithTransport_H */


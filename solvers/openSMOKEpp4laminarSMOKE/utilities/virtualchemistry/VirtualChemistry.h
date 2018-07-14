/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Alberto Cuoci, Giampaolo Maio, Benoit Fiorina                |
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
|   Copyright(C) 2018 Alberto Cuoci                                       |
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

#include "Eigen/Dense"
#include "math.h"
#include "LookupTable.h"
#include "Grammar_VirtualChemistry.h"

#ifndef OpenSMOKEpp_VirtualChemistry
#define OpenSMOKEpp_VirtualChemistry

namespace OpenSMOKE
{
	//!  A class to manage Virtual Chemistry Approach
	/*!
	This class provides the tools to manage Virtual Chemistry
	*/

	class VirtualChemistry
	{
	public:

		/**
		*@brief Default constructor
		*/
		VirtualChemistry( OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, bool is_active=true);
		VirtualChemistry(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML, OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Setup from a dictionary
		*@param dictionary name
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		void SetTableMain(const boost::filesystem::path path_table);
		void SetTableCO(const boost::filesystem::path path_table);
		void SetTableNO(	const boost::filesystem::path path_table_1, const boost::filesystem::path path_table_2,
							const boost::filesystem::path path_table_3, const boost::filesystem::path path_table_4,
							const boost::filesystem::path path_table_5);

		void SetReactions(const bool flag);

		void SetVersion(const unsigned int value);
		void SetFuel(const std::string name, const double mw);
		void SetOxidizer(const std::string name, const double mw);
		void SetInert(const std::string name, const double mw);

		void SetTransportProperties(const double mu0, const double T0, const double Beta0_, const double Pr0_);
		void SetLewisNumber(const unsigned int i, const double Le);

		void SetOnTheFlyOptimization(const std::string flag);
		void SetParameters(const Eigen::VectorXd& parameters, const std::vector<bool>& flag);
		void GetParameters(Eigen::VectorXd& parameters, const std::vector<bool>& flag) const;
		void GetOriginalParameters(Eigen::VectorXd& parameters);

		double MWMix(double* Y);
		double MW(const unsigned int j, double* Y);

		void MoleFractions(const double mw, double* Y, double* X);
		double CpMix(const double T, const double P_Pa, double* Y);
		void CpSpecies(const double T, const double P_Pa, double* Cp);
		double HMix(const double T, const double P_Pa, double* Y);
		double Qdot(const double T, const double P_Pa, double* Omega);

		double DynamicViscosity(const double T);
		double ThermalConductivity(const double mu, const double cp);
		void MassDiffusionCoefficients(const double lambda, const double rho, const double Cp, double* Dmix);

		void FormationRates(const double cTot, const double MW, const double T, double* Y, double* Omega);

		bool is_active() const { return (NSVector_.size() != 0); }

		std::string on_the_fly_optimization() const { return on_the_fly_optimization_; }
		
		unsigned int ns() const { return ns_; }
		unsigned int ns_main() const { return ns_main_; }

	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map

		std::string fuel_name_;
		std::string oxidizer_name_;
		std::string inert_name_;

		double fuel_mw_;
		double oxidizer_mw_;
		double inert_mw_;

		unsigned int fuel_index_;
		unsigned int oxidizer_index_;
		unsigned int inert_index_;

		unsigned int I_index_;
		unsigned int P1_index_;
		unsigned int P2_index_;
		unsigned int P3_index_;
		unsigned int P4_index_;

		unsigned int CO_index_;
		unsigned int V1_index_;
		unsigned int V2_index_;

		unsigned int NO_index_;
		unsigned int W1_index_;
		unsigned int W2_index_;
		unsigned int W3_index_;

		unsigned int iVersion_;
		bool iReactions_;
		bool iSubMechanism_CO_;
		bool iSubMechanism_NO_;

		unsigned int ns_;
		unsigned int ns_main_;
		Eigen::VectorXd NSVector_;
		Eigen::VectorXd MW_;

		double mu0_;	// [kg/m/s]
		double T0_;		// [K]
		double Beta0_;
		double Pr0_;
		Eigen::VectorXd Le_;

		LookupTable	table_main_;
		LookupTable table_co_;
		LookupTable table_no_1_;
		LookupTable table_no_2_;
		LookupTable table_no_3_;
		LookupTable table_no_4_;
		LookupTable table_no_5_;

		std::string on_the_fly_optimization_;

		double A1_;
		double A2_;
		double E1_;
		double E2_;
		double nuF_1_;
		double nuOX_1_;
		double nuI_2_;

		// Interpolation (table main)
		double alpha1_;
		double alpha2_;
		double alpha3_;
		double alpha4_;

		double A3_NO_;
		double A4_NO_;
		double A5_NO_;
		double A6_NO_;
		double A7_NO_;

		double E3_NO_;
		double E4_NO_;
		double E5_NO_;
		double E6_NO_;
		double E7_NO_;

		double alpha_NO_;
		double gamma_NO_;
		double beta_NO_;

		double nuF_3_NO_;
		double nuOX_3_NO_;
		double nuW1_4_NO_;
		double nuNO_5f_NO_;
		double nuW2_5f_NO_;
		double nuNO_5b_NO_;
		double nuW2_5b_NO_;
		double nuW3_6_NO_;
		double nuW3_7_NO_;

		double Kc5_NO_;
	};
}

#include "VirtualChemistry.hpp"

#endif /* OpenSMOKEpp_VirtualChemistry */

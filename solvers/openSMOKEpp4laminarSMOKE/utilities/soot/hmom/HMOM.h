/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Benedetta Franzelli, Agnes Livia Bodor, Alberto Cuoci        |
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
|   Copyright(C) 2016 B. Franzelli, A.L. Bodor, A. Cuoci                  |
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
#include "Grammar_HMOM.h"

#ifndef OpenSMOKE_HMOM_H
#define	OpenSMOKE_HMOM_H

namespace OpenSMOKE
{
	//!  A class for modeling soot formation through the Hybrid Method of Moments
	/*!
	A class for modeling soot formation through the Hybrid Method of Moments
	*/
	class HMOM
	{

		enum SootPlanckCoefficient { SOOT_PLANCK_COEFFICIENT_NONE, SOOT_PLANCK_COEFFICIENT_SMOOKE, SOOT_PLANCK_COEFFICIENT_KENT, SOOT_PLANCK_COEFFICIENT_SAZHIN };
	
	public:

		/**
		*@brief Default constructor
		*/
		HMOM();

		/**
		*@brief Default destructor
		*/
		~HMOM();

		/**
		*@brief Setup from a dictionary
		*@param dictionary name
		*/
		void SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary);

		/**
		*@brief Sets the species to be assumed as PAHs
		*@param pah_species PAH species
		*/
		void SetPAH(const std::vector<std::string> pah_species);

		/**
		*@brief Enables soot nucleation
		*@param flag nucleation model (0=off, 1=on)
		*/
		void SetNucleation(const int flag);

		/**
		*@brief Enables soot surface growth
		*@param flag surface growth model (0=off, 1=on)
		*/
		void SetSurfaceGrowth(const int flag);

		/**
		*@brief Enables soot condensation
		*@param flag condensation model (0=off, 1=on)
		*/
		void SetCondensation(const int flag);

		/**
		*@brief Enables soot oxidation
		*@param flag oxidation model (0=off, 1=on)
		*/
		void SetOxidation(const int flag);

		/**
		*@brief Enables soot coagulation
		*@param flag coagulation model (0=off, 1=on)
		*/
		void SetCoagulation(const int flag);

		/**
		*@brief Enables soot coagulation (continous)
		*@param flag continous coagulation model (0=off, 1=on)
		*/
		void SetCoagulationContinous(const int flag);
		
		/**
		*@brief Enables thermophoretic effect
		*@param flag thermophoretic model (0=off, 1=on)
		*/
		void SetThermophoreticModel(const int flag);

		/**
		*@brief Sets the number of carbon atoms in the primary particle
		*@param n_carbon_pah the number of carbon atoms in the primary particle
		*/
		void SetNumberCarbonPAH(const int n_carbon_pah);

		/**
		*@brief Sets the fractal diametr model
		*@param volume_to_surface the fractal diametr model (0 or 1, default equal to 1)
		*/
		void SetFractalDiameterModel(const int volume_to_surface);

		/**
		*@brief Sets the collision diametr model
		*@param volume_to_surface the collision diametr model (1 or 2, default equal to 2)
		*/
		void SetCollisionDiameterModel(const int dc_model);
		
		/**
		*@brief Sets if the PAH consumption rate is turned on or off
		*@param flag if true, consumption of PAH is turned on
		*/
		void SetPAHConsumption(const bool flag);

		/**
		*@brief Sets the normalized moments
		*@param M00_normalized normalized M00 [mol/m3]
		*@param M10_normalized normalized M10 [mol/m3]
		*@param M01_normalized normalized M01 [mol/m3]
		*@param N0_normalized normalized N0 [mol/m3]
		*/
		void SetNormalizedMoments(	const double M00_normalized, const double M10_normalized,
									const double M01_normalized, const double N0_normalized);

		/**
		*@brief Sets temperature and pressure
		*@param T temperature [K]
		*@param P_Pa pressure [Pa]
		*/
		void SetTemperatureAndPressure(const double T, const double P_Pa);

		/**
		*@brief Sets mass fractions of relevant species
		*@param mass_fraction_OH OH mass fraction
		*@param mass_fraction_H H mass fraction
		*/
		void SetMassFractions(const double mass_fraction_OH, const double mass_fraction_H);
		
		/**
		*@brief Sets concentrations of relevant species
		*@param units units of concentrations (kmol/m3 or mol/cm3)
		*@param conc_OH concentration of OH [mol/m3]
		*@param conc_OH concentration of H [mol/m3]
		*@param conc_OH concentration of H2O [mol/m3]
		*@param conc_OH concentration of H2 [mol/m3]
		*@param conc_OH concentration of C2H2 [mol/m3]
		*@param conc_OH concentration of O2 [mol/m3]
		*@param conc_OH concentration of PAH [mol/m3]
		*/
		void SetConcentrations(	const std::string units,
								const double conc_OH, const double conc_H, const double conc_H2O,
								const double conc_H2, const double conc_C2H2, const double conc_O2,
								const double conc_PAH);

		/**
		*@brief Sets the laminar viscosity of the mixture
		*@param viscosity the laminar viscosity [kg/m/s]
		*/
		void SetViscosity(const double viscosity);

		/**
		*@brief Sets the radiative heat transfer from soot
		*@param flag true if the radiative heat transfer from soot is turned on
		*/
		void SetRadiativeHeatTransfer(const bool flag);

		/**
		*@brief Sets the law for calculating the soot Planck mean absorption coefficient
		*@param soot_planck_coefficient the law for calculating the soot Planck mean absorption coefficient
		*/
		void SetPlanckAbsorptionCoefficient(const SootPlanckCoefficient soot_planck_coefficient);

		/**
		*@brief Sets the law for calculating the soot Planck mean absorption coefficient
		*@param soot_planck_coefficient the law for calculating the soot Planck
		mean absorption coefficient: Smooke (default) | Kent | Sazhin
		*/
		void SetPlanckAbsorptionCoefficient(const std::string soot_planck_coefficient);

		/**
		*@brief Sets the Schmidt number for moments
		*@param value Schmidt number for moments
		*/
		void SetSchmidtNumber(const double value);

		/**
		*@brief Calculates the source terms for moment equations
		*/
		void CalculateSourceMoments();

		/**
		*@brief Returns the (normalized) M00 source term in [mol/m3/s]
		*/
		double SourceM00() const { return source_all_(0); }

		/**
		*@brief Returns the (normalized) M10 source term in [mol/m3/s]
		*/
		double SourceM10() const { return source_all_(1); }

		/**
		*@brief Returns the (normalized) M01 source term in [mol/m3/s]
		*/
		double SourceM01() const { return source_all_(2); }

		/**
		*@brief Returns the (normalized) N0 source term in [mol/m3/s]
		*/
		double SourceN0() const { return source_all_(3); }

		/**
		*@brief Returns the soot volume fraction [-]
		*/
		double SootVolumeFraction() const;

		/**
		*@brief Returns the soot particle diameter [m]
		*/
		double SootParticleDiameter() const;

		/**
		*@brief Returns the soot collisional particle diameter [m]
		*/
		double SootCollisionParticleDiameter() const;

		/**
		*@brief Returns the number of primary particles [-]
		*/
		double SootNumberOfPrimaryParticles() const;

		/**
		*@brief Returns the soot particle number density [#/m3]
		*/
		double SootParticleNumberDensity() const;

		/**
		*@brief Returns nucleated particle volume [m3]
		*/
		double V0() const { return V0_;  }

		/**
		*@brief Returns nucleated particle surface [m2]
		*/
		double S0() const { return S0_; }

		/**
		*@brief Returns true if the HMOM is turned on
		*/
		bool is_active() const { return is_active_;  }

		/**
		*@brief Returns the species to be considered as PAH
		*/
		const std::vector<std::string>& pah_species() const { return pah_species_; }

		/**
		*@brief Returns the number of moments
		*/
		unsigned int n_moments() const { return n_moments_; }

		/**
		*@brief Returns the dimerization rate [mol/m3/s]
		*/
		double dimerization_rate() const { return dimerization_rate_; }

		/**
		*@brief Returns the nucleation model
		*/
		int nucleation_model() const { return nucleation_model_; }

		/**
		*@brief Returns the surface growth model
		*/
		int surface_growth_model() const { return surface_growth_model_; }

		/**
		*@brief Returns the oxidation model
		*/
		int oxidation_model() const { return oxidation_model_; }

		/**
		*@brief Returns the condensation model
		*/
		int condensation_model() const { return condensation_model_; }

		/**
		*@brief Returns the coagulation model
		*/
		int coagulation_model() const { return coagulation_model_; }

		/**
		*@brief Returns the continous coagulation model
		*/
		int continous_coagulation_model() const { return coagulation_continous_model_; }

		/**
		*@brief Returns the thermophoretic model
		*/
		int thermophoretic_model() const { return thermophoretic_model_; }

		/**
		*@brief Returns the (normalized) source terms for moment equations calculated internally in [mol/m3/s]
		*/
		const Eigen::VectorXd& sources() const { return source_all_; }

		/**
		*@brief Returns the (normalized) nucleation source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_nucleation() const { return source_nucleation_; }

		/**
		*@brief Returns the (normalized) surface growth source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_growth() const { return source_growth_; }

		/**
		*@brief Returns the (normalized) oxidation source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_oxidation() const { return source_oxidation_; }

		/**
		*@brief Returns the (normalized) condensation source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_condensation() const { return source_condensation_; }

		/**
		*@brief Returns the (normalized) coagulation (overall) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_overall() const { return source_coagulation_all_; }

		/**
		*@brief Returns the (normalized) coagulation (discrete) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_discrete() const { return source_coagulation_; }

		/**
		*@brief Returns the (normalized) coagulation (discrete, small-small) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_discrete_ss() const { return source_coagulation_ss_; }

		/**
		*@brief Returns the (normalized) coagulation (discrete, small-large) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_discrete_sl() const { return source_coagulation_sl_; }

		/**
		*@brief Returns the (normalized) coagulation (discrete, large-large) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_discrete_ll() const { return source_coagulation_ll_; }

		/**
		*@brief Returns the (normalized) coagulation (continous) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_continous() const { return source_coagulation_continous_; }

		/**
		*@brief Returns the (normalized) coagulation (continous, small-small) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_continous_ss() const { return source_coagulation_continous_ss_; }

		/**
		*@brief Returns the (normalized) coagulation (continous, small-large) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_continous_sl() const { return source_coagulation_continous_sl_; }

		/**
		*@brief Returns the (normalized) coagulation (continous, large-large) source terms for moment equations [mol/m3/s]
		*/
		const Eigen::VectorXd& sources_coagulation_continous_ll() const { return source_coagulation_continous_ll_; }

		/**
		*@brief Returns if PAH consumption is turned on or off
		*/
		bool PAHConsumption() const { return pah_consumption_; }
		
		/**
		*@brief Returns the PAH consumption rate [mol/m3/s]
		*/
		double PAHConsumptionRate() const;

		/**
		*@brief Returns the soot Planck mean absorption coefficient
		*@param T temperature of gaseous mixture [K]
		*@param fv soot volume fraction
		*/
		double planck_coefficient(const double T, const double fv) const;

		/**
		*@brief Returns true is soot has to be accounted for in radiative heat transfer
		*/
		bool radiative_heat_transfer() const { return radiative_heat_transfer_; }

		/**
		*@brief Returns the Schmidt number for moments
		*/
		double schmidt_number() const { return schmidt_number_; }

	private:

		/**
		*@brief Allocates internal memory
		*/
		void MemoryAllocation();

		/**
		*@brief Returns the moment (i,j)
		*@param i index
		*@param j index
		*/
		double GetMoment(const double i, const double j) const;

		/**
		*@brief Returns the moment (i,j)
		*@param i index
		*@param j index
		*/
		double GetMissingMoment(const double i, const double j) const;

		/**
		*@brief Reconstruct moments (internally)
		*/
		void GetMoments();

		/**
		*@brief Calculates the BetaC coefficient
		*/
		double GetBetaC();

		/**
		*@brief Calculates the dimer concentration
		*/
		void DimerConcentration();

		/**
		*@brief Calculates the nucleation source terms
		*/
		void SootNucleationM4();

		/**
		*@brief Calculates the kinetic constants for surface growth and oxidation
		*/
		void SootKineticConstants();

		/**
		*@brief Calculates the surface growth source terms
		*/
		void SootSurfaceGrowthM4();

		/**
		*@brief Calculates the oxidation source terms
		*/
		void SootOxidationM4();

		/**
		*@brief Calculates the condensation source terms
		*/
		void SootCondensationM4();

		/**
		*@brief Calculates the coagulation source terms
		*/
		void SootCoagulationM4();

		/**
		*@brief Calculates the coagulation source terms (small/small)
		*/
		void SootCoagulationSmallSmallM4();

		/**
		*@brief Calculates the coagulation source terms (small/large)
		*/
		void SootCoagulationSmallLargeM4();

		/**
		*@brief Calculates the coagulation source terms (large/large)
		*/
		void SootCoagulationLargeLargeM4();

		/**
		*@brief Calculates the coagulation (continous) source terms
		*/
		void SootCoagulationContinousM4();

		/**
		*@brief Calculates the coagulation (continous) source terms (small/small)
		*/
		void SootCoagulationContinousSmallSmallM4(const double lambda);

		/**
		*@brief Calculates the coagulation (continous) source terms (small/large)
		*/
		void SootCoagulationContinousSmallLargeM4(const double lambda);

		/**
		*@brief Calculates the coagulation (continous) source terms (large/large)
		*/
		void SootCoagulationContinousLargeLargeM4(const double lambda);

	private:

		double M00_normalized_;		//!< normalized moment [mol/m3]
		double M10_normalized_;		//!< normalized moment [mol/m3]
		double M01_normalized_;		//!< normalized moment [mol/m3]
		double N0_normalized_;		//!< normalized moment [mol/m3]

		double M00_;				//!< moment M00 [#/m3]
		double M10_;				//!< moment M10 [#]
		double M01_;				//!< moment M01 [#/m]
		double N0_;					//!< moment N0 [#/m3]
			
		double NL_;					//!< number density associated to large particles (2nd mode) [#/m3]
		double NLVL_;				//!< mean volume of large particles (2nd mode) multiplied by NL [#]
		double NLSL_;				//!< mean surface of large particles (2nd mode) multiplied by NL [#/m]

		double T_;					//!< temperature [K]
		double P_Pa_;				//!< pressure [Pa]

		double conc_OH_;			//!< concentration [mol/cm3]
		double conc_H_;				//!< concentration [mol/cm3]
		double conc_H2O_;			//!< concentration [mol/cm3]
		double conc_H2_;			//!< concentration [mol/cm3]
		double conc_C2H2_;			//!< concentration [mol/cm3]
		double conc_O2_;			//!< concentration [mol/cm3]
		double conc_PAH_;			//!< concentration [mol/cm3]

		double mass_fraction_H_;	//!< mass fraction
		double mass_fraction_OH_;	//!< mass fraction

		double viscosity_;			//!< laminar viscosity [kg/m/s]

		double conc_DIMER_;			//!< concentration of dimers [???]
		double kox_;				//!< oxidation kinetic constant [???]
		double ksg_;				//!< surface growth kinetic constant [???]
		double betaN_;				//!< [???]

		Eigen::VectorXd source_nucleation_;						//!< nucleation source terms [mol/m3/s]
		Eigen::VectorXd source_growth_;							//!< surface growth source terms [mol/m3/s]
		Eigen::VectorXd source_oxidation_;						//!< oxidation source terms [mol/m3/s]
		Eigen::VectorXd source_condensation_;					//!< condensation source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_;					//!< coagulation source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_ss_;					//!< coagulation (small/small) source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_ll_;					//!< coagulation (small/large) source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_sl_;					//!< coagulation (large/large) source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_continous_;			//!< coagulation continous source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_continous_ss_;		//!< coagulation continous (small/small) source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_continous_ll_;		//!< coagulation continous (small/large) source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_continous_sl_;		//!< coagulation continous (large/large) source terms [mol/m3/s]
		Eigen::VectorXd source_coagulation_all_;				//!< coagulation (overall) source terms [mol/m3/s]
		Eigen::VectorXd source_all_;							//!< overall source terms [mol/m3/s]


	private:

		bool is_active_;						//!< true if HMOM is turned on

		unsigned int n_moments_;				//!< number of moments			
		int nucleation_model_;					//!< nucleation model
		int condensation_model_;				//!< condensation model
		int surface_growth_model_;				//!< surface growth model
		int oxidation_model_;					//!< oxidation model
		int coagulation_model_;					//!< coagulation model
		int coagulation_continous_model_;		//!< coagulation (continous) model
		int thermophoretic_model_;				//!< thermophoretic model

		double dimerization_rate_;				//!< dimerization rate [mol/m3/s]
		double pah_volume_;						//!< PAH volume [m3]
		double dimer_volume_;					//!< dimer volume [m3]
		double dimer_surface_;					//!< dimer surface [m2]
		double V0_;								//!< nucleated particle volume [m3]
		double S0_;								//!< nucleated particle surface [m2]
		double VC2_;							//!< [???]

		double Av_fractal_;						//!< [???]
		double As_fractal_;						//!< [???]
		double K_fractal_;						//!< [???]

		double D_collisional_;					//!< [???]
		double Av_collisional_;					//!< [???]
		double As_collisional_;					//!< [???]
		double K_collisional_;					//!< [???]

		std::vector<std::string> pah_species_;	//!< names of PAH species
		bool pah_consumption_;					//!< PAH consumption

		SootPlanckCoefficient soot_planck_coefficient_;
		bool radiative_heat_transfer_;

		double schmidt_number_;

	private:

		static const double pi;
		static const double WC;
		static const double AvogadroNumber;
		static const double K_diam;
		static const double K_spher;
		static const double rho_soot;
		static const double Cfm;
		static const double betaN_TV;
		static const double Rgas;
		static const double KB;
		static const double surface_density;

		static const double A1f;
		static const double n1f;
		static const double E1f;
		static const double A1b;
		static const double n1b;
		static const double E1b;
		static const double A2f;
		static const double n2f;
		static const double E2f;
		static const double A2b;
		static const double n2b;
		static const double E2b;
		static const double A3f;
		static const double n3f;
		static const double E3f;
		static const double A3b;
		static const double n3b;
		static const double E3b;
		static const double A4;
		static const double n4;
		static const double E4;
		static const double A5;
		static const double n5;
		static const double E5;
		static const double eff6;
	};
}

#include "HMOM.hpp"

#endif	/* OpenSMOKE_HMOM_H */
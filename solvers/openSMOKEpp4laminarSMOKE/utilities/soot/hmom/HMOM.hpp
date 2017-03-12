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

namespace OpenSMOKE
{
	const double HMOM::surface_density = 1.7e19;				// [#/m2]
	const double HMOM::pi = 3.141592653589793;					// [-]
	const double HMOM::WC = 12.e-3;								// [kg/mol]
	const double HMOM::AvogadroNumber = 6.022e23;				// [#/mol]
	const double HMOM::K_diam = std::pow(6. / pi, 1. / 3.);			// [-]
	const double HMOM::K_spher = std::pow(36.*pi, 1. / 3.);		// [-]
	const double HMOM::rho_soot = 1800.;						// [kg/m3]
	const double HMOM::Rgas = 8.3143;							// [J/mol/K]
	const double HMOM::KB = Rgas / AvogadroNumber;				// [J/K]
	const double HMOM::Cfm = std::sqrt(pi*KB / 2.0 / rho_soot);
	const double HMOM::betaN_TV = 2.2 * 4 * std::sqrt(2.0)*std::pow(K_diam, 2.)*Cfm;

	// Kinetics
	const double HMOM::A1f = 6.72e1;
	const double HMOM::n1f = 3.33;
	const double HMOM::E1f = 6.09  * 1e3 / Rgas;
	const double HMOM::A1b = 6.44e-1;
	const double HMOM::n1b = 3.79;
	const double HMOM::E1b = 27.96  * 1e3 / Rgas;
	const double HMOM::A2f = 1.00e8;
	const double HMOM::n2f = 1.80;
	const double HMOM::E2f = 68.42  * 1e3 / Rgas;
	const double HMOM::A2b = 8.68e4;
	const double HMOM::n2b = 2.36;
	const double HMOM::E2b = 25.46  * 1e3 / Rgas;
	const double HMOM::A3f = 1.13e16;
	const double HMOM::n3f = -0.06;
	const double HMOM::E3f = 476.05  * 1e3 / Rgas;
	const double HMOM::A3b = 4.17e13;
	const double HMOM::n3b = 0.15;
	const double HMOM::E3b = 0.00  * 1e3 / Rgas;
	const double HMOM::A4 = 2.52e9;
	const double HMOM::n4 = 1.10;
	const double HMOM::E4 = 17.13  * 1e3 / Rgas;
	const double HMOM::A5 = 2.20e12;
	const double HMOM::n5 = 0.00;
	const double HMOM::E5 = 31.38  * 1e3 / Rgas;
	const double HMOM::eff6 = 0.13;

	HMOM::HMOM()
	{
		// HMOM is turned on
		is_active_ = true;

		// Total number of moments
		n_moments_ = 4;
		
		// PAH species
		pah_species_.resize(1);
		pah_species_[0] = "C10H8";

		// Physical/Chemical phenomena
		nucleation_model_ = 1;
		condensation_model_ = 1;
		surface_growth_model_ = 1;
		oxidation_model_ = 1;
		coagulation_model_ = 1;
		coagulation_continous_model_ = 1;
		thermophoretic_model_ = 1;
		schmidt_number_ = 50.;
		radiative_heat_transfer_ = true;
		soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_SMOOKE;

		// Nuclei particles
		const int n_carbon_pah = 16;
		SetNumberCarbonPAH(n_carbon_pah);

		// Fractal dimension model
		const int volume_to_surface = 1;
		SetFractalDiameterModel(volume_to_surface);

		// Collision diameter model
		const int dc_model = 2;
		SetCollisionDiameterModel(dc_model);
		
		// PAH consumption
		SetPAHConsumption(true);

		// Memory allocation
		MemoryAllocation();
	}

	void HMOM::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_HMOM grammar;
		dictionary.SetGrammar(grammar);

		if (dictionary.CheckOption("@HMOM") == true)
			dictionary.ReadBool("@HMOM", is_active_);

		if (dictionary.CheckOption("@NumberOfCarbonPAH") == true)
		{
			int flag;
			dictionary.ReadInt("@NumberOfCarbonPAH", flag);
			SetNumberCarbonPAH(flag);
		}

		if (dictionary.CheckOption("@FractalDiameterModel") == true)
		{
			int flag;
			dictionary.ReadInt("@FractalDiameterModel", flag);
			SetFractalDiameterModel(flag);
		}

		if (dictionary.CheckOption("@CollisionDiameterModel") == true)
		{
			int flag;
			dictionary.ReadInt("@CollisionDiameterModel", flag);
			SetCollisionDiameterModel(flag);
		}
		
		if (dictionary.CheckOption("@PAHConsumption") == true)
		{
			bool flag;
			dictionary.ReadBool("@PAHConsumption", flag);
			SetPAHConsumption(flag);
		}

		if (dictionary.CheckOption("@PAH") == true)
		{
			dictionary.ReadOption("@PAH", pah_species_);
			SetPAH(pah_species_);
		}

		if (dictionary.CheckOption("@NucleationModel") == true)
		{
			int flag;
			dictionary.ReadInt("@NucleationModel", flag);
			SetNucleation(flag);
		}

		if (dictionary.CheckOption("@SurfaceGrowthModel") == true)
		{
			int flag;
			dictionary.ReadInt("@SurfaceGrowthModel", flag);
			SetSurfaceGrowth(flag);
		}

		if (dictionary.CheckOption("@OxidationModel") == true)
		{
			int flag;
			dictionary.ReadInt("@OxidationModel", flag);
			SetOxidation(flag);
		}

		if (dictionary.CheckOption("@CondensationModel") == true)
		{
			int flag;
			dictionary.ReadInt("@CondensationModel", flag);
			SetCondensation(flag);
		}

		if (dictionary.CheckOption("@CoagulationModel") == true)
		{
			int flag;
			dictionary.ReadInt("@CoagulationModel", flag);
			SetCoagulation(flag);
		}

		if (dictionary.CheckOption("@ContinousCoagulationModel") == true)
		{
			int flag;
			dictionary.ReadInt("@ContinousCoagulationModel", flag);
			SetCoagulationContinous(flag);
		}

		if (dictionary.CheckOption("@ThermophoreticModel") == true)
		{
			int flag;
			dictionary.ReadInt("@ThermophoreticModel", flag);
			SetThermophoreticModel(flag);
		}

		if (dictionary.CheckOption("@RadiativeHeatTransfer") == true)
			dictionary.ReadBool("@RadiativeHeatTransfer", radiative_heat_transfer_);

		if (dictionary.CheckOption("@PlanckCoefficient") == true)
		{
			std::string flag;
			dictionary.ReadString("@PlanckCoefficient", flag);
			SetPlanckAbsorptionCoefficient(flag);
		}

		if (dictionary.CheckOption("@SchmidtNumber") == true)
		{
			double value;
			dictionary.ReadDouble("@SchmidtNumber", value);
			SetSchmidtNumber(value);
		}
	}

	void HMOM::MemoryAllocation()
	{
		source_nucleation_.resize(n_moments_);
		source_nucleation_.setZero();
		source_growth_.resize(n_moments_);
		source_growth_.setZero();
		source_oxidation_.resize(n_moments_);
		source_oxidation_.setZero();
		source_condensation_.resize(n_moments_);
		source_condensation_.setZero();
		source_coagulation_.resize(n_moments_);
		source_coagulation_.setZero();
		source_coagulation_ss_.resize(n_moments_);
		source_coagulation_ss_.setZero();
		source_coagulation_ll_.resize(n_moments_);
		source_coagulation_ll_.setZero();
		source_coagulation_sl_.resize(n_moments_);
		source_coagulation_sl_.setZero();
		source_coagulation_continous_.resize(n_moments_);
		source_coagulation_continous_.setZero();
		source_coagulation_continous_ss_.resize(n_moments_);
		source_coagulation_continous_ss_.setZero();
		source_coagulation_continous_ll_.resize(n_moments_);
		source_coagulation_continous_ll_.setZero();
		source_coagulation_continous_sl_.resize(n_moments_);
		source_coagulation_continous_sl_.setZero();
		source_coagulation_all_.resize(n_moments_);
		source_coagulation_all_.setZero();
		source_all_.resize(n_moments_);
		source_all_.setZero();

		NL_ = 0.;		// [#/m3]
		NLVL_ = 0.;		// [m3]
		NLSL_ = 0.;		// [#/m]
	}

	void HMOM::SetPAH(const std::vector<std::string> pah_species)
	{
		pah_species_ = pah_species;
	}

	void HMOM::SetNucleation(const int flag)
	{
		nucleation_model_ = flag;
	}

	void HMOM::SetSurfaceGrowth(const int flag)
	{
		surface_growth_model_ = flag;
	}

	void HMOM::SetCondensation(const int flag)
	{
		condensation_model_ = flag;
	}

	void HMOM::SetOxidation(const int flag)
	{
		oxidation_model_ = flag;
	}

	void HMOM::SetCoagulation(const int flag)
	{
		coagulation_model_ = flag;
	}

	void HMOM::SetCoagulationContinous(const int flag)
	{
		coagulation_continous_model_ = flag;
	}

	void HMOM::SetThermophoreticModel(const int flag)
	{
		thermophoretic_model_ = flag;
	}

	void HMOM::SetNumberCarbonPAH(const int n_carbon_pah)
	{
		pah_volume_ = (WC / rho_soot / AvogadroNumber) * double(n_carbon_pah);	// [m3]
		dimer_volume_ = 2. * pah_volume_;										// [m3]
		dimer_surface_ = K_spher * std::pow(dimer_volume_, 2. / 3.);			// [m2]
		V0_ = 2.* dimer_volume_;												// nucleated particle volume (eq. 16) [m3]
		S0_ = K_spher * std::pow(V0_, 2. / 3.);									// nucleated particle surface (eq. 17) [m2]
		VC2_ = (WC / rho_soot / AvogadroNumber)*2.;								// [m3]
	}

	void HMOM::SetCollisionDiameterModel(const int dc_model)
	{
		if (dc_model == 1)
		{
			D_collisional_ = 2.0;
			Av_collisional_ = 1.0 / 3.0;
			As_collisional_ = 0.0;
			K_collisional_ = K_diam;
		}
		else if (dc_model == 2)
		{
			D_collisional_ = 1.8;
			Av_collisional_ = 1.0 - 2. / D_collisional_;
			As_collisional_ = 3.0 / D_collisional_ - 1.;
			K_collisional_ = 6.0 / std::pow(36.*pi, 1. / D_collisional_);
		}
	}

	void HMOM::SetFractalDiameterModel(const int volume_to_surface)
	{
		double chi_fractal = 0.0;

		if (volume_to_surface == 1)	chi_fractal = -0.2043;
		else						chi_fractal = 0.0;

		Av_fractal_ = -2.0*chi_fractal - 1.;
		As_fractal_ = 3.0*chi_fractal;
		K_fractal_ = 2.0 / 3.0 * std::pow(1.0 / 36.0 / pi, chi_fractal);
	}
	
	void HMOM::SetPAHConsumption(const bool flag)
	{
		pah_consumption_ = flag;
	}

	void HMOM::SetPlanckAbsorptionCoefficient(const SootPlanckCoefficient soot_planck_coefficient)
	{
		soot_planck_coefficient_ = soot_planck_coefficient;
	}

	void HMOM::SetPlanckAbsorptionCoefficient(const std::string label)
	{
		if (label == "Smooke")
			soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_SMOOKE;
		else if (label == "Kent")
			soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_KENT;
		else if (label == "Sazhin")
			soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_SAZHIN;
		else if (label == "none")
			soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_NONE;
		else
			OpenSMOKE::FatalErrorMessage("@PlanckCoefficient: available options: none | Smooke (default) | Kent | Sazhin");
	}

	void HMOM::SetRadiativeHeatTransfer(const bool flag)
	{
		radiative_heat_transfer_ = flag;
	}

	void HMOM::SetSchmidtNumber(const double value)
	{
		schmidt_number_ = value;
	}
	
	void HMOM::SetNormalizedMoments(	const double M00_normalized, const double M10_normalized,
										const double M01_normalized, const double N0_normalized)
	{
		// Set normalized moments
		M00_normalized_ = M00_normalized;		// [mol/m3]
		M10_normalized_ = M10_normalized;		// [mol/m3]
		M01_normalized_ = M01_normalized;		// [mol/m3]
		N0_normalized_  = N0_normalized;		// [mol/m3]

		// Reconstruct moments
		GetMoments();
	}

	void HMOM::SetTemperatureAndPressure(const double T, const double P_Pa)
	{
		T_ = T;				// [K]
		P_Pa_ = P_Pa;		// [Pa]
	}

	void HMOM::SetMassFractions(const double mass_fraction_OH, const double mass_fraction_H)
	{
		mass_fraction_OH_ = mass_fraction_OH;	// [-]
		mass_fraction_H_ = mass_fraction_H;		// [-]
	}

	void HMOM::SetConcentrations(	const std::string units, const double conc_OH, const double conc_H, const double conc_H2O,
								 const double conc_H2, const double conc_C2H2, const double conc_O2,
								 const double conc_PAH)
	{
		if (units == "mol/cm3")
		{
			conc_OH_ = conc_OH;		// [mol/cm3]
			conc_H_ = conc_H;		// [mol/cm3]
			conc_H2O_ = conc_H2O;		// [mol/cm3]
			conc_H2_ = conc_H2;		// [mol/cm3]
			conc_C2H2_ = conc_C2H2;		// [mol/cm3]
			conc_O2_ = conc_O2;		// [mol/cm3]
			conc_PAH_ = conc_PAH;		// [mol/cm3]
		}
		else if (units == "kmol/m3")
		{
			conc_OH_ = conc_OH / 1.e3;		// [mol/cm3]
			conc_H_ = conc_H / 1.e3;		// [mol/cm3]
			conc_H2O_ = conc_H2O / 1.e3;		// [mol/cm3]
			conc_H2_ = conc_H2 / 1.e3;		// [mol/cm3]
			conc_C2H2_ = conc_C2H2 / 1.e3;		// [mol/cm3]
			conc_O2_ = conc_O2 / 1.e3;		// [mol/cm3]
			conc_PAH_ = conc_PAH / 1.e3;		// [mol/cm3]
		}
	}

	void HMOM::SetViscosity(const double viscosity)
	{
		viscosity_ = viscosity;		// [kg/m/s]
	}

	void HMOM::CalculateSourceMoments()
	{
		// Dimer concentration
		DimerConcentration();

		// Source terms: nucleation
		SootNucleationM4();

		// Source terms: surface growth and oxidation
		if ((surface_growth_model_ > 0) || (oxidation_model_ > 0))
		{
			SootKineticConstants();

			if (surface_growth_model_ > 0)
				SootSurfaceGrowthM4();

			if (oxidation_model_ > 0)
				SootOxidationM4();
		}

		// Source terms: condensation
		if (condensation_model_ > 0)
			SootCondensationM4();

		// Source terms: coagulation
		if (coagulation_model_ > 0)
			SootCoagulationM4();

		// Source terms: coagulation (continous)
		if (coagulation_continous_model_ > 0)
			SootCoagulationContinousM4();

		// Source terms: overall
		for (unsigned int i = 0; i < n_moments_; i++)
		{
			if (std::isnan(source_nucleation_(i)) == true)
				source_nucleation_(i) = 0.;
			if (std::isnan(source_oxidation_(i)) == true)
				source_oxidation_(i) = 0.;
			if (std::isnan(source_condensation_(i)) == true)
				source_condensation_(i) = 0.;
			if (std::isnan(source_growth_(i)) == true)
				source_growth_(i) = 0.;

			if (std::isnan(source_coagulation_ss_(i)) == true)
				source_coagulation_ss_(i) = 0.;
			if (std::isnan(source_coagulation_sl_(i)) == true)
				source_coagulation_sl_(i) = 0.;
			if (std::isnan(source_coagulation_ll_(i)) == true)
				source_coagulation_ll_(i) = 0.;

			if (std::isnan(source_coagulation_continous_ss_(i)) == true)
				source_coagulation_continous_ss_(i) = 0.;
			if (std::isnan(source_coagulation_continous_sl_(i)) == true)
				source_coagulation_continous_sl_(i) = 0.;
			if (std::isnan(source_coagulation_continous_ll_(i)) == true)
				source_coagulation_continous_ll_(i) = 0.;

			source_coagulation_(i) =	source_coagulation_ss_(i) +
										source_coagulation_sl_(i) +
										source_coagulation_ll_(i);

			source_coagulation_continous_(i) =	source_coagulation_continous_ss_(i) +
												source_coagulation_continous_sl_(i) +
												source_coagulation_continous_ll_(i);

			if (source_coagulation_(i) == 0.)
			{
				source_coagulation_all_(i) = source_coagulation_continous_(i);
			}
			else
			{
				if (source_coagulation_continous_(i) == 0.)
				{
					source_coagulation_all_(i) = source_coagulation_(i);
				}
				else
				{
					source_coagulation_all_(i) = source_coagulation_continous_(i)*source_coagulation_(i) /
												(source_coagulation_continous_(i) + source_coagulation_(i));
				}
			}

			source_all_(i) = source_nucleation_(i) + source_growth_(i) + source_oxidation_(i) +
				source_coagulation_all_(i) + source_condensation_(i);
		}
	}

	void HMOM::DimerConcentration()
	{
		betaN_ = betaN_TV * std::sqrt(T_) * std::pow(dimer_volume_, 1./6.);

		const double sticking_coefficient = 2e-3;
		const double betaC = GetBetaC();
		const double KfmPAH = betaN_TV * std::sqrt(T_) * std::pow(pah_volume_, 1. / 6.);
		const double aromatics_conc = conc_PAH_*1.e6;
		dimerization_rate_ = 0.5* KfmPAH * std::pow(aromatics_conc, 2.) * sticking_coefficient;	// [mol/m3/s]
		const double delta = std::pow(betaC, 2.) + 4.0 * betaN_ * dimerization_rate_ * std::pow(AvogadroNumber, 2.0);

		conc_DIMER_ = (std::sqrt(delta) - betaC) / (2.0*betaN_*AvogadroNumber);
	}

	double HMOM::PAHConsumptionRate() const
	{
		return 2.*dimerization_rate_*AvogadroNumber;		// [mol/m3/s]
	}

	void HMOM::SootKineticConstants()
	{
		const double k1f = A1f * std::pow(T_, n1f) * std::exp(-E1f / T_);
		const double k2f = A2f * std::pow(T_, n2f) * std::exp(-E2f / T_);
		const double k3f = A3f * std::pow(T_, n3f) * std::exp(-E3f / T_);
		const double k1b = A1b * std::pow(T_, n1b) * std::exp(-E1b / T_);
		const double k2b = A2b * std::pow(T_, n2b) * std::exp(-E2b / T_);
		const double k3b = A3b * std::pow(T_, n3b) * std::exp(-E3b / T_);
		const double k4 = A4  * std::pow(T_, n4)  * std::exp(-E4 / T_);

		const double ratio = k1b*conc_H2O_ + k2b*conc_H2_ + k3b*conc_H_ + k4*conc_C2H2_;
		double conc_sootStar = 0.;
		
		if (ratio > 0.)
			conc_sootStar = (k1f*conc_OH_ + k2f*conc_H_ + k3f) / ratio;

		if (mass_fraction_H_ < 2.e-9)
		{
			if (mass_fraction_OH_ < 2.e-8)
			{
				conc_sootStar = 0.0;
			}
		}
		conc_sootStar = conc_sootStar / (1.0 + conc_sootStar);
		conc_sootStar = std::max(conc_sootStar, 0.0);

		ksg_ = A4 * std::pow(T_, n4) * std::exp(-E4 / T_) * conc_C2H2_ * conc_sootStar;

		const double k6 = 8.94  * eff6 * std::sqrt(T_) * AvogadroNumber;
		kox_ = A5 * std::pow(T_, n5) * std::exp(-E5 / T_) * conc_O2_ * conc_sootStar +
			(0.5 / surface_density) * k6 * conc_OH_;
	}

	void HMOM::SootSurfaceGrowthM4()
	{
		source_growth_(0) = 0.;
		source_growth_(1) = ksg_ * surface_density * VC2_ * GetMoment(0., 1.) / AvogadroNumber / V0_;
		source_growth_(2) = ksg_ * surface_density * VC2_ * K_fractal_ * GetMoment(Av_fractal_, As_fractal_ + 2.) / AvogadroNumber / S0_;
		source_growth_(3) = -ksg_ * surface_density * S0_ * N0_ / AvogadroNumber;
	}

	void HMOM::SootOxidationM4()
	{
		const double coefficient = kox_ * surface_density * VC2_;
		const double M01 = GetMissingMoment(0., 1.);
		const double M_12 = GetMissingMoment(-1., 2.);

		source_oxidation_(0) = -coefficient * N0_ * S0_ / V0_ / AvogadroNumber;								// M00
		source_oxidation_(1) = -coefficient / V0_ / AvogadroNumber * M01;								// M10
		source_oxidation_(2) = -coefficient / S0_ / AvogadroNumber * (std::pow(S0_, 2.)*N0_ / V0_ + 2./3.*( M_12 - S0_*N0_));		// M01
		source_oxidation_(3) = -coefficient * N0_ * S0_ / V0_ / AvogadroNumber;								// N0
	}

	void HMOM::SootCondensationM4()
	{
		const double D_DIM = K_diam * std::pow(dimer_volume_, 1.0 / 3.0);		// [m]
		const double D_NUCL = K_diam * std::pow(V0_, 1.0 / 3.0);	// [m]

		source_condensation_(0) = 0.0;

		source_condensation_(1) = std::pow(D_DIM, 2.) * (GetMoment(0., 0.) + 0.5*dimer_volume_*GetMoment(-1., 0.)) +
			std::pow(K_collisional_, 2.0) * (GetMoment(2.0*Av_collisional_, 2.0*As_collisional_) + 0.5*dimer_volume_*GetMoment(2.0*Av_collisional_ - 1.0, 2.0*As_collisional_)) +
			2.0*D_DIM*K_collisional_ * (GetMoment(Av_collisional_, As_collisional_) + 0.5*dimer_volume_*GetMoment(Av_collisional_ - 1.0, As_collisional_));

		source_condensation_(1) = Cfm * std::sqrt(dimer_volume_) * conc_DIMER_ * std::sqrt(T_) * source_condensation_(1) / V0_;

		source_condensation_(2) = std::pow(D_DIM, 2.) *   (GetMoment(Av_fractal_, As_fractal_ + 1.0) + 0.5*dimer_volume_*GetMoment(Av_fractal_ - 1., As_fractal_ + 1.)) +
			std::pow(K_collisional_, 2.0) *   (GetMoment(Av_fractal_ + 2.0*Av_collisional_, As_fractal_ + 1. + 2.*As_collisional_)
				+ 0.5*dimer_volume_*GetMoment(Av_fractal_ + 2.0*Av_collisional_ - 1.0, As_fractal_ + 1. + 2.0*As_collisional_)) +
			2.0*D_DIM*K_collisional_ * (GetMoment(Av_fractal_ + Av_collisional_, As_fractal_ + 1. + As_collisional_)
				+ 0.5*dimer_volume_*GetMoment(Av_fractal_ + Av_collisional_ - 1.0, As_fractal_ + 1. + As_collisional_));
		source_condensation_(2) = Cfm * std::sqrt(dimer_volume_) * conc_DIMER_ * std::sqrt(T_) * source_condensation_(2) / S0_ * K_fractal_;

		source_condensation_(3) = -Cfm * 1. / std::sqrt(dimer_volume_) * conc_DIMER_ * std::sqrt(T_) * (1.0 + 0.5*dimer_volume_ / V0_) * std::pow(D_DIM + D_NUCL, 2.) * N0_;
	}

	double HMOM::GetBetaC()
	{
		double betaC = std::pow(K_collisional_, 2.0) * std::pow(dimer_volume_, -3. / 6.) * GetMoment(2.0*Av_collisional_, 2.0*As_collisional_) +
			2.0*K_diam*K_collisional_ * std::pow(dimer_volume_, -1. / 6.) * GetMoment(Av_collisional_, As_collisional_) +
			std::pow(K_diam, 2.0) * std::pow(dimer_volume_, 1. / 6.)  * GetMoment(0., 0.) +
			0.5*std::pow(K_collisional_, 2.) * std::pow(dimer_volume_, 3. / 6.)  * GetMoment(2.0*Av_collisional_ - 1.0, 2.0*As_collisional_) +
			2.0*K_diam*K_collisional_ * std::pow(dimer_volume_, 5. / 6.)  * GetMoment(Av_collisional_ - 1.0, As_collisional_) +
			0.5*std::pow(K_diam, 2.) * std::pow(dimer_volume_, 7. / 6.)  * GetMoment(-1., 0.);

		betaC = Cfm*std::sqrt(T_)*betaC;

		return betaC;
	}

	void HMOM::SootNucleationM4()
	{
		// Nucleation source term (equation 15)
		// Mueller et al., Comb. & Flame 156, pp.1143-1155 (2009)

		double nucleation_N0 = 0.;
		if (conc_DIMER_ >= 1.e-7)
			nucleation_N0 = 0.5 * betaN_ * std::pow(conc_DIMER_, 2.) * AvogadroNumber;

		for (unsigned int i = 0; i < n_moments_; i++)
			source_nucleation_(i) = nucleation_N0;
	}

	void HMOM::SootCoagulationM4()
	{
		SootCoagulationSmallSmallM4();
		SootCoagulationSmallLargeM4();
		SootCoagulationLargeLargeM4();
	}

	void HMOM::SootCoagulationSmallSmallM4()
	{
		const double DcNUCL = K_diam * std::pow(V0_, 1. / 3.);
		const double S00 = K_spher * std::pow(2.0*V0_, 2. / 3.);
		const double beta00 = 2.20*Cfm*std::sqrt(2.0 / V0_)*std::pow(2.*DcNUCL, 2.) * std::sqrt(T_);

		source_coagulation_ss_(0) = -0.5* beta00 * std::pow(N0_, 2.) / AvogadroNumber;
		source_coagulation_ss_(1) = 0.0;
		source_coagulation_ss_(2) = 0.5* beta00 * S00 / S0_ * std::pow(N0_, 2.) / AvogadroNumber + 2.0* source_coagulation_ss_(0);
		source_coagulation_ss_(3) = 2.0*source_coagulation_ss_(0);
	}

	void HMOM::SootCoagulationSmallLargeM4()
	{
		const double DcNUCL = K_diam * std::pow(V0_, 1. / 3.);

		{
			const double psi0 = 	std::pow(V0_, -1. / 2.)*(std::pow(K_collisional_, 2.) * GetMissingMoment(2.0*Av_collisional_ - 0.5, 2.0*As_collisional_) +
						std::pow(DcNUCL, 2.) * GetMissingMoment(-1. / 2., 0.) +
						2.0*DcNUCL*K_collisional_*GetMissingMoment(Av_collisional_ - 0.5, As_collisional_));

			double psi1 = 	std::pow(V0_, -1. / 2.)*(std::pow(K_collisional_, 2.) * GetMissingMoment(2.0*Av_collisional_ + 0.5, 2.0*As_collisional_) +
					std::pow(DcNUCL, 2.) * GetMissingMoment(1. / 2., 0.) +
					2.0*DcNUCL*K_collisional_*GetMissingMoment(Av_collisional_ + 0.5, As_collisional_));

			psi1 = psi1 + V0_*psi0;

			source_coagulation_sl_(0) = 0.;
			if (psi0*psi1>=0.)
				source_coagulation_sl_(0) = -2.2 * Cfm * std::sqrt(T_) * N0_ / AvogadroNumber * std::sqrt(psi0*psi1);
		}

		source_coagulation_sl_(1) = 0.0;

		{
			const double psi0 = 	std::pow(V0_, -1. / 2.)*(std::pow(K_collisional_, 2.) * GetMissingMoment(2.0*Av_collisional_ - 0.5 + Av_fractal_, 2.0*As_collisional_ + As_fractal_ + 1.) +
						std::pow(DcNUCL, 2.) * GetMissingMoment(-1. / 2. + Av_fractal_, As_fractal_ + 1.) +
						2.*DcNUCL*K_collisional_*GetMissingMoment(Av_collisional_ - 0.5 + Av_fractal_, As_collisional_ + As_fractal_ + 1.));

			double psi1 = 	std::pow(V0_, -1. / 2.)*(std::pow(K_collisional_, 2.) * GetMissingMoment(2.0*Av_collisional_ + 0.5 + Av_fractal_, 2.0*As_collisional_ + As_fractal_ + 1.) +
					std::pow(DcNUCL, 2.) * GetMissingMoment(1. / 2. + Av_fractal_, As_fractal_ + 1.) +
					2 * DcNUCL*K_collisional_*GetMissingMoment(Av_collisional_ + 0.5 + Av_fractal_, As_collisional_ + As_fractal_ + 1.));
			
			psi1 = psi1 + V0_*psi0;

			source_coagulation_sl_(2) = source_coagulation_sl_(0);
			if (psi0*psi1>=0.)
				source_coagulation_sl_(2) += 	2.2*Cfm*std::sqrt(T_)*N0_*V0_*K_fractal_ / S0_ / AvogadroNumber*std::sqrt(psi0*psi1);
		}

		source_coagulation_sl_(3) = source_coagulation_sl_(0);
	}

	void HMOM::SootCoagulationLargeLargeM4()
	{
		source_coagulation_ll_(0) = 0.;
		source_coagulation_ll_(1) = 0.;
		source_coagulation_ll_(2) = 0.;
		source_coagulation_ll_(3) = 0.;

		const double psi0 = 	2. * std::pow(K_collisional_, 2.)*(GetMissingMoment(2.*Av_collisional_ - 0.5, 2.*As_collisional_)*GetMissingMoment(-1. / 2., 0.) +
					GetMissingMoment(Av_collisional_ - 0.5, As_collisional_)* GetMissingMoment(Av_collisional_ - 0.5, As_collisional_));
		const double psi1 = 	2. * std::pow(K_collisional_, 2.)*(GetMissingMoment(2.*Av_collisional_ + 0.5, 2.*As_collisional_)*GetMissingMoment(-1. / 2., 0.) +
					GetMissingMoment(Av_collisional_ + 0.5, As_collisional_)* GetMissingMoment(Av_collisional_ - 0.5, As_collisional_) +
					GetMissingMoment(2.0*Av_collisional_ - 0.5, 2.0*As_collisional_)*GetMissingMoment(1. / 2., 0.));

		if (psi0*psi1>=0.)
			source_coagulation_ll_(0) = -0.5 * 2.2 * Cfm * std::sqrt(T_) / AvogadroNumber * std::sqrt(psi1*psi0);			
	}

	void HMOM::SootCoagulationContinousM4()
	{
		double lambda = 8.2057e-5 / std::sqrt(2.) / std::pow(200.e-12, 2.) / AvogadroNumber;
		lambda = 1.257 * lambda * T_ / (P_Pa_ / 101325.);

		SootCoagulationContinousSmallSmallM4(lambda);
		SootCoagulationContinousSmallLargeM4(lambda);
		SootCoagulationContinousLargeLargeM4(lambda);
	}

	void HMOM::SootCoagulationContinousSmallSmallM4(const double lambda)
	{
		const double DcNUCL = K_diam * std::pow(V0_, 1. / 3.);
		const double S00 = K_spher * std::pow(2.0*V0_, 2. / 3.);
		const double CC0 = 1. + lambda / DcNUCL;

		const double beta00 = 2.0*KB*T_ / 3.0 / viscosity_*(2 * CC0 / DcNUCL)*(2.0*DcNUCL);

		source_coagulation_continous_ss_(0) = -0.5* beta00 * std::pow(N0_, 2.) / AvogadroNumber;
		source_coagulation_continous_ss_(1) = 0.;
		source_coagulation_continous_ss_(2) = 0.5* beta00 * S00 / S0_ * std::pow(N0_, 2.) / AvogadroNumber + 2.0* source_coagulation_continous_ss_(0);
		source_coagulation_continous_ss_(3) = 2.*source_coagulation_continous_ss_(0);
	}

	void HMOM::SootCoagulationContinousSmallLargeM4(const double lambda)
	{
		const double DcNUCL = K_diam * std::pow(V0_, 1. / 3.);
		const double betai0 = 2.0*KB*T_ / 3.0 / viscosity_;

		source_coagulation_continous_sl_(0) = (2.0 + lambda / DcNUCL)*GetMissingMoment(0., 0.) +
			(DcNUCL + lambda) / K_collisional_ * GetMissingMoment(-1.0*Av_collisional_, -1.0*As_collisional_) +
			(1.0 + lambda / DcNUCL)*K_collisional_ / DcNUCL * GetMissingMoment(1.0*Av_collisional_, 1.0*As_collisional_) +
			lambda * DcNUCL* std::pow(K_collisional_, 2.) * GetMissingMoment(-2.*Av_collisional_, -2.*As_collisional_);

		source_coagulation_continous_sl_(0) = -N0_ / AvogadroNumber * betai0 * source_coagulation_continous_sl_(0);

		source_coagulation_continous_sl_(1) = 0.;

		source_coagulation_continous_sl_(2) = (2.0 + lambda / DcNUCL)*GetMissingMoment(Av_fractal_, As_fractal_ + 1.0) +
			(DcNUCL + lambda) / K_collisional_ * GetMissingMoment(-1.0*Av_collisional_ + Av_fractal_, -1.0*As_collisional_ + As_fractal_ + 1.0) +
			(1.0 + lambda / DcNUCL)*K_collisional_ / DcNUCL * GetMissingMoment(1.0*Av_collisional_ + Av_fractal_, 1.0*As_collisional_ + As_fractal_ + 1.0) +
			lambda * DcNUCL* std::pow(K_collisional_, 2.) * GetMissingMoment(-2.0*Av_collisional_ + Av_fractal_, -2.0*As_collisional_ + As_fractal_ + 1.0);

		source_coagulation_continous_sl_(2) = N0_*V0_*K_fractal_ / S0_ / AvogadroNumber * betai0* source_coagulation_continous_sl_(2) + source_coagulation_continous_sl_(0);

		source_coagulation_continous_sl_(3) = source_coagulation_continous_sl_(0);
	}

	void HMOM::SootCoagulationContinousLargeLargeM4(const double lambda)
	{
		const double DcNUCL = K_diam * std::pow(V0_, 1. / 3.);
		const double betai0 = 2.0*KB*T_ / 3.0 / viscosity_;

		source_coagulation_continous_ll_(0) = GetMissingMoment(0., 0.) * GetMissingMoment(0., 0.) +
			lambda / K_collisional_ * (GetMissingMoment(0., 0.) * GetMissingMoment(-1.*Av_collisional_, -1.*As_collisional_) +
				GetMissingMoment(Av_collisional_, As_collisional_) * GetMissingMoment(-2.*Av_collisional_, -2.*As_collisional_) +
				GetMissingMoment(Av_collisional_, As_collisional_) * GetMissingMoment(-1.*Av_collisional_, -1.*As_collisional_));

		source_coagulation_continous_ll_(0) = -0.5 * betai0 / AvogadroNumber * source_coagulation_continous_ll_(0);

		source_coagulation_continous_ll_(1) = 0.0;
		source_coagulation_continous_ll_(2) = 0.0;
		source_coagulation_continous_ll_(3) = 0.0;
	}

	void HMOM::GetMoments()
	{
		// Denormalize moments
		M00_ = M00_normalized_*AvogadroNumber;			// number density [#/m3]
		M10_ = M10_normalized_*V0_*AvogadroNumber;		// total soot volume [#]
		M01_ = M01_normalized_*S0_*AvogadroNumber;		// total soot surface area [#/m]
		N0_ = N0_normalized_*AvogadroNumber;			// weight of delta function [#/m3]

		// Properties of large particles (equations 8,9 and 10)
		// Mueller et al., Comb. & Flame 156, pp.1143-1155 (2009)
		// NL = M00-N0
		// VL = (M10-N0*V0)/NL
		// SL = (M01-N0*S0)/NL

		NL_ = M00_ - N0_;						//!< number density associated to large particles (2nd mode) [#/m3]
		NLVL_ = M10_ - V0_*N0_;		//!< mean volume of large particles (2nd mode) multiplied by NL [#]
		NLSL_ = M01_ - S0_*N0_;		//!< mean surface of large particles (2nd mode) multiplied by NL [#/m]
	}

	double HMOM::GetMoment(const double i, const double j) const
	{
		double moment = 0.;
		if (NL_ < 1.e-25 || NLVL_ < 1.e-25 || NLSL_ < 1.e-25)
		{
			if ((M00_ != 0.) && (M10_ != 0.) && (M01_ != 0.))
				moment = M00_ * std::pow(M10_ / M00_, i) * std::pow(M01_ / M00_, j);
		}
		else
		{
			// Moments (equation 7)
			// Mueller et al., Comb. & Flame 156, pp.1143-1155 (2009)
			// M(x,y) = N0*V0^x*S0^y + NL*VL^x*SL^y

			moment =	N0_*std::pow(V0_, i)*std::pow(S0_, j) +
					NL_*std::pow(NLVL_ / NL_, i)*std::pow(NLSL_ / NL_, j);
		}

		return moment;
	}

	double HMOM::GetMissingMoment(const double i, const double j) const
	{
		const double coefficient = std::pow(V0_, i) * std::pow(S0_, j);

		if (NL_ < 1.e-25 || NLVL_ < 1.e-25 || NLSL_ < 1.e-25)
			return (1.e-25*coefficient);
		else
			return GetMoment(i, j) - N0_*coefficient;
	}

	double HMOM::SootVolumeFraction() const
	{
		return GetMoment(1., 0.);
	}

	double HMOM::SootParticleDiameter() const
	{
		// Diameter of primary particles (including large aggregate) (equation 1)
		// Mueller et al., Comb. & Flame 156, pp.1143-1155 (2009)
		// dp = 6*V/S

		return 6.*GetMoment(1., -1) / (GetMoment(0., 0.)+1e-32);
	}

	double HMOM::SootCollisionParticleDiameter() const
	{
		return K_collisional_*GetMoment(Av_collisional_, As_collisional_) / (GetMoment(0., 0.)+1e-32);
	}

	double HMOM::SootParticleNumberDensity() const
	{
		return GetMoment(0., 0.);
	}

	double HMOM::SootNumberOfPrimaryParticles() const
	{
		// Primary particle number np (equation 2)
		// Mueller et al., Comb. & Flame 156, pp.1143-1155 (2009)
		// np = 1/(36pi)*V^(-2)*S(3)

		return std::pow(K_spher, -3.)*GetMoment(-2., 3.) / (GetMoment(0., 0.)+1e-32);
	}

	double HMOM::planck_coefficient(const double T, const double fv) const
	{
		if (soot_planck_coefficient_ == SOOT_PLANCK_COEFFICIENT_SMOOKE)
			return (1307.*fv*T);							// [1/m]	(Smooke et al. Combustion and Flame 2009)
		else if (soot_planck_coefficient_ == SOOT_PLANCK_COEFFICIENT_KENT)
			return (2262.*fv*T);							// [1/m]	(Kent al. Combustion and Flame 1990)
		else if (soot_planck_coefficient_ == SOOT_PLANCK_COEFFICIENT_SAZHIN)
			return (1232.*(1. + 4.8e-4*(T - 2000.)));		// [1/m]	(Sazhin, Fluent 1994)
		else
			return 0.;
	}

	HMOM::~HMOM()
	{
	}

}
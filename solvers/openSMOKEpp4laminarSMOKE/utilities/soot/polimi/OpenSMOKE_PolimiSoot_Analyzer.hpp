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

namespace OpenSMOKE
{
	#include <sstream>

	PolimiSoot_Analyzer::PolimiSoot_Analyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML) :
		thermo_(*thermodynamicsMapXML)
	{
		//	Default values
		bin_label_			= "BIN";
		bin_minimum_		= "BIN5";
		bin_minimum_fractal_dimension_ = "BIN12";
		bin_index_zero_		= 10;
		bin_index_final_	= 20;
		bin_density_zero_   = 1500.;
		bin_density_final_  = 1700.;
		Df_					= 1.8;
		physical_diffusion_ = true;
		physical_diffusion_reduction_coefficient_ = 1.;

		thermophoretic_effect_ = true;
		thermophoretic_effect_amplification_factor_ = 1.;
		thermophoretic_effect_included_in_correction_ = false;
		thermophoretic_effect_in_enthalpy_fluxes_ = DO_NOT_EXCLUDE_SOOT_CONTRIBUTION;
		thermophoretic_effect_smoothing_time_ = 0.;
		radiative_heat_transfer_ = true;
		soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_SMOOKE;

		write_psdf_ = true;
		threshold_for_psdf_ = 1e-11;
	}

	PolimiSoot_Analyzer::PolimiSoot_Analyzer(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML, OpenSMOKE::OpenSMOKE_Dictionary& dictionary) :
		thermo_(*thermodynamicsMapXML)
	{
		SetupFromDictionary(dictionary);
	}

	void PolimiSoot_Analyzer::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		// Default values
		bin_index_zero_ = 10;
		bin_index_final_ = 20;
		bin_density_zero_ = 1500.;
		bin_density_final_ = 1700.;
		Df_ = 1.8;
		bin_label_ = "BIN";
		bin_minimum_ = "BIN5";
		bin_minimum_fractal_dimension_ = "BIN12";
		physical_diffusion_ = true;
		physical_diffusion_reduction_coefficient_ = 1.;
		
		// These values will be overwritten by user-definitions (thermophoresis)
		thermophoretic_effect_ = true;
		thermophoretic_effect_amplification_factor_ = 1.;
		thermophoretic_effect_included_in_correction_ = false;
		thermophoretic_effect_in_enthalpy_fluxes_ = DO_NOT_EXCLUDE_SOOT_CONTRIBUTION;
		thermophoretic_effect_smoothing_time_ = 0.;
		
		// These values will be overwritten by user-definitions (radiative heat transfer)
		radiative_heat_transfer_ = true;
		soot_planck_coefficient_ = SOOT_PLANCK_COEFFICIENT_SMOOKE;
		
		// Read from dictionary
		{
			Grammar_PolimiSoot_Analyzer grammar;
			dictionary.SetGrammar(grammar);

			if (dictionary.CheckOption("@FractalDimension") == true)
				dictionary.ReadDouble("@FractalDimension", Df_);

			if (dictionary.CheckOption("@SootLabel") == true)
				dictionary.ReadString("@SootLabel", bin_label_);

			if (dictionary.CheckOption("@SootMinimumSection") == true)
			{
				int minimum_section;
				dictionary.ReadInt("@SootMinimumSection", minimum_section);

				std::stringstream number;
				number << minimum_section;
				bin_minimum_ = bin_label_ + number.str();
			}

			if (dictionary.CheckOption("@SootMinimumFractalDimension") == true)
			{
				int minimum_section;
				dictionary.ReadInt("@SootMinimumFractalDimension", minimum_section);

				std::stringstream number;
				number << minimum_section;
				bin_minimum_fractal_dimension_ = bin_label_ + number.str();
			}

			if (dictionary.CheckOption("@ThermophoreticEffect") == true)
				dictionary.ReadBool("@ThermophoreticEffect", thermophoretic_effect_);

			if (dictionary.CheckOption("@ThermophoreticEffectAmplificationFactor") == true)
				dictionary.ReadDouble("@ThermophoreticEffectAmplificationFactor", thermophoretic_effect_amplification_factor_);			

			if (dictionary.CheckOption("@ThermophoreticEffectInCorrectionVelocity") == true)
				dictionary.ReadBool("@ThermophoreticEffectInCorrectionVelocity", thermophoretic_effect_included_in_correction_);

			if (dictionary.CheckOption("@ThermophoreticEffectSmoothingTime") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@ThermophoreticEffectSmoothingTime", thermophoretic_effect_smoothing_time_, units);

				if (units == "s")			thermophoretic_effect_smoothing_time_ = thermophoretic_effect_smoothing_time_;
				else if (units == "ms")		thermophoretic_effect_smoothing_time_ /= 1000.;
				else if (units == "min")	thermophoretic_effect_smoothing_time_ *= 60.;
				else OpenSMOKE::FatalErrorMessage("@ThermophoreticEffectSmoothingTime: check units of time");
			}

			if (dictionary.CheckOption("@ThermophoreticEffectInEnthalpyFluxes") == true)
			{
				std::string flag;
				dictionary.ReadString("@ThermophoreticEffectInEnthalpyFluxes", flag);

				if (flag == "DoNotExclude")
					thermophoretic_effect_in_enthalpy_fluxes_ = DO_NOT_EXCLUDE_SOOT_CONTRIBUTION;
				else if (flag == "Exclude")
					thermophoretic_effect_in_enthalpy_fluxes_ = EXCLUDE_TOTAL_SOOT_CONTRIBUTION;
				else if (flag == "ExcludeOnlyThermophoreticEffect")
					thermophoretic_effect_in_enthalpy_fluxes_ = EXCLUDE_THERMOPHORETIC_SOOT_CONTRIBUTION;
				else
					OpenSMOKE::FatalErrorMessage("@ThermophoreticEffectInEnthalpyFluxes: available options: DoNotExclude | Exclude | ExcludeOnlyThermophoreticEffect");
			}

			if (dictionary.CheckOption("@RadiativeHeatTransfer") == true)
				dictionary.ReadBool("@RadiativeHeatTransfer", radiative_heat_transfer_);

			if (dictionary.CheckOption("@PhysicalDiffusion") == true)
				dictionary.ReadBool("@PhysicalDiffusion", physical_diffusion_);

			if (dictionary.CheckOption("@PhysicalDiffusionReductionCoefficient") == true)
			{
				dictionary.ReadDouble("@PhysicalDiffusionReductionCoefficient", physical_diffusion_reduction_coefficient_);
			}

			if (dictionary.CheckOption("@PlanckCoefficient") == true)
			{
				std::string flag;
				dictionary.ReadString("@PlanckCoefficient", flag);
				SetPlanckAbsorptionCoefficient(flag);
			}

			if (dictionary.CheckOption("@WritePSDF") == true)
				dictionary.ReadBool("@WritePSDF", write_psdf_);

			if (dictionary.CheckOption("@ThresholdForPSDF") == true)
				dictionary.ReadDouble("@ThresholdForPSDF", threshold_for_psdf_);
		}

		// Setup
		Setup();
	}

	void PolimiSoot_Analyzer::SetLabel(const std::string label)
	{
		bin_label_ = label;
	}

	void PolimiSoot_Analyzer::SetFractalDiameter(const double Df)
	{
		Df_ = Df;
	}
	void PolimiSoot_Analyzer::SetMinimumSection(const int minimum_section)
	{
		std::stringstream number;
		number << minimum_section;
		bin_minimum_ = bin_label_ + number.str();
	}
	void PolimiSoot_Analyzer::SetMinimumSection(const std::string bin_minimum)
	{
		bin_minimum_ = bin_minimum;
	}

	void PolimiSoot_Analyzer::SetDensity(const int bin_index_zero, const int bin_index_final, const double bin_density_zero, const double bin_density_final)
	{
		bin_index_zero_ = bin_index_zero;
		bin_index_final_ = bin_index_final;
		bin_density_zero_ = bin_density_zero;
		bin_density_final_ = bin_density_final;
	}

	void PolimiSoot_Analyzer::SetPlanckAbsorptionCoefficient(const SootPlanckCoefficient soot_planck_coefficient)
	{
		soot_planck_coefficient_ = soot_planck_coefficient;
	}

	void PolimiSoot_Analyzer::SetPlanckAbsorptionCoefficient(const std::string label)
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

	void PolimiSoot_Analyzer::Setup()
	{
		iC_ = thermo_.IndexOfElement("C") - 1;
		iH_ = thermo_.IndexOfElement("H") - 1;
		iO_ = thermo_.IndexOfElement("O") - 1;

		nspecies_ = thermo_.NumberOfSpecies();

		// Check is soot sections are really available in the kinetic mechanism
		{
			unsigned int count = 0;
			for (unsigned int i = 0; i < nspecies_; i++)
			if (thermo_.NamesOfSpecies()[i].compare(0, bin_label_.size(), bin_label_) == 0)
				count++;

			if (count == 0)
				OpenSMOKE::FatalErrorMessage("No soot sections (" + bin_label_ + ") were found in the provided kinetic mechanism");
		}

		// Check existence of minimum section describing soot
		double min_mw_soot = 1.e32;
		{
			bool iFound = false;
			for (unsigned int i = 0; i < nspecies_; i++)
				if (thermo_.NamesOfSpecies()[i].compare(0, bin_minimum_.size(), bin_minimum_) == 0)
				{
					iFound = true;
					if (thermo_.MW(i) < min_mw_soot)	min_mw_soot = thermo_.MW(i);
				}

			if (iFound == false)
				OpenSMOKE::FatalErrorMessage("No minimum soot section (" + bin_minimum_ + ") was found in the provided kinetic mechanism");
		}

		unsigned int bin_minimum_fractal_dimension_index = 0;
		for (unsigned int i = 0; i < nspecies_; i++)
		if (thermo_.NamesOfSpecies()[i].compare(0, bin_label_.size(), bin_label_) == 0)
		{
			const double nc = thermo_.atomic_composition()(i, iC_);
			const double nh = thermo_.atomic_composition()(i, iH_);
			const double no = thermo_.atomic_composition()(i, iO_);

			// BIN density
			{
				const double index = std::log(nc / 24.) / std::log(2.) + 1.;

				if (index <= bin_index_zero_)
				{
					bin_density_.push_back(bin_density_zero_);
				}
				else
				{
					const double m = (bin_density_final_ - bin_density_zero_) / (bin_index_final_ - bin_index_zero_);
					const double c = bin_density_zero_ - m*bin_index_zero_;
					bin_density_.push_back(c + m*index);
				}
			}

			if (thermo_.NamesOfSpecies()[i].compare(0, bin_minimum_fractal_dimension_.size(), bin_minimum_fractal_dimension_) == 0 && bin_minimum_fractal_dimension_index == 0)
				bin_minimum_fractal_dimension_index = bin_indices_.size();

			bin_indices_.push_back(i);						// Index of bin in the gas phase kinetic scheme [-]
			bin_names_.push_back(thermo_.NamesOfSpecies()[i]);			// Name of bin in the gas phase kinetic scheme [-]
			bin_mw_.push_back(thermo_.MW(i));					// Molecular weight [kg/kmol]
			bin_m_.push_back(thermo_.MW(i) / PhysicalConstants::Nav_kmol);		// Mass of particle [kg]

			bin_ds_.push_back(std::pow(6. / PhysicalConstants::pi*thermo_.MW(i) / (bin_density_[bin_density_.size() - 1] / 1000.)
				/ (PhysicalConstants::Nav_mol), 1. / 3.)*1.e-2);	// Diameter of particle [m]

			bin_V_.push_back(PhysicalConstants::pi / 6.*std::pow(bin_ds_[bin_ds_.size() - 1], 3.));	// Volume of particle [m3]
			bin_c_.push_back(nc);	// C
			bin_h_.push_back(nh);// H
			bin_o_.push_back(no);// O
			bin_h_over_c_.push_back(nh / nc);	// Ratio H/C
			bin_o_over_c_.push_back(no / nc);// Ratio O/C

			if (nh > 0)	bin_o_over_h_.push_back(no / nh);	// Ratio O/H
			else    	bin_o_over_h_.push_back(0.);	// Ratio O/H

			if ( thermo_.MW(i) >=  min_mw_soot)	
			{
				int index = bin_indices_.size()-1;
				bin_indices_large_.push_back(index);
				bin_indices_large_global_.push_back(i);
				bin_density_large_.push_back(bin_density_[index]);
				bin_V_large_.push_back(bin_V_[index]);
			}
			else
			{
				int index = bin_indices_.size()-1;
				bin_indices_small_.push_back(index);
				bin_indices_small_global_.push_back(i);
				bin_density_small_.push_back(bin_density_[index]);
				bin_V_small_.push_back(bin_V_[index]);
			}

			// Collisional diameter and diameter
			{
				bool iCollisional = false;
				const double index = std::log(nc / 24.) / std::log(2.) + 1;

				if (index > 12)
					iCollisional = true;

				if (iCollisional == true)
				{
					bin_np_.push_back(bin_m_[bin_m_.size() - 1] / bin_m_[bin_minimum_fractal_dimension_index]);
					bin_dc_.push_back(std::sqrt(5. / 3.)*bin_ds_[bin_minimum_fractal_dimension_index] * std::pow(bin_np_[bin_np_.size() - 1] / std::pow(1. + 2. / Df_, Df_ / 2.), 1. / Df_));
					bin_d_.push_back(bin_dc_[bin_dc_.size() - 1]);

					if (thermo_.MW(i) >= min_mw_soot)
					{
						bin_indices_large_aggregates_.push_back(bin_indices_.size() - 1);
						bin_indices_large_aggregates_global_.push_back(i);
					}
				}
				else
				{
					bin_np_.push_back(0.);
					bin_dc_.push_back(0.);
					bin_d_.push_back(bin_ds_[bin_ds_.size() - 1]);

					if (thermo_.MW(i) >= min_mw_soot)
					{
						bin_indices_large_spherical_.push_back(bin_indices_.size() - 1);
						bin_indices_large_spherical_global_.push_back(i);
					}
				}
			}
		}

		// Physical diffusivities (correction factors)
		{
			const double bin_to_cut = 10.;
			const int binReference = std::min_element(bin_mw_.begin(), bin_mw_.end()) - bin_mw_.begin();
			const double MWReference = bin_mw_[binReference];

			if (binReference >= 0)
			{
				bin_diffusivity_correction_factors_.resize(bin_indices_.size());
				for (unsigned int i = 0; i < bin_indices_.size(); i++)
				{
					double MWratio = std::min(bin_mw_[i] / MWReference, std::pow(2., bin_to_cut));
					bin_diffusivity_correction_factors_[i] = std::pow(MWratio, -0.681);
				}
			}

			bin_diffusivity_reference_species_ = bin_indices_[binReference];
		}

		// Write summary
		WriteBinData();

		// Final operations
		{
			// Memory allocation
			bin_omega_.resize(bin_indices_.size());
			bin_x_.resize(bin_indices_.size());
			bin_fv_.resize(bin_indices_.size());
			bin_rho_.resize(bin_indices_.size());
			bin_N_.resize(bin_indices_.size());

			//if (iBin_ == true)
			{
				for (unsigned int i = 0; i < bin_indices_.size(); i++)
				{
					bool iNew = true;
					for (unsigned int k = 0; k < bin_baskets_.size(); k++)
					if (bin_c_[i] == bin_baskets_[k])
					{
						iNew = false;
						break;
					}
					if (iNew == true)
						bin_baskets_.push_back(bin_c_[i]);
				}

				std::sort(bin_baskets_.begin(), bin_baskets_.end());
				bin_baskets_indices_.resize(bin_baskets_.size());


				for (unsigned int i = 0; i < bin_indices_.size(); i++)
				for (unsigned int k = 0; k < bin_baskets_.size(); k++)
				{
					if (bin_c_[i] == bin_baskets_[k])
						bin_baskets_indices_[k].push_back(i);
				}

				bin_baskets_d_.resize(bin_baskets_.size());
				bin_baskets_mw_.resize(bin_baskets_.size());
				bin_baskets_log10d_.resize(bin_baskets_.size());
				bin_baskets_dlog10d_.resize(bin_baskets_.size());
				dN_over_dlog10d_.resize(bin_baskets_.size());

				bin_baskets_V_.resize(bin_baskets_.size());
				bin_baskets_log10V_.resize(bin_baskets_.size());
				bin_baskets_dlog10V_.resize(bin_baskets_.size());
				dN_over_dlog10V_.resize(bin_baskets_.size());

				bin_baskets_m_.resize(bin_baskets_.size());
				bin_baskets_log10m_.resize(bin_baskets_.size());
				bin_baskets_dlog10m_.resize(bin_baskets_.size());
				dN_over_dlog10m_.resize(bin_baskets_.size());

				bin_baskets_N_.resize(bin_baskets_.size());
				bin_baskets_fv_.resize(bin_baskets_.size());
				bin_baskets_rho_.resize(bin_baskets_.size());
				bin_baskets_x_.resize(bin_baskets_.size());
				bin_baskets_omega_.resize(bin_baskets_.size());
			}

			for (unsigned int k = 0; k < bin_baskets_.size(); k++)
			{
				bin_baskets_mw_[k] = 0.;
				bin_baskets_d_[k] = 0.;
				bin_baskets_m_[k] = 0.;
				bin_baskets_V_[k] = 0.;

				for (unsigned int i = 0; i < bin_baskets_indices_[k].size(); i++)
				{
					int j = bin_baskets_indices_[k][i];
					bin_baskets_d_[k] += bin_d_[j];
					bin_baskets_mw_[k] += bin_mw_[j];
					bin_baskets_V_[k] += bin_V_[j];
					bin_baskets_m_[k] += bin_m_[j];
				}

				bin_baskets_d_[k] /= double(bin_baskets_indices_[k].size());
				bin_baskets_mw_[k] /= double(bin_baskets_indices_[k].size());
				bin_baskets_V_[k] /= double(bin_baskets_indices_[k].size());
				bin_baskets_m_[k] /= double(bin_baskets_indices_[k].size());

				bin_baskets_log10d_[k] = std::log10(bin_baskets_d_[k]);
				bin_baskets_log10V_[k] = std::log10(bin_baskets_V_[k]);
				bin_baskets_log10m_[k] = std::log10(bin_baskets_m_[k]);
			}

			// Intervals
			bin_baskets_dlog10d_[0] = ((bin_baskets_log10d_[1] + bin_baskets_log10d_[0]) / 2. - bin_baskets_log10d_[0])*2.;
			bin_baskets_dlog10V_[0] = ((bin_baskets_log10V_[1] + bin_baskets_log10V_[0]) / 2. - bin_baskets_log10V_[0])*2.;
			bin_baskets_dlog10m_[0] = ((bin_baskets_log10m_[1] + bin_baskets_log10m_[0]) / 2. - bin_baskets_log10m_[0])*2.;

			for (unsigned int k = 1; k < bin_baskets_.size() - 1; k++)
			{
				bin_baskets_dlog10d_[k] = (bin_baskets_log10d_[k + 1] + bin_baskets_log10d_[k]) / 2. - (bin_baskets_log10d_[k] + bin_baskets_log10d_[k - 1]) / 2.;
				bin_baskets_dlog10V_[k] = (bin_baskets_log10V_[k + 1] + bin_baskets_log10V_[k]) / 2. - (bin_baskets_log10V_[k] + bin_baskets_log10V_[k - 1]) / 2.;
				bin_baskets_dlog10m_[k] = (bin_baskets_log10m_[k + 1] + bin_baskets_log10m_[k]) / 2. - (bin_baskets_log10m_[k] + bin_baskets_log10m_[k - 1]) / 2.;
			}

			int k = bin_baskets_.size() - 1;
			bin_baskets_dlog10d_[k] = ((bin_baskets_log10d_[k] + bin_baskets_log10d_[k - 1]) / 2. - bin_baskets_log10d_[k - 1])*2.;
			bin_baskets_dlog10V_[k] = ((bin_baskets_log10V_[k] + bin_baskets_log10V_[k - 1]) / 2. - bin_baskets_log10V_[k - 1])*2.;
			bin_baskets_dlog10m_[k] = ((bin_baskets_log10m_[k] + bin_baskets_log10m_[k - 1]) / 2. - bin_baskets_log10m_[k - 1])*2.;
		}

		// Soot dimer
		for(unsigned int i=0;i<nspecies_;i++)
			if (thermo_.NamesOfSpecies()[i].compare(0, bin_minimum_.size(), bin_minimum_) == 0)
				soot_dimer_indices_global_.push_back(i);

		// PAH (340 nm, 1/2 aromatic rings)
		pah_1_2_rings_indices_global_.resize(0);
		if ( thermo_.IndexOfSpeciesWithoutError("C6H6") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H6")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C7H8") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C7H8")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("INDENE") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("INDENE")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C10H8") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C10H8")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C12H8") > 0 )		pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C12H8")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("BIPHENYL") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("BIPHENYL")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("FLUORENE") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("FLUORENE")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5C2H") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5C2H")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5C2H3") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5C2H3")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5C2H5") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5C2H5")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C6H5CH2C6H5") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C6H5CH2C6H5")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C10H7CH3") > 0 )	pah_1_2_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C10H7CH3")-1 );

		// PAH (340 nm, 3/4 aromatic rings)
		pah_3_4_rings_indices_global_.resize(0);
		if ( thermo_.IndexOfSpeciesWithoutError("C14H10") > 0 )		pah_3_4_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C14H10")-1 );
		if ( thermo_.IndexOfSpeciesWithoutError("C16H10") > 0 )		pah_3_4_rings_indices_global_.push_back( thermo_.IndexOfSpeciesWithoutError("C16H10")-1 );

		// PAH (340 nm, more than 4 aromatic rings)
		// Actually, they are equivalent to small BINs
		pah_more_than_4_rings_indices_global_.resize(0);
		pah_more_than_4_rings_indices_global_ = bin_indices_small_global_;	

		// Summary PAH data
		std::cout << "PAHs with 1/2 rings (340 nm)" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<pah_1_2_rings_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[pah_1_2_rings_indices_global_[i]] << " " << pah_1_2_rings_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "PAHs with 3/4 rings (400 nm)" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<pah_3_4_rings_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[pah_3_4_rings_indices_global_[i]] << " " << pah_3_4_rings_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "PAHs with more than 4 rings (500 nm)" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<pah_more_than_4_rings_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[pah_more_than_4_rings_indices_global_[i]] << " " << pah_more_than_4_rings_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "Soot dimers" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<soot_dimer_indices_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[soot_dimer_indices_global_[i]] << " " << soot_dimer_indices_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "Soot particles" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for(unsigned int i=0;i<bin_indices_large_global_.size();i++)
			std::cout << thermo_.NamesOfSpecies()[bin_indices_large_global_[i]] << " " << bin_indices_large_global_[i] << std::endl;	
		std::cout << std::endl;

		std::cout << "Soot spherical particles" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i<bin_indices_large_spherical_global_.size(); i++)
			std::cout << thermo_.NamesOfSpecies()[bin_indices_large_spherical_global_[i]] << " " << bin_indices_large_spherical_global_[i] << std::endl;
		std::cout << std::endl;

		std::cout << "Soot aggregates" << std::endl;
		std::cout << "------------------------------------------------" << std::endl;
		for (unsigned int i = 0; i<bin_indices_large_aggregates_global_.size(); i++)
			std::cout << thermo_.NamesOfSpecies()[bin_indices_large_aggregates_global_[i]] << " " << bin_indices_large_aggregates_global_[i] << std::endl;
		std::cout << std::endl;
	}

	double PolimiSoot_Analyzer::fv_large(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < bin_indices_large_global_.size(); i++)
		{
			int j = bin_indices_large_global_[i];
			sum += rhoGas*omegaGas(j) / bin_density_large_[i];
		}
		return sum;
	}

	double PolimiSoot_Analyzer::rho_large(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < bin_indices_large_global_.size(); i++)
		{
			int j = bin_indices_large_global_[i];
			sum += rhoGas*omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::fv_small(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < bin_indices_small_global_.size(); i++)
		{
			int j = bin_indices_small_global_[i];
			sum += rhoGas*omegaGas(j) / bin_density_small_[i];
		}
		return sum;
	}

	double PolimiSoot_Analyzer::fv(const double rhoGas, const Eigen::VectorXd &omegaGas) const
	{
		return (fv_small(rhoGas, omegaGas) + fv_large(rhoGas, omegaGas));
	}

	double PolimiSoot_Analyzer::omega_pah_1_2_rings(const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_1_2_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_1_2_rings_indices_global_[i];
			sum += omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::omega_pah_3_4_rings(const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_3_4_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_3_4_rings_indices_global_[i];
			sum += omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::omega_pah_more_than_4_rings(const Eigen::VectorXd &omegaGas) const
	{
		double sum = 0.;
		for (unsigned int i = 0; i < pah_more_than_4_rings_indices_global_.size(); i++)
		{
			const unsigned int j = pah_more_than_4_rings_indices_global_[i];
			sum += omegaGas(j);
		}
		return sum;
	}

	double PolimiSoot_Analyzer::planck_coefficient(const double rhoGas, const double T, const Eigen::VectorXd &omegaGas) const
	{
		if (soot_planck_coefficient_ == SOOT_PLANCK_COEFFICIENT_SMOOKE)
			return (1307.*fv_large(rhoGas, omegaGas)*T);							// [1/m]	(Smooke et al. Combustion and Flame 2009)
		else if (soot_planck_coefficient_ == SOOT_PLANCK_COEFFICIENT_KENT)
			return (2262.*fv_large(rhoGas, omegaGas)*T);							// [1/m]	(Kent al. Combustion and Flame 1990)
		else if (soot_planck_coefficient_ == SOOT_PLANCK_COEFFICIENT_SAZHIN)
			return (1232.*fv_large(rhoGas, omegaGas)*(1. + 4.8e-4*(T - 2000.)));	// [1/m]	(Sazhin, Fluent 1994)
		else
			return 0.;
	}

	void PolimiSoot_Analyzer::Analysis(const double T, const double P_Pa, const double rhoGas, const Eigen::VectorXd &omegaGas, const Eigen::VectorXd &xGas)
	{
		//if (iBin_ == true)
		{
			double small_eps = 1e-20;

			for (unsigned int i = 0; i < bin_indices_.size(); i++)
			{
				int j = bin_indices_[i];
				bin_omega_[i] = omegaGas(j);				// mass fraction
				bin_x_[i] = xGas(j);				// mole fraction
				bin_rho_[i] = rhoGas*omegaGas(j);			// density [kg_soot/m3]
				bin_fv_[i] = bin_rho_[i] / bin_density_[i];		// volume fraction [m3_soot/m3]
				bin_N_[i] = bin_fv_[i] / bin_V_[i];			// [1/m3]
			}

			fv_small_ = 0.;
			rho_small_ = 0.;
			N_small_ = 0.;
			omega_small_ = 0.;
			x_small_ = 0.;
			h_over_c_small_ = 0.;
			o_over_c_small_ = 0.;
			o_over_h_small_ = 0.;

			for (unsigned int i = 0; i < bin_indices_small_.size(); i++)
			{
				int j = bin_indices_small_[i];
				fv_small_ += bin_fv_[j];
				rho_small_ += bin_rho_[j];
				N_small_ += bin_N_[j];
				omega_small_ += bin_omega_[j];
				x_small_ += bin_x_[j];
				h_over_c_small_ += bin_omega_[j] * bin_h_over_c_[j];
				o_over_c_small_ += bin_omega_[j] * bin_o_over_c_[j];
				o_over_h_small_ += bin_omega_[j] * bin_o_over_h_[j];
			}

			double denominator = omega_small_;
			if (denominator < small_eps)	denominator = 1.e32;
			h_over_c_small_ /= denominator;
			o_over_c_small_ /= denominator;
			o_over_h_small_ /= denominator;

			// Analysis of soot
			{
				fv_large_ = 0.;
				rho_large_ = 0.;
				N_large_ = 0.;
				omega_large_ = 0.;
				x_large_ = 0.;
				h_over_c_large_ = 0.;
				o_over_c_large_ = 0.;
				o_over_h_large_ = 0.;

				for (unsigned int i = 0; i < bin_indices_large_.size(); i++)
				{
					int j = bin_indices_large_[i];
					fv_large_ += bin_fv_[j];
					rho_large_ += bin_rho_[j];
					N_large_ += bin_N_[j];
					omega_large_ += bin_omega_[j];
					x_large_ += bin_x_[j];
					h_over_c_large_ += bin_omega_[j] * bin_h_over_c_[j];
					o_over_c_large_ += bin_omega_[j] * bin_o_over_c_[j];
					o_over_h_large_ += bin_omega_[j] * bin_o_over_h_[j];
				}

				double denominator = omega_large_;
				if (denominator < small_eps)	denominator = 1.e32;
				h_over_c_large_ /= denominator;
				o_over_c_large_ /= denominator;
				o_over_h_large_ /= denominator;

				// Spherical particles
				fv_large_spherical_ = 0.;
				rho_large_spherical_ = 0.;
				N_large_spherical_ = 0.;
				omega_large_spherical_ = 0.;
				x_large_spherical_ = 0.;
				for (unsigned int i = 0; i < bin_indices_large_spherical_.size(); i++)
				{
					int j = bin_indices_large_spherical_[i];
					fv_large_spherical_ += bin_fv_[j];
					rho_large_spherical_ += bin_rho_[j];
					N_large_spherical_ += bin_N_[j];
					omega_large_spherical_ += bin_omega_[j];
					x_large_spherical_ += bin_x_[j];
				}

				// Aggregates
				fv_large_aggregates_ = 0.;
				rho_large_aggregates_ = 0.;
				N_large_aggregates_ = 0.;
				omega_large_aggregates_ = 0.;
				x_large_aggregates_ = 0.;
				for (unsigned int i = 0; i < bin_indices_large_aggregates_.size(); i++)
				{
					int j = bin_indices_large_aggregates_[i];
					fv_large_aggregates_ += bin_fv_[j];
					rho_large_aggregates_ += bin_rho_[j];
					N_large_aggregates_ += bin_N_[j];
					omega_large_aggregates_ += bin_omega_[j];
					x_large_aggregates_ += bin_x_[j];
				}
				
			}

			// Analysis of PAHs
			{
				omega_pah_1_2_rings_ = 0.;
				for (unsigned int i = 0; i < pah_1_2_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_1_2_rings_indices_global_[i];
					omega_pah_1_2_rings_ += omegaGas(j);
				}

				omega_pah_3_4_rings_ = 0.;
				for (unsigned int i = 0; i < pah_3_4_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_3_4_rings_indices_global_[i];
					omega_pah_3_4_rings_ += omegaGas(j);
				}

				omega_pah_more_than_4_rings_ = 0.;
				for (unsigned int i = 0; i < pah_more_than_4_rings_indices_global_.size(); i++)
				{
					const unsigned int j = pah_more_than_4_rings_indices_global_[i];
					omega_pah_more_than_4_rings_ += omegaGas(j);
				}
			}

		}
	}

	void PolimiSoot_Analyzer::Distribution()
	{
		//if (iBin_ == true)
		{
			for (unsigned int k = 0; k < bin_baskets_.size(); k++)
			{
				bin_baskets_N_[k] = 0.;
				bin_baskets_fv_[k] = 0.;
				bin_baskets_rho_[k] = 0.;
				bin_baskets_omega_[k] = 0.;
				bin_baskets_x_[k] = 0.;

				for (unsigned int i = 0; i < bin_baskets_indices_[k].size(); i++)
				{
					int j = bin_baskets_indices_[k][i];
					bin_baskets_N_[k] += bin_N_[j];
					bin_baskets_fv_[k] += bin_fv_[j];
					bin_baskets_rho_[k] += bin_rho_[j];
					bin_baskets_omega_[k] += bin_omega_[j];
					bin_baskets_x_[k] += bin_x_[j];
				}

				dN_over_dlog10d_[k] = bin_baskets_N_[k] / bin_baskets_dlog10d_[k];
				dN_over_dlog10V_[k] = bin_baskets_N_[k] / bin_baskets_dlog10V_[k];
				dN_over_dlog10m_[k] = bin_baskets_N_[k] / bin_baskets_dlog10m_[k];
			}
		}

		//if (iBin_ == true)
		{
			double small_eps = 1e-16;

			double m0 = 0.;
			double m1 = 0.;
			double m2 = 0.;
			double m3 = 0.;
			for (unsigned int k = 0; k < bin_baskets_.size(); k++)
			{
				m0 += bin_baskets_N_[k];
				m1 += bin_baskets_N_[k] * bin_baskets_d_[k];
				m2 += bin_baskets_N_[k] * bin_baskets_d_[k] * bin_baskets_d_[k];
				m3 += bin_baskets_N_[k] * bin_baskets_d_[k] * bin_baskets_d_[k] * bin_baskets_d_[k];
			}

			if (m0 < small_eps)	m0 = 1e32;
			if (m2 < small_eps)	m2 = 1e32;

			// Mean diameters [m]
			dmean_N_large_ = m1 / m0;
			d32_N_large_ = m3 / m2;

			// Variance [m2]
			dvariance_N_ = 0.;
			for (unsigned int k = 0; k < bin_baskets_.size(); k++)
				dvariance_N_ += bin_baskets_N_[k] * std::pow(bin_baskets_d_[k] - dmean_N_large_, 2.);
			dvariance_N_ /= m0;

			// Standard deviation [m]
			dstd_N_ = std::sqrt(dvariance_N_);
		}
	}

	void PolimiSoot_Analyzer::WriteDistribution(std::ofstream& fSootDistribution, const double t, const double x, const double y, const double z, const double T)
	{
		for (unsigned int k = 0; k < bin_baskets_.size(); k++)
		{
			fSootDistribution << std::scientific << std::setw(20) << std::left << t;
			fSootDistribution << std::scientific << std::setw(20) << std::left << x;
			fSootDistribution << std::scientific << std::setw(20) << std::left << y;
			fSootDistribution << std::scientific << std::setw(20) << std::left << z;
			fSootDistribution << std::scientific << std::setw(20) << std::left << T;

			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_d_[k] * 1e9;		// [nm]
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_m_[k] * 1e12;		// [ng]
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_V_[k] * 1e27;		// [nm3]

			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_d_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_m_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_V_[k]);

			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_dlog10d_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_dlog10m_[k]);
			fSootDistribution << std::scientific << std::setw(20) << std::left << std::log10(bin_baskets_dlog10V_[k]);

			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_fv_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_x_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_omega_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_rho_[k];
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_N_[k] / 1.e6; // [#/cm3]
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_N_[k] / bin_baskets_dlog10d_[k] / 1.e6; // [#/cm3]

			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_fv_[k] / (fv_large_ + fv_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_x_[k] / (x_large_ + x_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_omega_[k] / (omega_large_ + omega_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_rho_[k] / (rho_large_ + rho_small_);
			fSootDistribution << std::scientific << std::setw(20) << std::left << bin_baskets_N_[k] / (N_large_ + N_small_);

			fSootDistribution << std::endl;
		}
	}

	void PolimiSoot_Analyzer::WriteDistributionLabel(std::ofstream &fSoot)
	{
		fSoot	<< std::setw(20) << std::left << "time[s](1)"
			<< std::setw(20) << std::left << "x[m](2)"
			<< std::setw(20) << std::left << "y[m](3)"
			<< std::setw(20) << std::left << "z[m](4)"
			<< std::setw(20) << std::left << "T[K](5)";

		fSoot	<< std::setw(20) << std::left << "d[nm](6)"
				<< std::setw(20) << std::left << "m[ng](7)"
				<< std::setw(20) << std::left << "V[nm3](8)"
				<< std::setw(20) << std::left << "log10d(9)"
				<< std::setw(20) << std::left << "log10m(10)"
				<< std::setw(20) << std::left << "log10V(11)"
				<< std::setw(20) << std::left << "dlog10d(12)"
				<< std::setw(20) << std::left << "dlog10m(13)"
				<< std::setw(20) << std::left << "dlog10V(14)";

		fSoot	<< std::setw(20) << std::left << "fv(15)"
				<< std::setw(20) << std::left << "x(16)"
				<< std::setw(20) << std::left << "y(17)"
				<< std::setw(20) << std::left << "rho[kg/m3](18)"
				<< std::setw(20) << std::left << "N[#/cm3](19)"
				<< std::setw(20) << std::left << "dN/dlgd10[cm-3](20)";

		fSoot	<< std::setw(20) << std::left << "fvN(21)"
				<< std::setw(20) << std::left << "xN(22)"
				<< std::setw(20) << std::left << "yN(23)"
				<< std::setw(20) << std::left << "rhoN[](24)"
				<< std::setw(20) << std::left << "NN[](25)";

		fSoot	<< std::endl << std::endl;
	}

	void PolimiSoot_Analyzer::WriteBinData()
	{
		std::ofstream fOutput;
		fOutput.open("PolimiSootModel.out", std::ios::out);
		fOutput.setf(std::ios::scientific);

		fOutput << std::setw(5) << std::left << "#";
		fOutput << std::setw(5) << std::left << "#";
		fOutput << std::setw(10) << std::left << "Name";
		fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << "MW[kg/kmol]";	// [kg/kmol]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "d[nm]";			// [nm]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "dsph[nm]";		// [nm]
		fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << "dcol[nm]";		// [nm]
		fOutput << std::setw(16) << std::scientific << std::left << "m[ng]";							// [ng] 
		fOutput << std::setw(16) << std::scientific << std::left << "V[nm3]";							// [nm3]
		fOutput << std::setw(12) << std::scientific << std::left << "np[-]";							// [-]
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "C";	
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "H";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << "O";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "H/C";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "O/C";
		fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << "O/H";
		fOutput << std::setw(18) << std::scientific << std::left << "Teta(corr.fact.)";					// [-]
		fOutput << std::endl;

		for (unsigned int i = 0; i < bin_indices_small_.size(); i++)
		{
			int j = bin_indices_small_[i];
			fOutput << std::setw(5) << std::left << j + 1;
			fOutput << std::setw(5) << std::left << i + 1;
			fOutput << std::setw(10) << std::left << bin_names_[j].c_str();
			fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << bin_mw_[j];			// [kg/kmol]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_d_[j]  * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
			fOutput << std::setw(16) << std::scientific << std::left << bin_m_[j] * 1e12;						// [ng] 
			fOutput << std::setw(16) << std::scientific << std::left << bin_V_[j] * 1e27;						// [nm3]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(2) << bin_np_[j];			// [-]
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_h_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_o_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_h_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_h_[j];
			fOutput << std::setw(18) << std::scientific << std::left << bin_diffusivity_correction_factors_[j];
			fOutput << std::endl;
		}
		fOutput << std::endl;

		for (unsigned int i = 0; i < bin_indices_large_.size(); i++)
		{
			int j = bin_indices_large_[i];
			fOutput << std::setw(5) << std::left << j + 1;
			fOutput << std::setw(5) << std::left << i + 1;
			fOutput << std::setw(10) << std::left << bin_names_[j].c_str();
			fOutput << std::setw(16) << std::fixed << std::left << std::setprecision(2) << bin_mw_[j];			// [kg/kmol]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_d_[j] * 1e9;		// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_ds_[j] * 1e9;	// [nm]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(5) << bin_dc_[j] * 1e9;	// [nm]
			fOutput << std::setw(16) << std::scientific << std::left << bin_m_[j] * 1e12;						// [ng] 
			fOutput << std::setw(16) << std::scientific << std::left << bin_V_[j] * 1e27;						// [nm3]
			fOutput << std::setw(12) << std::fixed << std::left << std::setprecision(2) << bin_np_[j];			// [-]
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_h_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(0) << bin_o_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_h_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_c_[j];
			fOutput << std::setw(10) << std::fixed << std::left << std::setprecision(4) << bin_o_over_h_[j];
			fOutput << std::setw(18) << std::scientific << std::left << bin_diffusivity_correction_factors_[j];
			fOutput << std::endl;
		}

		fOutput.close();
	}

	void PolimiSoot_Analyzer::WriteLabelIntegralDataFile(std::ofstream& fOutput, unsigned int& count, const unsigned int width)
	{
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "fv(L)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "x(L)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "w(L)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "rho(L)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "N(L)[#/cm3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "H/C(L)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "O/C(L)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "O/H(L)[-]", count);

		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "fv(S)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "x(S)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "w(S)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "rho(S)[kg/m3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "N(S)[#/cm3]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "H/C(S)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "O/C(S)[-]", count);
		OpenSMOKE::PrintTagOnASCIILabel(width, fOutput, "O/H(S)[-]", count);
	}

	void PolimiSoot_Analyzer::WriteIntegralDataFile(std::ofstream& fOutput, const unsigned int width)
	{
		// Large sections (soot)
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << fv_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << x_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << omega_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << rho_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << N_large() / 1.e6;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << h_over_c_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << o_over_c_large();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << o_over_h_large();

		// Small sections (PAH)
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << fv_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << x_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << omega_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << rho_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << N_small() / 1.e6;
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << h_over_c_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << o_over_c_small();
		fOutput << std::scientific << std::setprecision(9) << std::setw(width) << o_over_h_small();
	}
}

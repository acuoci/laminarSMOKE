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

#include "math/OpenSMOKEUtilities.h"
#include "TransportPropertiesMap_CHEMKIN.h"
#include "preprocessing/CollisionIntegralMatrices.hpp"

#if OPENSMOKE_USE_MKL == 1
#include "mkl.h"
#endif

namespace OpenSMOKE
{
	const double TransportPropertiesMap_CHEMKIN::threshold_ = 1.e-14;

	TransportPropertiesMap_CHEMKIN::TransportPropertiesMap_CHEMKIN(const unsigned int nSpecies)
	{
		this->nspecies_ = nSpecies;
		temperature_lambda_must_be_recalculated_ = true;
		temperature_eta_must_be_recalculated_ = true;
		temperature_gamma_must_be_recalculated_ = true;
		temperature_teta_must_be_recalculated_ = true;
		pressure_gamma_must_be_recalculated_ = true;

		this->T_ = this->P_ = 0.;
		this->T_old_ = this->P_old_ = 0.;

		this->species_bundling_ = false;
		is_lennard_jones_available_ = false;

		MemoryAllocation();
	}

	TransportPropertiesMap_CHEMKIN::TransportPropertiesMap_CHEMKIN(rapidxml::xml_document<>& doc)
	{
		temperature_lambda_must_be_recalculated_ = true;
		temperature_eta_must_be_recalculated_ = true;
		temperature_gamma_must_be_recalculated_ = true;
		temperature_teta_must_be_recalculated_ = true;
		pressure_gamma_must_be_recalculated_ = true;
		this->T_ = this->P_ = 0.;
		this->T_old_ = this->P_old_ = 0.;

		this->species_bundling_ = false;

		ImportSpeciesFromXMLFile(doc);
		ImportViscosityModelFromXMLFile(doc);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(doc);
	}

	TransportPropertiesMap_CHEMKIN::TransportPropertiesMap_CHEMKIN(const TransportPropertiesMap_CHEMKIN& rhs)
	{
		CopyFromMap(rhs);
	}

	TransportPropertiesMap_CHEMKIN::~TransportPropertiesMap_CHEMKIN()
	{
		delete[] this->M;
		delete[] this->fittingLambda;
		delete[] this->fittingEta;
		delete[] this->fittingTeta;
		delete[] this->fittingGamma;

		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			delete[] MWRatio1over4;
			delete[] phi_eta_sup;
			delete[] phi_eta_inf;
			delete[] sqrtEta;
			delete[] usqrtEta;
		}
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			delete[] sqrtMWRatio_inf;
			delete[] sqrtMWRatio_sup;
			delete[] sqrtMW;
		}

		delete[] sum_diffusion_coefficients;
		delete[] x_corrected;

		if (this->species_bundling_ == true)
		{
			delete[] this->bundling_fittingGammaSelfDiffusion_;
			delete[] this->bundling_fittingGamma_;
			delete[] this->bundling_gammaSpecies_;
			delete[] this->bundling_gammaSpeciesSelfDiffusion_;
		}

		iThermalDiffusionRatios_.clear();
	}

	void TransportPropertiesMap_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_old_ = this->T_;
		this->T_ = T;

		const double epsilon = 0.;
		if (std::fabs(this->T_ - this->T_old_) / this->T_ > epsilon)
		{
			temperature_lambda_must_be_recalculated_ = true;
			temperature_eta_must_be_recalculated_ = true;
			temperature_gamma_must_be_recalculated_ = true;
			temperature_teta_must_be_recalculated_ = true;
		}
	}

	void TransportPropertiesMap_CHEMKIN::SetPressure(const double& P)
	{
		this->P_old_ = this->P_;
		this->P_ = P;

		if (std::fabs(this->P_ - this->P_old_) / this->P_ > 1.e-14)
		{
			pressure_gamma_must_be_recalculated_ = true;
		}
	}

	void TransportPropertiesMap_CHEMKIN::CopyFromMap(const TransportPropertiesMap_CHEMKIN& rhs)
	{
		this->nspecies_ = rhs.nspecies_;
		this->species_bundling_ = rhs.species_bundling_;

		MemoryAllocation();

		this->count_species_thermal_diffusion_ratios_ = rhs.count_species_thermal_diffusion_ratios_;
		this->sumK = rhs.sumK;
		this->iThermalDiffusionRatios_ = rhs.iThermalDiffusionRatios_;

		for (unsigned int i = 0; i < this->nspecies_; i++)
			this->M[i] = rhs.M[i];

		for (unsigned int i = 0; i < this->nspecies_ * 4; i++)
			this->fittingLambda[i] = rhs.fittingLambda[i];

		for (unsigned int i = 0; i < this->nspecies_ * 4; i++)
			this->fittingEta[i] = rhs.fittingEta[i];

		for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2 * 4; i++)
			this->fittingGamma[i] = rhs.fittingGamma[i];

		this->fittingTeta = new double[4 * count_species_thermal_diffusion_ratios_*this->nspecies_];
		for (unsigned int i = 0; i < this->count_species_thermal_diffusion_ratios_*(this->nspecies_) * 4; i++)
			this->fittingTeta[i] = rhs.fittingTeta[i];

		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->MWRatio1over4[i] = rhs.MWRatio1over4[i];

			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->phi_eta_sup[i] = rhs.phi_eta_sup[i];

			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->phi_eta_inf[i] = rhs.phi_eta_inf[i];

			for (unsigned int i = 0; i < this->nspecies_; i++)
				this->sqrtEta[i] = rhs.sqrtEta[i];

			for (unsigned int i = 0; i < this->nspecies_; i++)
				this->usqrtEta[i] = rhs.usqrtEta[i];
		}
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->sqrtMWRatio_inf[i] = rhs.sqrtMWRatio_inf[i];

			for (unsigned int i = 0; i < this->nspecies_*(this->nspecies_ - 1) / 2; i++)
				this->sqrtMWRatio_sup[i] = rhs.sqrtMWRatio_sup[i];

			for (unsigned int i = 0; i < this->nspecies_; i++)
				this->sqrtMW[i] = rhs.sqrtMW[i];
		}

		if (this->species_bundling_ == true)
		{
			this->bundling_number_groups_ = rhs.bundling_number_groups_;
			this->bundling_reference_species_ = rhs.bundling_reference_species_;
			this->bundling_groups_ = rhs.bundling_groups_;
			this->bundling_species_group_ = rhs.bundling_species_group_;

			this->bundling_fittingGammaSelfDiffusion_ = new double[4 * this->nspecies_];
			this->bundling_fittingGamma_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2 * 4];
			this->bundling_gammaSpecies_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2];
			this->bundling_gammaSpeciesSelfDiffusion_ = new double[this->bundling_number_groups_];

			for (unsigned int i = 0; i < 4 * this->nspecies_; i++)
				this->bundling_fittingGammaSelfDiffusion_[i] = rhs.bundling_fittingGammaSelfDiffusion_[i];

			for (unsigned int i = 0; i < this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2 * 4; i++)
				this->bundling_fittingGamma_[i] = rhs.bundling_fittingGamma_[i];

			for (unsigned int i = 0; i < this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2; i++)
				this->bundling_gammaSpecies_[i] = rhs.bundling_gammaSpecies_[i];

			for (unsigned int i = 0; i < this->bundling_number_groups_; i++)
				this->bundling_gammaSpeciesSelfDiffusion_[i] = rhs.bundling_gammaSpeciesSelfDiffusion_[i];
		}

		this->tetaSpecies_.resize(this->nspecies_*iThermalDiffusionRatios_.size());
		this->tetaSpecies_.setZero();
	}

	void TransportPropertiesMap_CHEMKIN::MemoryAllocation()
	{
		// viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING;
		// viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_MATHUR_SAXENA;
		// viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE; (default)

		sum_threshold_ = 1. + threshold_*double(this->nspecies_);

		M = new double[this->nspecies_];
		fittingLambda = new double[this->nspecies_ * 4];
		fittingEta = new double[this->nspecies_ * 4];
		fittingGamma = new double[this->nspecies_*(this->nspecies_ - 1) / 2 * 4];

		// fittingTeta is not inizialized here, because the number of species for which the thermal diffusion ratios is 
		// important is not yet known

		sumK.resize(this->nspecies_);
		sumK.setZero();

		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			MWRatio1over4 = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			phi_eta_sup = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			phi_eta_inf = new double[this->nspecies_*(this->nspecies_ - 1) / 2];

			sqrtEta = new double[this->nspecies_];
			usqrtEta = new double[this->nspecies_];
		}
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			sqrtMWRatio_inf = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			sqrtMWRatio_sup = new double[this->nspecies_*(this->nspecies_ - 1) / 2];
			sqrtMW = new double[this->nspecies_];
		}

		sum_diffusion_coefficients = new double[this->nspecies_];
		x_corrected = new double[this->nspecies_];

		this->lambdaSpecies_.resize(this->nspecies_);
		this->lambdaSpecies_.setZero();

		this->etaSpecies_.resize(this->nspecies_);
		this->etaSpecies_.setZero();

		this->gammaSpecies_.resize(this->nspecies_*(this->nspecies_ - 1) / 2);
		this->gammaSpecies_.setZero();

		// Thermal diffusion ratios are allocated after reading the properties
	}

	void TransportPropertiesMap_CHEMKIN::SetCoefficients(const unsigned int j, const double* coefficients)
	{
		// Molecular weight
		M[j] = coefficients[0];

		// Thermal conductivity
		{
			unsigned int i = j * 4;
			fittingLambda[i++] = coefficients[1];
			fittingLambda[i++] = coefficients[2];
			fittingLambda[i++] = coefficients[3];
			fittingLambda[i++] = coefficients[4];
		}

		// Dynamic viscosity
		{
			unsigned int i = j * 4;
			fittingEta[i++] = coefficients[5];
			fittingEta[i++] = coefficients[6];
			fittingEta[i++] = coefficients[7];
			fittingEta[i++] = coefficients[8];
		}

		// Mass diffusion coefficients
		{
			unsigned int countFitting = 0;
			for (unsigned int i = 1; i <= j; i++)
				countFitting += this->nspecies_ - i;
			countFitting *= 4;
			unsigned int countCoefficients = 9 + 4 * (j + 1);
			for (unsigned int k = j + 1; k < this->nspecies_; k++)
			for (unsigned int i = 1; i <= 4; i++)
				fittingGamma[countFitting++] = coefficients[countCoefficients++];
		}

		// Thermal diffusion coefficients
		{
			if (M[j] <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)
				iThermalDiffusionRatios_.push_back(j + 1);

			unsigned int i = j * 4 * this->nspecies_;
			unsigned int count = 8 + 4 * this->nspecies_ + 1;

			for (unsigned int k = 0; k < this->nspecies_; k++)
			for (unsigned int m = 1; m <= 4; m++)
				fittingTeta[i++] = coefficients[count++];
		}

		if (j == this->nspecies_ - 1)
			CompleteInitialization();
	}

	void TransportPropertiesMap_CHEMKIN::ImportCoefficientsFromASCIIFile(std::istream& fInput)
	{
		{
			fInput >> count_species_thermal_diffusion_ratios_;
			fittingTeta = new double[4 * count_species_thermal_diffusion_ratios_*this->nspecies_];
		}

		std::string tag;
		for (unsigned int j = 0; j < this->nspecies_; j++)
		{
			// Molecular weight
			fInput >> M[j];

			// Thermal conductivity
			{
				unsigned int i = j * 4;
				fInput >> fittingLambda[i++];
				fInput >> fittingLambda[i++];
				fInput >> fittingLambda[i++];
				fInput >> fittingLambda[i++];
			}

			// Dynamic viscosity
			{
				unsigned int i = j * 4;
				fInput >> fittingEta[i++];
				fInput >> fittingEta[i++];
				fInput >> fittingEta[i++];
				fInput >> fittingEta[i++];
			}

			// Mass diffusion coefficients
			{
				unsigned int countFitting = 0;
				for (unsigned int i = 1; i <= j; i++)
					countFitting += this->nspecies_ - i;
				countFitting *= 4;
				for (unsigned int k = j + 1; k < this->nspecies_; k++)
				for (unsigned int i = 1; i <= 4; i++)
					fInput >> fittingGamma[countFitting++];
			}

			// Thermal diffusion coefficients
			{
				if (M[j] <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)
				{
					iThermalDiffusionRatios_.push_back(j + 1);

					std::size_t i = 4 * this->nspecies_*(iThermalDiffusionRatios_.size() - 1);

					for (unsigned int k = 0; k < this->nspecies_; k++)
					for (unsigned int m = 1; m <= 4; m++)
						fInput >> fittingTeta[i++];
				}
			}
		}

		// Check that everything is read properly
		fInput >> tag;
		CheckForFatalError(tag == "E");

		// Complete the initialization
		CompleteInitialization();
	}

	void TransportPropertiesMap_CHEMKIN::ImportLennardJonesCoefficientsFromASCIIFile(std::istream& fInput)
	{
		is_lennard_jones_available_ = true;
		mu_.resize(this->nspecies_);
		sigma_.resize(this->nspecies_);
		epsilon_over_kb_.resize(this->nspecies_);

		for (unsigned int j = 0; j < this->nspecies_; j++)
		{
			fInput >> mu_[j];					// [kg]
			fInput >> sigma_[j];				// [m]
			fInput >> epsilon_over_kb_[j];		// [K]
		}
	}

	void TransportPropertiesMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		std::cout << " * Reading transport properties from XML file..." << std::endl;

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

		// Transport coefficients (from fitting)
		rapidxml::xml_node<>* transport_node = opensmoke_node->first_node("Transport");
		if (transport_node != 0)
		{
			std::string transport_type = transport_node->first_attribute("type")->value();
			if (transport_type == "CHEMKIN")
			{
				std::stringstream fInput;
				fInput << transport_node->value();

				ImportCoefficientsFromASCIIFile(fInput);
			}
			else
				ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportCoefficientsFromXMLFile", "Transport type not supported!");
		}
		else
			ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportCoefficientsFromXMLFile", "Transport tag was not found!");

		// Lennard-Jones parametes (not stricly needed)
		rapidxml::xml_node<>* lennard_jones_node = opensmoke_node->first_node("Lennard-Jones");
		if (lennard_jones_node != 0)
		{
			std::stringstream fInput;
			fInput << lennard_jones_node->value();
			ImportLennardJonesCoefficientsFromASCIIFile(fInput);
		}
		else
			ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportLennardJonesCoefficientsFromASCIIFile", "Lennard-Jones coefficients unavailable!");
	}

	unsigned int index_from_mtrix_to_vector(const unsigned int n, const unsigned int i, const unsigned int j)
	{
		unsigned int index = 0;
		for (unsigned int k = 1; k <= i; k++)
			index += (n - k);

		index += j - i - 1;

		return index;
	}

	void TransportPropertiesMap_CHEMKIN::ImportSpeciesBundlingFromXMLFile(rapidxml::xml_document<>& doc, const double epsilon)
	{
		std::cout << " * Reading species bundling from XML file..." << std::endl;

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

		rapidxml::xml_node<>* transport_node = opensmoke_node->first_node("SpeciesBundling");
		if (transport_node != 0)
		{
			// Read self diffusion coefficients
			{
				rapidxml::xml_node<>* selfdiffusion_node = transport_node->first_node("SelfDiffusion");
			
				std::stringstream fInput;
				fInput << selfdiffusion_node->value();

				this->bundling_fittingGammaSelfDiffusion_ = new double[4 * this->nspecies_];

				for (unsigned int j = 0; j < this->nspecies_; j++)
				{
					{
						unsigned int i = j * 4;
						fInput >> this->bundling_fittingGammaSelfDiffusion_[i++];
						fInput >> this->bundling_fittingGammaSelfDiffusion_[i++];
						fInput >> this->bundling_fittingGammaSelfDiffusion_[i++];
						fInput >> this->bundling_fittingGammaSelfDiffusion_[i++];
					}
				}
			}

			bool iFound = false;
			for (rapidxml::xml_node<> *bundling_node = transport_node->first_node("Bundling"); bundling_node; bundling_node = bundling_node->next_sibling())
			{
				std::string epsilon_string = bundling_node->first_attribute("epsilon")->value();
				double epsilon_double = std::atof(epsilon_string.c_str());

				if (iFound == false && epsilon_double == epsilon)
				{
					this->species_bundling_ = true;

					iFound = true;
					std::cout << "   Maximum error (epsilon): " << epsilon_double << std::endl;

					// Number of groups
					{
						rapidxml::xml_node<>* number_of_groups_node = bundling_node->first_node("NumberOfGroups");
						this->bundling_number_groups_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_groups_node->value())));
					}

					// Reference species
					{
						rapidxml::xml_node<>* current_node = bundling_node->first_node("ReferenceSpecies");
						std::stringstream fInput;
						fInput << current_node->value();

						unsigned int n;
						fInput >> n;
						this->bundling_reference_species_.resize(n);
						for (unsigned int i = 0; i<n; i++)
							fInput >> this->bundling_reference_species_[i];
					}

					// Species in groups
					{
						rapidxml::xml_node<>* current_node = bundling_node->first_node("SpeciesInGroups");
						std::stringstream fInput;
						fInput << current_node->value();

						this->bundling_groups_.resize(this->bundling_number_groups_);
						for (unsigned k = 0; k < this->bundling_number_groups_; k++)
						{
							unsigned int n;
							fInput >> n;
							this->bundling_groups_[k].resize(n);
							for (unsigned int i = 0; i<n; i++)
								fInput >> this->bundling_groups_[k][i];
						}
					}

					// Group of species
					{
						rapidxml::xml_node<>* current_node = bundling_node->first_node("GroupOfSpecies");
						std::stringstream fInput;
						fInput << current_node->value();

						unsigned int n;
						fInput >> n;
						this->bundling_species_group_.resize(n);
						for (unsigned int i = 0; i<n; i++)
							fInput >> this->bundling_species_group_[i];
					}

					std::cout << "   Number of groups: " << this->bundling_number_groups_ << "/" << this->nspecies_ << std::endl;
					std::cout << "   List of reference species: " << std::endl;
					for (unsigned int i = 0; i < this->bundling_number_groups_; i++)
						std::cout << "    - " << this->bundling_reference_species_[i] << " (" << this->bundling_groups_[i].size() << ")" << std::endl;

					// Diffusion coefficients
					bundling_sum_diffusion_coefficients_.resize(this->bundling_number_groups_);
					bundling_sum_x_groups_.resize(this->bundling_number_groups_);

					this->bundling_fittingGamma_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2 * 4];
					this->bundling_gammaSpecies_ = new double[this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2];
					this->bundling_gammaSpeciesSelfDiffusion_ = new double[this->bundling_number_groups_];

					unsigned int countFitting = 0;
					for (unsigned int i = 0; i < this->bundling_number_groups_; i++)
					for (unsigned int j = i + 1; j < this->bundling_number_groups_; j++)
						{
							unsigned int index_i = this->bundling_reference_species_[i];
							unsigned int index_j = this->bundling_reference_species_[j];
							if (this->bundling_reference_species_[i] > this->bundling_reference_species_[j])
							{
								index_j = this->bundling_reference_species_[i];
								index_i = this->bundling_reference_species_[j];
							}

							unsigned int index = index_from_mtrix_to_vector(this->nspecies_, index_i, index_j) * 4;

							for (unsigned int i = 0; i < 4; i++)
								this->bundling_fittingGamma_[countFitting + i] = fittingGamma[index + i];
							countFitting += 4;

						}
				}
			}

			if (iFound == false)
				ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportSpeciesBundlingFromXMLFile", "The requested epsilon for species bundling was not found");
		}
		else
			ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportSpeciesBundlingFromXMLFile", "The kinetic mechanism was preprocessed without enabling species bundling");
	}

	void TransportPropertiesMap_CHEMKIN::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		rapidxml::xml_node<>* names_of_species_node = opensmoke_node->first_node("NamesOfSpecies");
		try
		{
			// Recognize number of species
			this->nspecies_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));
			
			// Recognize indices of relevant species
			std::vector<std::string> names_;
			names_.resize(this->nspecies_);
			std::stringstream names_of_species_string;
			names_of_species_string << names_of_species_node->value();
			for (unsigned int i = 0; i<this->nspecies_; i++)
				names_of_species_string >> names_[i];

			index_H2O_ = index_CO2_ = index_CO_ = index_CH4_ = 0;
			for (unsigned int i = 0; i < this->nspecies_; i++)
			{
				if (names_[i] == "H2O" || names_[i] == "h2o")	index_H2O_ = i + 1;
				if (names_[i] == "CO2" || names_[i] == "co2")	index_CO2_ = i + 1;
				if (names_[i] == "CO"  || names_[i] == "co")	index_CO_  = i + 1;
				if (names_[i] == "CH4" || names_[i] == "ch4")	index_CH4_ = i + 1;
			}
		}
		catch (...)
		{
			ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the list of species.");
		}
	}

	void TransportPropertiesMap_CHEMKIN::ImportViscosityModelFromXMLFile(rapidxml::xml_document<>& doc)
	{
		viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE;

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* viscosity_model_node = opensmoke_node->first_node("ViscosityModel");
		if (viscosity_model_node != 0)
		{
			const std::string dummy = boost::trim_copy(std::string(viscosity_model_node->value()));

			std::cout << " * User defined viscosity model: " << dummy << std::endl;

			if (dummy == "Wilke")				viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE;
			else if (dummy == "Herning")		viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING;
			else if (dummy == "MathurSaxena")	viscosity_model = PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_MATHUR_SAXENA;
			else ErrorMessage("TransportPropertiesMap_CHEMKIN::ImportViscosityModelFromXMLFile", "Error in reading the viscosity model.");
		}
	}

	void TransportPropertiesMap_CHEMKIN::CompleteInitialization()
	{
		if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			unsigned int i=0;
			for(unsigned int k=0;k<this->nspecies_;k++)
				for(unsigned int j=k+1;j<this->nspecies_;j++)
				{
					phi_eta_sup[i] = 1./std::sqrt( 8.*(1.+M[k]/M[j]) );
					phi_eta_inf[i] = 1./std::sqrt( 8.*(1.+M[j]/M[k]) );
					i++;
				}

			i=0;
			for(unsigned int k=0;k<this->nspecies_;k++)
				for(unsigned int j=k+1;j<this->nspecies_;j++)
					MWRatio1over4[i++] = std::sqrt(std::sqrt(M[j]/M[k]));
		}

		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			for(unsigned int k=0;k<this->nspecies_;k++)
				sqrtMW[k] = std::sqrt(M[k]);

			unsigned int i=0;
			for(unsigned int k=0;k<this->nspecies_;k++)
				for(unsigned int j=k+1;j<this->nspecies_;j++)
				{
					sqrtMWRatio_sup[i] = std::sqrt(M[j]/M[k]);
					sqrtMWRatio_inf[i++] = std::sqrt(M[k]/M[j]);
				}
		}
		
		this->tetaSpecies_.resize(this->nspecies_*boost::lexical_cast<unsigned int>(iThermalDiffusionRatios_.size()));
		this->tetaSpecies_.setZero();
	}

	inline void TransportPropertiesMap_CHEMKIN::lambda()
	{
        if (temperature_lambda_must_be_recalculated_ == true)
        {
			const double logT=std::log(T_);
			const double logT2=logT*logT;
			const double logT3=logT*logT2;

			#if OPENSMOKE_USE_MKL == 0

				const double* k = fittingLambda;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
							sum+= (*k++)*logT;
							sum+= (*k++)*logT2;
							sum+= (*k++)*logT3;
					this->lambdaSpecies_(j)=std::exp(sum);
				}

			#elif OPENSMOKE_USE_MKL == 1
	
				const double* k = fittingLambda;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
							sum+= (*k++)*logT;
							sum+= (*k++)*logT2;
							sum+= (*k++)*logT3;
					this->lambdaSpecies_(j)=sum;
				}
				vdExp( this->nspecies_, this->lambdaSpecies_.data(), this->lambdaSpecies_.data() );
	
			#endif

            temperature_lambda_must_be_recalculated_ = false;
        }
	}

	inline void TransportPropertiesMap_CHEMKIN::eta()
	{
        if (temperature_eta_must_be_recalculated_ == true)
        {
			const double logT=std::log(T_);
			const double logT2=logT*logT;
			const double logT3=logT*logT2;

			#if OPENSMOKE_USE_MKL == 0

				const double* k = fittingEta;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
						   sum+= (*k++)*logT;
						   sum+= (*k++)*logT2;
						   sum+= (*k++)*logT3;
					this->etaSpecies_(j)=std::exp(sum);
				}

			#elif OPENSMOKE_USE_MKL == 1
	
				const double* k = fittingEta;
				for (unsigned int j=0;j<this->nspecies_;j++)
				{
					double sum = (*k++);
						   sum+= (*k++)*logT;
						   sum+= (*k++)*logT2;
						   sum+= (*k++)*logT3;
					this->etaSpecies_(j)=sum;
				}
				vdExp( this->nspecies_, this->etaSpecies_.data(), this->etaSpecies_.data() );
	
			#endif

			temperature_eta_must_be_recalculated_ = false;
		}
	}
	
	inline void TransportPropertiesMap_CHEMKIN::teta()
	{
		if (temperature_teta_must_be_recalculated_ == true)
		{
			const double T2=T_*T_;
			const double T3=T_*T2;
		
			// Thermal diffusion ratios evaluation
			for (unsigned int i=1;i<=iThermalDiffusionRatios_.size();i++)
			{
				unsigned int j = 4*(i-1)*this->nspecies_;
				unsigned int m = (i-1)*this->nspecies_;

				double* ptFitting = &fittingTeta[j];

				for (unsigned int k=1;k<=this->nspecies_;k++)
				{
					double sum  = (*ptFitting++);
						   sum += (*ptFitting++)*T_;
						   sum += (*ptFitting++)*T2;
						   sum += (*ptFitting++)*T3;
					this->tetaSpecies_[m++] = sum; 

				}
			}
                
			temperature_teta_must_be_recalculated_ = false;
		}
	}

	inline void TransportPropertiesMap_CHEMKIN::gamma()
	{
		if ( temperature_gamma_must_be_recalculated_ == false && pressure_gamma_must_be_recalculated_  == true )
		{
            const double multiplier = P_/P_old_;
			double* D = this->gammaSpecies_.data();
			for (unsigned int j=0;j<this->nspecies_;j++)
			for (unsigned int k=j+1;k<this->nspecies_;k++)
			{
				*D++ *= multiplier;
			}
			
			pressure_gamma_must_be_recalculated_ = false;
		}
            
        if ( temperature_gamma_must_be_recalculated_ == true || pressure_gamma_must_be_recalculated_ == true )
        {
			// Only the upper hals of this matrix is evaluated (the main diagonalis not evaluated)
			// Indeed the matrix is symmetric and the mixture rule is able to exploit this kind of symmetry

			const double P_bar = P_/100000.;
			const double logT=std::log(T_);
			const double logT2=logT*logT;
			const double logT3=logT*logT2;

			#if OPENSMOKE_USE_MKL == 0

			const double* d = fittingGamma;
			double* D = this->gammaSpecies_.data();
			for (unsigned int j=1;j<=this->nspecies_;j++)
				for (unsigned int k=j+1;k<=this->nspecies_;k++)
				{
					double sum = (*d++);
							sum+= (*d++)*logT;
							sum+= (*d++)*logT2;
							sum+= (*d++)*logT3;
					*D++ = P_bar/std::exp(sum);
				}

			#elif OPENSMOKE_USE_MKL == 1
	
			const double lnP_bar = std::log(P_bar);
			const double* d = fittingGamma;
			double* D = this->gammaSpecies_.data();
			for (unsigned int j=0;j<this->nspecies_;j++)
				for (unsigned int k=j+1;k<this->nspecies_;k++)
				{
					double sum = (*d++);
							sum+= (*d++)*logT;
							sum+= (*d++)*logT2;
							sum+= (*d++)*logT3;
					*D++ = lnP_bar-sum;
				}

			vdExp( this->nspecies_*(this->nspecies_-1)/2, this->gammaSpecies_.data(), this->gammaSpecies_.data() );

			#endif

            temperature_gamma_must_be_recalculated_ = false;
            pressure_gamma_must_be_recalculated_ = false;
        }
	}

	inline void TransportPropertiesMap_CHEMKIN::bundling_gamma()
	{
		if (temperature_gamma_must_be_recalculated_ == false && pressure_gamma_must_be_recalculated_ == true )
		{
			const double multiplier = P_ / P_old_;

			double* D = this->bundling_gammaSpecies_;
			for (unsigned int j = 0; j<this->bundling_number_groups_; j++)
			for (unsigned int k = j + 1; k<this->bundling_number_groups_; k++)
			{
				*D++ *= multiplier;
			}

			double* Dself = this->bundling_gammaSpeciesSelfDiffusion_;
			for (unsigned int j = 0; j < this->bundling_number_groups_; j++)
			{
				*Dself++ *= multiplier;
			}

			pressure_gamma_must_be_recalculated_ = false;
		}

		if ( temperature_gamma_must_be_recalculated_ == true || pressure_gamma_must_be_recalculated_ == true)
		{
			// Only the upper hals of this matrix is evaluated (the main diagonalis not evaluated)
			// Indeed the matrix is symmetric and the mixture rule is able to exploit this kind of symmetry

			const double P_bar = P_ / 100000.;
			const double logT = std::log(T_);
			const double logT2 = logT*logT;
			const double logT3 = logT*logT2;

			// Self diffusion coefficients
			double* Dself = this->bundling_gammaSpeciesSelfDiffusion_;
			for (unsigned int j = 0; j <this->bundling_number_groups_; j++)
			{
				unsigned int index = this->bundling_reference_species_[j];
				const double* d = &this->bundling_fittingGammaSelfDiffusion_[index * 4];

				double sum = (*d++);
				sum += (*d++)*logT;
				sum += (*d++)*logT2;
				sum += (*d++)*logT3;

				*Dself++ = P_bar / std::exp(sum);
			}


			#if OPENSMOKE_USE_MKL == 0

			const double* d = this->bundling_fittingGamma_;
			double* D = this->bundling_gammaSpecies_;
			for (unsigned int j = 1; j <= this->bundling_number_groups_; j++)
			for (unsigned int k = j + 1; k <= this->bundling_number_groups_; k++)
			{
				double sum = (*d++);
				sum += (*d++)*logT;
				sum += (*d++)*logT2;
				sum += (*d++)*logT3;
				*D++ = P_bar / std::exp(sum);
			}

			#elif OPENSMOKE_USE_MKL == 1

			const double lnP_bar = std::log(P_bar);
			const double* d = this->bundling_fittingGamma_;
			double* D = this->bundling_gammaSpecies_;
			for (unsigned int j = 0; j<this->bundling_number_groups_; j++)
			for (unsigned int k = j + 1; k<this->bundling_number_groups_; k++)
			{
				double sum = (*d++);
				sum += (*d++)*logT;
				sum += (*d++)*logT2;
				sum += (*d++)*logT3;
				*D++ = lnP_bar - sum;
			}

			vdExp(this->bundling_number_groups_*(this->bundling_number_groups_ - 1) / 2, this->bundling_gammaSpecies_, this->bundling_gammaSpecies_);

			#endif

			temperature_gamma_must_be_recalculated_ = false;
			pressure_gamma_must_be_recalculated_ = false;
		}
	}

	double TransportPropertiesMap_CHEMKIN::lambdaMix(const double* moleFractions)
	{
		// Calcolo della conducibilita della miscela
		// Formula di Mathur, Todor, Saxena - Molecular Physics 52:569 (1967)

		const double sum1 = Dot(this->nspecies_, moleFractions, this->lambdaSpecies_.data());
		const double sum2 = UDot(this->nspecies_, moleFractions, this->lambdaSpecies_.data());
		const double lambdamix = 0.50 * ( sum1 + 1./sum2 );
		
		return lambdamix;
	}

	double TransportPropertiesMap_CHEMKIN::etaMix(const double* moleFractions)
	{
		for (unsigned int k = 0; k < this->nspecies_; k++)
			sumK(k) = moleFractions[k];

		// Wilke - Journal of Chemical Physics 18:517 (1950)
		// Modified by Bird, Stewart, Lightfoot - Transport phenomena (1960)
		// Available in: Reid, Prausnitz, Poling - The properties of gases and liquids, p. 407
		if(viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE)
		{
			#if OPENSMOKE_USE_MKL == 0
			
				for(unsigned int k=0;k<this->nspecies_;k++)
				{
					sqrtEta[k] = std::sqrt(this->etaSpecies_(k));
					usqrtEta[k] = 1./sqrtEta[k];
				}

			#elif OPENSMOKE_USE_MKL == 1
			
				vdSqrt(this->nspecies_, this->etaSpecies_.data(), sqrtEta);
				vdInv(this->nspecies_, sqrtEta, usqrtEta);
			
			#endif

			const double* ptMWRatio1over4=MWRatio1over4;
			const double* ptphi_eta_sup=phi_eta_sup;
			const double* ptphi_eta_inf=phi_eta_inf;

			for (unsigned int k = 0; k < this->nspecies_; k++)
				for (unsigned int j = k + 1; j < this->nspecies_; j++)
				{
					double delta_phi = sqrtEta[k] * usqrtEta[j] * (*ptMWRatio1over4++);	// F.(49)
					sumK(k) += moleFractions[j] * (*ptphi_eta_sup++)*(1. + delta_phi)*(1. + delta_phi);
					sumK(j) += moleFractions[k] * (*ptphi_eta_inf++)*(1. + 1. / delta_phi)*(1. + 1. / delta_phi);
				}

			double etamix = 0.;
			for(unsigned int k=0;k<this->nspecies_;k++)
				etamix += moleFractions[k]*this->etaSpecies_(k)/sumK(k);				// F.(48)

			return etamix;
		}

		// Herning and Zipperer
		// Available in: Reid, Prausnitz, Poling - The properties of gases and liquids, p. 410
		// This is a simplified formula, less accurate, but much faster
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING)
		{
			double sum = 0;
			for (unsigned int k = 0; k < this->nspecies_; k++)
				sum += moleFractions[k] * sqrtMW[k];

			double etamix = 0.;
			for (unsigned int k = 0; k < this->nspecies_; k++)
				etamix += moleFractions[k] * this->etaSpecies_(k) * sqrtMW[k];				// F.(48)
			etamix /= sum;

			return etamix;
		}
		// Mathur and Saxena
		// Molecular Physics 52:569 (1967)
		else if (viscosity_model == PhysicalConstants::OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_MATHUR_SAXENA)
		{
			const double sum1 = Dot(this->nspecies_, moleFractions, this->etaSpecies_.data());
			const double sum2 = UDot(this->nspecies_, moleFractions, this->etaSpecies_.data());
			const double etamix = 0.50 * (sum1 + 1. / sum2);
			return etamix;
		}
	
		return 0.;
	}

	void TransportPropertiesMap_CHEMKIN::gammaMix(double* gammamix, const double* moleFractions)
	{	
		// Reset
		for(unsigned int k=0;k<this->nspecies_;k++)
			sum_diffusion_coefficients[k] = 0.;

		// Adjust mole fractions
		for(unsigned int i=0;i<this->nspecies_;i++)
			x_corrected[i] = (moleFractions[i]+threshold_)/sum_threshold_;

		double MWmix = 0.;
		for(unsigned int i=0;i<this->nspecies_;i++)
			MWmix += x_corrected[i]*M[i];
		
		// a. Evaluating Mass Diffusion coefficients (mixture averaged)
		const double *d = this->gammaSpecies_.data();
		for(unsigned int k=0;k<this->nspecies_;k++)
			for(unsigned int j=k+1;j<this->nspecies_;j++)
			{
				sum_diffusion_coefficients[j] += x_corrected[k] * (*d);
				sum_diffusion_coefficients[k] += x_corrected[j] * (*d);
				d++;
			}
			
		// b. Evaluating Mass Diffusion coefficients (mixture averaged)
		for(unsigned int k=0;k<this->nspecies_;k++)
			gammamix[k] = (MWmix - x_corrected[k]*M[k]) / (MWmix*sum_diffusion_coefficients[k]);
	}

	void TransportPropertiesMap_CHEMKIN::bundling_gammaMix(double* gammamix, const double* moleFractions)
	{
		// Reset
		for (unsigned int k = 0; k<this->bundling_number_groups_; k++)
			this->bundling_sum_diffusion_coefficients_[k] = 0.;

		// Adjust mole fractions
		for (unsigned int i = 0; i<this->nspecies_; i++)
			x_corrected[i] = (moleFractions[i] + threshold_) / sum_threshold_;

		double MWmix = 0.;
		for (unsigned int i = 0; i<this->nspecies_; i++)
			MWmix += x_corrected[i] * M[i];

		for (unsigned int k = 0; k < this->bundling_number_groups_; k++)
		{
			this->bundling_sum_x_groups_[k] = 0.;
			for (unsigned int j = 0; j < this->bundling_groups_[k].size(); j++)
				this->bundling_sum_x_groups_[k] += x_corrected[this->bundling_groups_[k][j]];
		}
			
		// a. Evaluating Mass Diffusion coefficients (mixture averaged)
		const double *d = this->bundling_gammaSpecies_;
		for (unsigned int k = 0; k<this->bundling_number_groups_; k++)
		for (unsigned int j = k + 1; j<this->bundling_number_groups_; j++)
		{
			this->bundling_sum_diffusion_coefficients_[j] += this->bundling_sum_x_groups_[k] * (*d);
			this->bundling_sum_diffusion_coefficients_[k] += this->bundling_sum_x_groups_[j] * (*d);
			d++;
		}

		// b. Evaluating Mass Diffusion coefficients (mixture averaged)
		for (unsigned int k = 0; k < this->nspecies_; k++)
		{
			const unsigned int group = this->bundling_species_group_[k];
			const double correction = (this->bundling_sum_x_groups_[group] - x_corrected[k])*this->bundling_gammaSpeciesSelfDiffusion_[group];

			gammamix[k] =	(MWmix - x_corrected[k] * M[k]) / 
							(MWmix* (this->bundling_sum_diffusion_coefficients_[group] + correction));
		}
	}

	void TransportPropertiesMap_CHEMKIN::tetaMix(double* tetamix, const double* moleFractions)
	{
		for (unsigned int j=0;j<this->nspecies_;j++)
			tetamix[j] = 0.;

		for (unsigned int k=0;k<iThermalDiffusionRatios_.size();k++)
		{
			unsigned int iSpecies = iThermalDiffusionRatios_[k]-1;
			unsigned int i = k*this->nspecies_;
			for (unsigned int j=0;j<this->nspecies_;j++)
			{
				tetamix[iSpecies] += this->tetaSpecies_[i++]*moleFractions[iSpecies]*moleFractions[j];
			}
		}
	}

	double TransportPropertiesMap_CHEMKIN::kPlanckMix(const double* moleFractions)
	{
		const double uT = 1000. / this->T_;
		const double T = this->T_;

		// 1. Water [1/m/bar]
		const double K_H2O = -0.23093 + uT*(-1.1239 + uT*(9.4153 + uT*(-2.9988 + uT*(0.51382 + uT*(-1.8684e-5)))));

		// 2. Carbon Dioxide [1/m/bar]
		// From: Huaqiang CHU, Mingyan GU, Huaichun ZHOU, Fengshan LIU
		//       Calculations of narrow-band transimissities and the Planck mean absorption coefficients 
		//       of real gases using line - by - line and statistical narrow - band models
		//       Front. Energy 2014, 8(1): 41–48, DOI 10.1007/s11708-013-0292-4
		const double K_CO2 = (-40.129 + T*(3.40989e-1 + T*(-5.39606e-4 + T*(3.60606e-7 + T*(-1.11598e-10 + T*1.31847e-14))))) / 1.01325;

		// 3. Carbon monoxide [1/m/bar]
		double K_CO;
		if (T < 750.)	K_CO = 4.7869 + T*(-0.06953 + T*(2.95775e-4 + T*(-4.25732e-7 + T*2.02894e-10)));
		else            K_CO = 10.09 + T*(-0.01183 + T*(4.7753e-6 + T*(-5.87209e-10 + T*-2.5334e-14)));

		// 4. Methane [1/m/bar]
		const double K_CH4 = 6.6334 + T*(-0.0035686 + T*(1.6682e-08 + T*(2.5611e-10 - 2.6558e-14*T)));

		// Mole fractions
		const double x_H2O = (index_H2O_>0) ? moleFractions[index_H2O_-1] : 0.;
		const double x_CO2 = (index_CO2_>0) ? moleFractions[index_CO2_-1] : 0.;
		const double x_CO  = (index_CO_>0)  ? moleFractions[index_CO_-1]  : 0.;
		const double x_CH4 = (index_CH4_>0) ? moleFractions[index_CH4_-1] : 0.;

		// Absorption coefficients
		const double as_H2O = K_H2O*x_H2O;	// [1/m/bar]
		const double as_CO2 = K_CO2*x_CO2;	// [1/m/bar]
		const double as_CO  = K_CO*x_CO;	// [1/m/bar]
		const double as_CH4 = K_CH4*x_CH4;	// [1/m/bar]

		// Global absorption coefficient
		const double kPlanck = (as_H2O + as_CO2 + as_CO + as_CH4) * (this->P_ / 1.e5);	// [1/m]

		return kPlanck;
	}

	double TransportPropertiesMap_CHEMKIN::kCollision(const unsigned int i, const unsigned int j, const double T)
	{
		const double sigma = 0.50*(sigma_[i] + sigma_[j]);
		const double epsilon_over_kb = std::sqrt(epsilon_over_kb_[i]*epsilon_over_kb_[j]);
		const double mu = mu_[i]*mu_[j]/(mu_[i] + mu_[j]);
		const double TStar = T / epsilon_over_kb;

		return	PhysicalConstants::alphaCollision *
				(sigma*sigma)*Omega11(TStar)*std::sqrt(T/mu);
	}

	double TransportPropertiesMap_CHEMKIN::Omega11(const double TStar)
	{
		return	1.16145*std::pow(TStar, -0.14874) +
				0.52487*std::exp(-0.7732*TStar) +
				2.16178*std::exp(-2.437887*TStar);
	}

	double TransportPropertiesMap_CHEMKIN::Omega11(const double TStar, const double DStar)
	{
		return OpenSMOKE::CollisionIntegral11(TStar, DStar);
	}

	void TransportPropertiesMap_CHEMKIN::Test(const int nLoops, const double& T, int* index)
	{
		double lambdamix, etamix;
		Eigen::VectorXd gammamix(this->nspecies_);
		Eigen::VectorXd tetamix(this->nspecies_);
		
		// Composition (mole fractions)
		Eigen::VectorXd x(this->nspecies_);
		for(unsigned int i=0;i<this->nspecies_;i++)
			x(i) = 1./double(this->nspecies_);

		// Loops
		unsigned int speciesLoopsLambda	= nLoops*125;
		unsigned int speciesLoopsEta		= nLoops*125;
		unsigned int speciesLoopsGamma	= nLoops*1;
		unsigned int speciesLoopsTeta    = nLoops*25;
		unsigned int mixtureLoopsLambda	= nLoops*100;
		unsigned int mixtureLoopsEta		= nLoops*1;
		unsigned int mixtureLoopsGamma	= nLoops*1;
		unsigned int mixtureLoopsTeta	= nLoops*5;


		// Times
		double speciesLambdaTime, speciesEtaTime, speciesGammaTime, speciesTetaTime;
		double mixtureLambdaTime, mixtureEtaTime, mixtureGammaTime, mixtureTetaTime;
                
                SetTemperature(T);
                SetPressure(1e5);

		{
			std::cout << "Species Lambda..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsLambda;k++)
			{
				lambda();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesLambdaTime = tEnd - tStart;
		}

		{
			std::cout << "Species Eta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsEta;k++)
			{
				eta();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesEtaTime = tEnd - tStart;
		}

		{
			std::cout << "Species Gamma..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsGamma;k++)
			{
				gamma();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesGammaTime = tEnd - tStart;
		}

		{
			std::cout << "Species Teta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=speciesLoopsTeta;k++)
			{
				teta();
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesTetaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Lambda..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsLambda;k++)
			{
				lambdamix = lambdaMix(x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureLambdaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Eta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsEta;k++)
			{
				etamix = etaMix(x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureEtaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Gamma..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsGamma;k++)
			{
				gammaMix(gammamix.data(), x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureGammaTime = tEnd - tStart;
		}

		{
			std::cout << "Mixture Teta..." << std::endl;
			double tStart = OpenSMOKEGetCpuTime();
			for(unsigned int k=1;k<=mixtureLoopsTeta;k++)
			{
				tetaMix(tetamix.data(), x.data());
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureTetaTime = tEnd - tStart;
		}

		double lambdaTime	= mixtureLambdaTime/double(mixtureLoopsLambda)	+ speciesLambdaTime/double(speciesLoopsLambda);
		double etaTime		= mixtureEtaTime/double(mixtureLoopsEta)		+ speciesEtaTime/double(speciesLoopsEta);
		double gammaTime	= mixtureGammaTime/double(mixtureLoopsGamma)	+ speciesGammaTime/double(speciesLoopsGamma);
		double tetaTime		= mixtureTetaTime/double(mixtureLoopsTeta)		+ speciesTetaTime/double(speciesLoopsTeta);
		double totTime		= lambdaTime + etaTime + gammaTime + tetaTime;

		std::ofstream fBenchmark;
		fBenchmark.open("Benchmark.plus", std::ios::out);
		fBenchmark.setf(std::ios::scientific);

		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;
		fBenchmark << "                                       PROPERTY VALUES                                                   " << std::endl;
		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;

		fBenchmark << "Thermal conductivities     ";
		fBenchmark << this->lambdaSpecies_(index[1]-1)	<< " " << this->lambdaSpecies_(index[2]-1) << " " << this->lambdaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Dynamic viscosities        ";
		fBenchmark << this->etaSpecies_(index[1]-1)		<< " " << this->etaSpecies_(index[2]-1) << " " << this->etaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Mass diffusivities         ";
		fBenchmark << this->gammaSpecies_(index[1]-1)	<< " " << this->gammaSpecies_(index[2]-1) << " " << this->gammaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Therm. diff ratios         ";
		fBenchmark << this->tetaSpecies_(index[1]-1)	<< " " << this->tetaSpecies_(index[2]-1) << " " << this->tetaSpecies_(index[3]-1) << std::endl; 

		fBenchmark << "Thermal conductivity (mix) ";
		fBenchmark << lambdamix << std::endl;

		fBenchmark << "Dynamic viscosity (mix)    ";
		fBenchmark << etamix << std::endl;

		fBenchmark << "Mass diffusivities (mix)   ";
		fBenchmark << gammamix(index[1]-1) << " " << gammamix(index[2]-1) << " " << gammamix(index[3]-1) << std::endl;

		fBenchmark << "Th. diff. ratios (mix)     ";
		fBenchmark << tetamix(index[1]-1) << " " << tetamix(index[2]-1) << " " << tetamix(index[3]-1) << std::endl;
	
		fBenchmark << std::endl;


		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;
		fBenchmark << "                                      CPU TIME DETAILS                                                   " << std::endl;
		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;

		fBenchmark << "Thermal conductivity (T)       ";
		fBenchmark << speciesLoopsLambda	<< " " << speciesLambdaTime << " " << speciesLambdaTime/double(speciesLoopsLambda)*1000.	<< std::endl;
		
		fBenchmark << "Dynamic Viscosity (T)          ";
		fBenchmark << speciesLoopsEta		<< " " << speciesEtaTime	<< " " << speciesEtaTime/double(speciesLoopsEta)*1000.			<< std::endl;
		
		fBenchmark << "Mass Diffusivities (T)         ";		
		fBenchmark << speciesLoopsGamma		<< " " << speciesGammaTime	<< " " << speciesGammaTime/double(speciesLoopsGamma)*1000.		<< std::endl;
		
		fBenchmark << "Thermal diffusion ratios (T)   ";					
		fBenchmark << speciesLoopsTeta		<< " " << speciesTetaTime	<< " " << speciesTetaTime/double(speciesLoopsTeta)*1000.		<< std::endl;
		
		fBenchmark << "Thermal conductivity (mix)     ";	
		fBenchmark << mixtureLoopsLambda	<< " " << mixtureLambdaTime << " " << mixtureLambdaTime/double(mixtureLoopsLambda)*1000.	<< std::endl;
		
		fBenchmark << "Dynamic Viscosity (mix)        ";
		fBenchmark << mixtureLoopsEta		<< " " << mixtureEtaTime	<< " " << mixtureEtaTime/double(mixtureLoopsEta)*1000.			<< std::endl;
		
		fBenchmark << "Mass Diffusivities (mix)       ";
		fBenchmark << mixtureLoopsGamma		<< " " << mixtureGammaTime	<< " " << mixtureGammaTime/double(mixtureLoopsGamma)*1000.		<< std::endl;
		
		fBenchmark << "Thermal diffusion ratios (mix) ";	
		fBenchmark << mixtureLoopsTeta		<< " " << mixtureTetaTime	<< " " << mixtureTetaTime/double(mixtureLoopsTeta)*1000.		<< std::endl;
		
		fBenchmark << std::endl;


		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;
		fBenchmark << "                                      CPU TIME SUMMARY                                                   " << std::endl;
		fBenchmark << "---------------------------------------------------------------------------------------------------------" << std::endl;

		fBenchmark << "Thermal conductivity (tot)      ";		
		fBenchmark << lambdaTime*1000.		<< " ms - " << lambdaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "Dynamic Viscosity (tot)         ";
		fBenchmark << etaTime*1000.			<< " ms - " << etaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "Mass Diffusivities (tot)        ";
		fBenchmark << gammaTime*1000.		<< " ms - " << gammaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "Thermal diffusion ratios (tot)  ";
		fBenchmark << tetaTime*1000.		<< " ms - " << tetaTime/totTime*100 << " %" << std::endl;

		fBenchmark << "All the properties  (tot)       ";
		fBenchmark << totTime*1000.			<< " ms - " << totTime/totTime*100 << " %" << std::endl;

		fBenchmark.close();
	}
}

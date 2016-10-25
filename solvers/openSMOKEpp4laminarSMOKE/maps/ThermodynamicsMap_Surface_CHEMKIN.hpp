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

#include "math/OpenSMOKEUtilities.h"

namespace OpenSMOKE
{
	template<typename map> 
	ThermodynamicsMap_Surface_CHEMKIN<map>::ThermodynamicsMap_Surface_CHEMKIN(const unsigned int nSpecies, const unsigned int nPoints)
	{
		this->nspecies_ = nSpecies;
		this->npoints_ = nPoints;
		MemoryAllocation();
	}

	template<typename map> 
	ThermodynamicsMap_Surface_CHEMKIN<map>::ThermodynamicsMap_Surface_CHEMKIN(rapidxml::xml_document<>& doc, const unsigned int nPoints)
	{
		this->npoints_ = nPoints;
	
		ImportSpeciesFromXMLFile(doc);
		this->ImportElementsFromXMLFile(doc);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(doc);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::MemoryAllocation()
	{
		Cp_LT = new double[5*this->nspecies_];	
		Cp_HT = new double[5*this->nspecies_];
		DH_LT = new double[6*this->nspecies_];
		DH_HT = new double[6*this->nspecies_];
		DS_LT = new double[6*this->nspecies_];
		DS_HT = new double[6*this->nspecies_];
		TL = new double[this->nspecies_];	
		TH = new double[this->nspecies_];
		TM = new double[this->nspecies_];

		ChangeDimensions(this->nspecies_, &this->MW_, true);
		ChangeDimensions(this->nspecies_, &this->uMW_, true);

		ChangeDimensions(this->nspecies_*this->npoints_, &species_cp_over_R_, true);
		ChangeDimensions(this->nspecies_*this->npoints_, &species_h_over_RT_, true);
		ChangeDimensions(this->nspecies_*this->npoints_, &species_g_over_RT_, true);
		ChangeDimensions(this->nspecies_*this->npoints_, &species_s_over_R_, true);

		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::SetTemperature(const map& T)
	{
		this->T_ = T;
		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::SetPressure(const map& P)
	{
		this->P_ = P;
	}

	
	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::SetCoefficients(const unsigned k, const double* coefficients)
	{
		const double one_third = 1./3.;

		// Specific heat: high temperature
		{
			unsigned int i = k*5;
			Cp_HT[i++] = coefficients[0];
			Cp_HT[i++] = coefficients[1];
			Cp_HT[i++] = coefficients[2];
			Cp_HT[i++] = coefficients[3];
			Cp_HT[i++] = coefficients[4];
		}

		// Specific heat: low temperature
		{
			unsigned int i = k*5;
			Cp_LT[i++] = coefficients[7];
			Cp_LT[i++] = coefficients[8];
			Cp_LT[i++] = coefficients[9];
			Cp_LT[i++] = coefficients[10];
			Cp_LT[i++] = coefficients[11];
		}

		// Enthalpy: high temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DH_HT[i++] = Cp_HT[j++];
			DH_HT[i++] = 0.50 *Cp_HT[j++];
			DH_HT[i++] = one_third*Cp_HT[j++];
			DH_HT[i++] = 0.25 *Cp_HT[j++];
			DH_HT[i++] = 0.20 *Cp_HT[j++];
			DH_HT[i++] = coefficients[5];
		}
	
		// Enthalpy: low temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DH_LT[i++] = Cp_LT[j++];
			DH_LT[i++] = 0.50 *Cp_LT[j++];
			DH_LT[i++] = one_third*Cp_LT[j++];
			DH_LT[i++] = 0.25 *Cp_LT[j++];
			DH_LT[i++] = 0.20 *Cp_LT[j++];
			DH_LT[i++] = coefficients[12];
		}

		// Entropy: high temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DS_HT[i++] = Cp_HT[j++];
			DS_HT[i++] = Cp_HT[j++];
			DS_HT[i++] = 0.50 *Cp_HT[j++];
			DS_HT[i++] = one_third*Cp_HT[j++];
			DS_HT[i++] = 0.25 *Cp_HT[j++];
			DS_HT[i++] = coefficients[6];
		}

		// Entropy: low temperature
		{
			unsigned int j = k*5;
			unsigned int i = k*6;
			DS_LT[i++] = Cp_LT[j++];
			DS_LT[i++] = Cp_LT[j++];
			DS_LT[i++] = 0.50 *Cp_LT[j++];
			DS_LT[i++] = one_third*Cp_LT[j++];
			DS_LT[i++] = 0.25 *Cp_LT[j++];
			DS_LT[i++] = coefficients[13];
		}

		// Temperature limits
		{
			TL[k] = coefficients[14];
			TH[k] = coefficients[15];
			TM[k] = coefficients[16];
		}

		this->MW_[k+1] = coefficients[17];
		this->uMW_[k+1] = 1./coefficients[17];
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromASCIIFile(std::ifstream& fInput)
	{
		std::cout << " * Reading thermodynamic coefficients of species..." << std::endl;
		double coefficients[18];
		for(unsigned int i=0;i<this->nspecies_;i++)
		{
			for(unsigned int j=0;j<18;j++)
				fInput >> coefficients[j];
			SetCoefficients(i, coefficients);
		}
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		// Names of species
		{
			rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
			rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
			rapidxml::xml_node<>* names_of_species_node = opensmoke_node->first_node("NamesOfSpecies");
			
			rapidxml::xml_node<>* number_of_materials_node = opensmoke_node->first_node("NumberOfMaterials");
			rapidxml::xml_node<>* number_of_site_species_node = opensmoke_node->first_node("NumberOfSiteSpecies");
			rapidxml::xml_node<>* number_of_bulk_species_node = opensmoke_node->first_node("NumberOfBulkSpecies");

			rapidxml::xml_node<>* names_of_site_species_node = opensmoke_node->first_node("SiteSpecies");
			rapidxml::xml_node<>* names_of_bulk_species_node = opensmoke_node->first_node("BulkSpecies");

			try
			{
				this->nspecies_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));
				this->names_.resize(this->nspecies_);
				std::stringstream names_of_species_string;
				names_of_species_string << names_of_species_node->value();
				for(unsigned int i=0;i<this->nspecies_;i++)
					names_of_species_string >> this->names_[i]; 	
				
				// Materials
				number_of_materials_    = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_materials_node->value())));
				number_of_site_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_site_species_node->value())));
				number_of_bulk_species_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_bulk_species_node->value())));
				number_of_gas_species_  = this->nspecies_ - number_of_site_species_ - number_of_bulk_species_;

				// Site species
				{
					vector_names_site_species_.resize(number_of_site_species_);
					std::stringstream names_site_species_string;
					names_site_species_string << names_of_site_species_node->value();
					for(unsigned int i=0;i<number_of_site_species_;i++)
						names_site_species_string >> vector_names_site_species_[i]; 	
				}

				// Bulk species
				{
					vector_names_bulk_species_.resize(number_of_bulk_species_);
					std::stringstream names_bulk_species_string;
					names_bulk_species_string << names_of_bulk_species_node->value();
					for(unsigned int i=0;i<number_of_bulk_species_;i++)
						names_bulk_species_string >> vector_names_bulk_species_[i]; 	
				}

			}
			catch(...)
			{
				ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::ImportSpeciesFromXMLFile", "Error in reading the list of species.");
			}
		}

		{
			matrix_names_site_phases_.resize(number_of_materials_);
			matrix_densities_site_phases_.resize(number_of_materials_);
			matrix_names_site_species_.resize(number_of_materials_);
			matrix_occupancies_site_species_.resize(number_of_materials_);
			matrix_indices_site_species_.resize(number_of_materials_);

			matrix_names_bulk_phases_.resize(number_of_materials_);
			matrix_names_bulk_species_.resize(number_of_materials_);
			matrix_densities_bulk_species_.resize(number_of_materials_);
			matrix_indices_bulk_species_.resize(number_of_materials_);

			rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

			unsigned k = 0;
			for (rapidxml::xml_node<> *current_node = opensmoke_node->first_node("MaterialDescription"); current_node; current_node = current_node->next_sibling("MaterialDescription"))
			{
				std::cout << "Reading material description... " << std::endl;

				rapidxml::xml_attribute<> *pMaterialName = current_node->first_attribute("name");
				rapidxml::xml_attribute<> *pMaterialIndex = current_node->first_attribute("index");
					
				rapidxml::xml_node<>* number_of_sites_node = current_node->first_node("NumberOfSites");
				unsigned int number_of_sites_   = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_sites_node->value())));
				
				matrix_names_site_phases_[k].resize(number_of_sites_);
				matrix_densities_site_phases_[k].resize(number_of_sites_);
				matrix_names_site_species_[k].resize(number_of_sites_);
				matrix_occupancies_site_species_[k].resize(number_of_sites_);
				matrix_indices_site_species_[k].resize(number_of_sites_);

				unsigned j = 0;
				for (rapidxml::xml_node<> *site_node = current_node->first_node("Site"); site_node; site_node = site_node->next_sibling("Site"))
				{
					std::cout << "Reading site... " << std::endl;

					rapidxml::xml_attribute<> *pSiteName = site_node->first_attribute("name");
					rapidxml::xml_attribute<> *pSiteIndex = site_node->first_attribute("index");
					matrix_names_site_phases_[k][j] = pSiteName->value();

					unsigned int n;
					std::stringstream dummy;
					dummy << site_node->value();
					
					dummy >> matrix_densities_site_phases_[k][j];
					dummy >> n;
					matrix_names_site_species_[k][j].resize(n);
					matrix_occupancies_site_species_[k][j].resize(n);
					matrix_indices_site_species_[k][j].resize(n);
					for(unsigned int i=0;i<n;i++)
					{
						dummy >> matrix_names_site_species_[k][j][i]; 	
						dummy >> matrix_occupancies_site_species_[k][j][i]; 	
						dummy >> matrix_indices_site_species_[k][j][i]; 	
					}

					j++;
				}

				rapidxml::xml_node<>* number_of_bulks_node = current_node->first_node("NumberOfBulks");
				unsigned int number_of_bulks_   = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_bulks_node->value())));
				
				matrix_names_bulk_phases_[k].resize(number_of_bulks_);
				matrix_names_bulk_species_[k].resize(number_of_bulks_);
				matrix_densities_bulk_species_[k].resize(number_of_bulks_);
				matrix_indices_bulk_species_[k].resize(number_of_bulks_);
				j=0; 
				for (rapidxml::xml_node<> *bulk_node = current_node->first_node("Bulk"); bulk_node; bulk_node = bulk_node->next_sibling("Bulk"))
				{
					std::cout << "Reading bulk... " << std::endl;
					rapidxml::xml_attribute<> *pBulkName = bulk_node->first_attribute("name");
					rapidxml::xml_attribute<> *pBulkIndex = bulk_node->first_attribute("index");
					matrix_names_bulk_phases_[k][j] = pBulkName->value();

					unsigned int n;
					std::stringstream dummy;
					dummy << bulk_node->value();
					
					dummy >> n;
					matrix_names_bulk_species_[k][j].resize(n);
					matrix_densities_bulk_species_[k][j].resize(n);
					matrix_indices_bulk_species_[k][j].resize(n);
					for(unsigned int i=0;i<n;i++)
					{
						dummy >> matrix_names_bulk_species_[k][j][i]; 	
						dummy >> matrix_densities_bulk_species_[k][j][i]; 	
						dummy >> matrix_indices_bulk_species_[k][j][i]; 	
					}

					j++;
				}

				// Site phases belonging
				{
					vector_site_phases_belonging_.resize(number_of_site_species_);
					for(unsigned int i=0;i<number_of_site_species_;i++)
					{
						bool iFound = false;
						for(unsigned int k=0;k<number_of_site_phases(0);k++)
						{
							if (iFound == false)
							{
								for(unsigned int kk=0;kk<matrix_names_site_species_[0][k].size();kk++)
								{
									if (vector_names_site_species_[i] == matrix_names_site_species_[0][k][kk])
									{
										vector_site_phases_belonging_[i] = k;
										iFound = true;
										break;
									}
								}
							}
						}
					}
				}

				// Vector occupancies site species
				vector_occupancies_site_species_.resize(number_of_site_species_);
				for(unsigned int j=0;j<number_of_site_phases(0);j++)
				{
					for(unsigned int i=0;i<matrix_indices_site_species()[0][j].size();i++)
					{
						unsigned int index = matrix_indices_site_species_[0][j][i];
						vector_occupancies_site_species_[index] = matrix_occupancies_site_species_[0][j][i];
					}
				}

				// Vector densities bulk species
				vector_densities_bulk_species_.resize(number_of_bulk_species_);
				for (unsigned int j = 0; j<matrix_indices_bulk_species_[0].size(); j++)
				{
					for (unsigned int i = 0; i<matrix_indices_bulk_species_[0][j].size(); i++)
					{
						unsigned int index = matrix_indices_bulk_species_[0][j][i];
						vector_densities_bulk_species_[index] = matrix_densities_bulk_species_[0][j][i];
					}
				}

				/*
				for(unsigned int j=0;j<matrix_name_site_species_[k].size();j++)
					for(unsigned int i=0;i<matrix_name_site_species_[k][j].size();i++)
						std::cout << matrix_name_site_species_[k][j][i] << " " << matrix_occupancy_site_species_[k][j][i] << " " << matrix_index_site_species_[k][j][i] << std::endl;


				for(unsigned int j=0;j<matrix_name_bulk_species_[k].size();j++)
					for(unsigned int i=0;i<matrix_name_bulk_species_[k][j].size();i++)
						std::cout << matrix_name_bulk_species_[k][j][i] << " " << matrix_density_bulk_species_[k][j][i] << " " << matrix_index_bulk_species_[k][j][i] << std::endl;
						*/

				k++;
			}
		}
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		std::cout << " * Reading thermodynamic coefficients of species from XML file..." << std::endl;

		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

		rapidxml::xml_node<>* thermodynamics_node = opensmoke_node->first_node("Thermodynamics");
		if (thermodynamics_node != 0)
		{
			std::string thermodynamics_type = thermodynamics_node->first_attribute("type")->value();			
			if (thermodynamics_type == "NASA")
			{
				rapidxml::xml_node<>* nasa_coefficients_node = thermodynamics_node->first_node("NASA-coefficients");
				
				std::stringstream nasa_coefficients;
				nasa_coefficients << nasa_coefficients_node->value();
			
				double coefficients[18];
				for(unsigned int i=0;i<this->nspecies_;i++)
				{
					for(unsigned int j=0;j<18;j++)
						nasa_coefficients >> coefficients[j];
					SetCoefficients(i, coefficients);
				}
			}
			else
				ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile", "Thermodynamics type not supported!"); 
		}
		else
			ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::ImportCoefficientsFromXMLFile", "Thermodynamics tag was not found!");
	}

	template<typename map> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::MolecularWeight_From_MoleFractions(map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::MolecularWeight_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::MolecularWeight_From_MassFractions(map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::MolecularWeight_From_MassFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::MassFractions_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& y, map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::MassFractions_From_MoleFractions", "Function not available for surface thermodynamics!");
	}
	
	template<typename map> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::MoleFractions_From_MassFractions(OpenSMOKE::OpenSMOKEVectorDouble& x, map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::MoleFractions_From_MassFractions", "Function not available for surface thermodynamics!");
	}

	template<> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<double>::cp_over_R()
	{
		if (cp_must_be_recalculated_ == true)
		{
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;

			unsigned int j = 0;
			for (unsigned int k=1;k<=this->nspecies_;k++)
			{
				species_cp_over_R_[k] = (this->T_>TM[k-1]) ?		Cp_HT[j] + this->T_*Cp_HT[j+1] + T2*Cp_HT[j+2] + T3*Cp_HT[j+3] + T4*Cp_HT[j+4] :
															Cp_LT[j] + this->T_*Cp_LT[j+1] + T2*Cp_LT[j+2] + T3*Cp_LT[j+3] + T4*Cp_LT[j+4] ;
				j += 5;
			}

			cp_must_be_recalculated_ = false;
		}
	}

	template<typename map> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::cp_over_R()
	{
		if (cp_must_be_recalculated_ == true)
		{
			int unsigned count = 1;
			for (unsigned int i=1;i<=this->npoints_;i++)
			{
				const double T1 = this->T_[i];
				const double T2 = T1*T1;
				const double T3 = T2*T1;
				const double T4 = T3*T1;

				unsigned int j = 0;
				for (unsigned int k=1;k<=this->nspecies_;k++)
				{
					species_cp_over_R_[count++] = (T1>TM[k-1]) ?	Cp_HT[j] + T1*Cp_HT[j+1] + T2*Cp_HT[j+2] + T3*Cp_HT[j+3] + T4*Cp_HT[j+4] :
																	Cp_LT[j] + T1*Cp_LT[j+1] + T2*Cp_LT[j+2] + T3*Cp_LT[j+3] + T4*Cp_LT[j+4] ;
					j += 5;
				}
			}
			
			cp_must_be_recalculated_ = false;
		}
	}

	template<> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<double>::h_over_RT()
	{
		if (h_must_be_recalculated_ == true)
		{
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;
			const double uT = 1./this->T_;

			int j = 0;
			for (unsigned int k=1;k<=this->nspecies_;k++)
			{
				species_h_over_RT_[k] = (this->T_>TM[k-1]) ?	DH_HT[j] + this->T_*DH_HT[j+1] + T2*DH_HT[j+2] + T3*DH_HT[j+3] + T4*DH_HT[j+4] + uT*DH_HT[j+5]:
														DH_LT[j] + this->T_*DH_LT[j+1] + T2*DH_LT[j+2] + T3*DH_LT[j+3] + T4*DH_LT[j+4] + uT*DH_LT[j+5];
				j += 6;
			}

			h_must_be_recalculated_ = false;
		}
	}
	
	template<typename map>
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::h_over_RT()
	{
		if (h_must_be_recalculated_ == true)
		{
			unsigned int count = 1;
			for (unsigned int i=1;i<=this->npoints_;i++)
			{
				const double T1 = this->T_[i];
				const double T2 = T1*T1;
				const double T3 = T2*T1;
				const double T4 = T3*T1;
				const double uT1 = 1./T1;

				unsigned int j = 0;
				for (unsigned int k=1;k<=this->nspecies_;k++)
				{
					species_h_over_RT_[count++] = (T1>TM[k-1]) ?	DH_HT[j] + T1*DH_HT[j+1] + T2*DH_HT[j+2] + T3*DH_HT[j+3] + T4*DH_HT[j+4] + uT1*DH_HT[j+5]:
																	DH_LT[j] + T1*DH_LT[j+1] + T2*DH_LT[j+2] + T3*DH_LT[j+3] + T4*DH_LT[j+4] + uT1*DH_LT[j+5];
					j += 6;
				}
			}

			h_must_be_recalculated_ = false;
		}
	}
	
	template<> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<double>::s_over_R()
	{
		if (s_must_be_recalculated_ == true)
		{
			const double logT = log(this->T_);
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;

			unsigned int j = 0;
			for (unsigned int k=1;k<=this->nspecies_;k++)
			{
				species_s_over_R_[k] = (this->T_>TM[k-1]) ?	DS_HT[j]*logT + this->T_*DS_HT[j+1] + T2*DS_HT[j+2] + T3*DS_HT[j+3] + T4*DS_HT[j+4] + DS_HT[j+5] :
															DS_LT[j]*logT + this->T_*DS_LT[j+1] + T2*DS_LT[j+2] + T3*DS_LT[j+3] + T4*DS_LT[j+4] + DS_LT[j+5] ;
				j += 6;
			}

			s_must_be_recalculated_ = false;
		}
	}

	template<> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<double>::g_over_RT()
	{
		h_over_RT();
		s_over_R();

		Sub(species_h_over_RT_, species_s_over_R_, &species_g_over_RT_);
	}
	
	template<typename map> 
	inline void ThermodynamicsMap_Surface_CHEMKIN<map>::s_over_R()
	{
		if (s_must_be_recalculated_ == true)
		{
			unsigned int count = 1;
			for (unsigned int i=1;i<=this->npoints_;i++)
			{
				const double T1 = this->T_[i];
				const double T2 = T1*T1;
				const double T3 = T2*T1;
				const double T4 = T3*T1;
				const double logT1 = log(T1);

				unsigned int j = 0;
				for (unsigned int k=1;k<=this->nspecies_;k++)
				{
					species_s_over_R_[count++] = (T1>TM[k-1]) ?	DS_HT[j]*logT1 + T1*DS_HT[j+1] + T2*DS_HT[j+2] + T3*DS_HT[j+3] + T4*DS_HT[j+4] + DS_HT[j+5] :
																DS_LT[j]*logT1 + T1*DS_LT[j+1] + T2*DS_LT[j+2] + T3*DS_LT[j+3] + T4*DS_LT[j+4] + DS_LT[j+5] ;
					j += 6;
				}
			}
		
			s_must_be_recalculated_ = false;
		}
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::cpMolar_Mixture_From_MoleFractions(map& cpmix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::cpMolar_Mixture_From_MoleFractions", "Function not available for surface thermodynamics!");
	}
	
	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::hMolar_Mixture_From_MoleFractions(map& hmix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::hMolar_Mixture_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::sMolar_Mixture_From_MoleFractions(map& smix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::sMolar_Mixture_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::uMolar_Mixture_From_MoleFractions(map& umix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::uMolar_Mixture_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::gMolar_Mixture_From_MoleFractions(map& gmix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::gMolar_Mixture_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::aMolar_Mixture_From_MoleFractions(map& amix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::aMolar_Mixture_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::cpMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& cp_species)
	{
		cp_over_R();
		Product(PhysicalConstants::R_J_kmol, species_cp_over_R_, &cp_species);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::hMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& h_species)
	{
		h_over_RT();
		Product(PhysicalConstants::R_J_kmol*this->T_, species_h_over_RT_, &h_species);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::sMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& s_species)
	{
		s_over_R();
		Product(PhysicalConstants::R_J_kmol, species_s_over_R_, &s_species);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::sMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& s_species, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::sMolar_Species_MixtureAveraged_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::uMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& u_species)
	{
		h_over_RT();
		Add(species_h_over_RT_, -1., &u_species);
		Product(PhysicalConstants::R_J_kmol*this->T_, &u_species);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::gMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& g_species)
	{
		h_over_RT();
		s_over_R();
		Sub(species_h_over_RT_, species_s_over_R_, &g_species);
		Product(PhysicalConstants::R_J_kmol*this->T_, &g_species);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::gMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& g_species, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::gMolar_Species_MixtureAveraged_From_MoleFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::aMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& a_species)
	{
		h_over_RT();
		s_over_R();
		Sub(species_h_over_RT_, species_s_over_R_, &a_species);
		a_species -= 1.;
		Product(PhysicalConstants::R_J_kmol*this->T_, &a_species);
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::aMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& a_species, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		gMolar_Species_MixtureAveraged_From_MoleFractions(a_species, x);
		a_species -= PhysicalConstants::R_J_kmol*this->T_;
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::DerivativesOfConcentrationsWithRespectToMassFractions(const double cTot, const double MW, const OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEMatrixDouble* dc_over_omega)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::DerivativesOfConcentrationsWithRespectToMassFractions", "Function not available for surface thermodynamics!");
	}

	template<typename map> 
	void ThermodynamicsMap_Surface_CHEMKIN<map>::Test(const int nLoops, const map& T, int* index)
	{
		ErrorMessage("ThermodynamicsMap_Surface_CHEMKIN<map>::Test", "Function not available for surface thermodynamics!");
	}

}

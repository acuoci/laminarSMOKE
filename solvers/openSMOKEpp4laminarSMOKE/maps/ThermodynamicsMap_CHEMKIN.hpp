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
	ThermodynamicsMap_CHEMKIN::ThermodynamicsMap_CHEMKIN(const unsigned int nSpecies)
	{
		this->nspecies_ = nSpecies;
        this->verbose_output_ = true;
                
		MemoryAllocation();
	}

	ThermodynamicsMap_CHEMKIN::ThermodynamicsMap_CHEMKIN(rapidxml::xml_document<>& doc)
	{
		this->verbose_output_ = true;
	
		ImportSpeciesFromXMLFile(doc);
		this->ImportElementsFromXMLFile(doc);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(doc);
	}
        
	ThermodynamicsMap_CHEMKIN::ThermodynamicsMap_CHEMKIN(rapidxml::xml_document<>& doc, bool verbose)
	{
		this->verbose_output_ = verbose;
	
		ImportSpeciesFromXMLFile(doc);
		this->ImportElementsFromXMLFile(doc);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(doc);
	}
    
	ThermodynamicsMap_CHEMKIN::ThermodynamicsMap_CHEMKIN( const ThermodynamicsMap_CHEMKIN& rhs )
    {
        CopyFromMap(rhs);
    }
        
	void ThermodynamicsMap_CHEMKIN::CopyFromMap( const ThermodynamicsMap_CHEMKIN& rhs )
    {
        this->nspecies_ = rhs.nspecies_;
            
        MemoryAllocation();
             
		this->names_ =  rhs.names_;
        this->MW__   =  rhs.MW__;

        this->atomic_composition_ = rhs.atomic_composition_;
		this->elements_ = rhs.elements_; 
		this->verbose_output_ = rhs.verbose_output_;
            
        for (unsigned int i=0;i<5*this->nspecies_;i++)
            this->Cp_LT[i] = rhs.Cp_LT[i];

        for (unsigned int i=0;i<5*this->nspecies_;i++)
            this->Cp_HT[i] = rhs.Cp_HT[i];
            
        for (unsigned int i=0;i<6*this->nspecies_;i++)
            this->DH_LT[i] = rhs.DH_LT[i];

        for (unsigned int i=0;i<6*this->nspecies_;i++)
            this->DH_HT[i] = rhs.DH_HT[i];
            
        for (unsigned int i=0;i<6*this->nspecies_;i++)
            this->DS_LT[i] = rhs.DS_LT[i];

        for (unsigned int i=0;i<6*this->nspecies_;i++)
            this->DS_HT[i] = rhs.DS_HT[i];      
            
        for (unsigned int i=0;i<this->nspecies_;i++)
            this->TL[i] = rhs.TL[i];      
            
        for (unsigned int i=0;i<this->nspecies_;i++)
            this->TH[i] = rhs.TH[i];      

        for (unsigned int i=0;i<this->nspecies_;i++)
            this->TM[i] = rhs.TM[i]; 
    }
	  
	ThermodynamicsMap_CHEMKIN::~ThermodynamicsMap_CHEMKIN(void)
	{
		this->names_.clear();		
		this->elements_.clear();		
		delete[] Cp_LT;	
		delete[] Cp_HT;
		delete[] DH_LT;
		delete[] DH_HT;
		delete[] DS_LT;
		delete[] DS_HT;
		delete[] TL;	
		delete[] TH;
		delete[] TM;
	}

	void ThermodynamicsMap_CHEMKIN::MemoryAllocation()
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

		this->MW__.resize(this->nspecies_);

		species_cp_over_R__.resize(this->nspecies_);
		std::fill(species_cp_over_R__.begin(), species_cp_over_R__.end(), 0.);
		species_h_over_RT__.resize(this->nspecies_);
		std::fill(species_h_over_RT__.begin(), species_h_over_RT__.end(), 0.);
		species_g_over_RT__.resize(this->nspecies_);
		std::fill(species_g_over_RT__.begin(), species_g_over_RT__.end(), 0.);
		species_s_over_R__.resize(this->nspecies_);
		std::fill(species_s_over_R__.begin(), species_s_over_R__.end(), 0.);

		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}

	void ThermodynamicsMap_CHEMKIN::SetTemperature(const double& T)
	{
		this->T_ = T;
		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}
 
	void ThermodynamicsMap_CHEMKIN::SetPressure(const double& P)
	{
		this->P_ = P;
	}

	void ThermodynamicsMap_CHEMKIN::Change_a_HT(const unsigned int species, const unsigned int j, const double value)
	{

		const unsigned int i1 = species * 5 + (j - 1);
		const unsigned int i2 = species * 6 + (j - 1);
			
		if (j == 1)
		{
			Cp_HT[i1] = value;
			DH_HT[i2] = value;
			DS_HT[i2] = value;
		}
		else if (j == 2)
		{
			Cp_HT[i1] = value;
			DH_HT[i2] = value/2.;
			DS_HT[i2] = value;
		}
		else if (j == 3)
		{
			Cp_HT[i1] = value;
			DH_HT[i2] = value/3.;
			DS_HT[i2] = value/2.;
		}
		else if (j == 4)
		{
			Cp_HT[i1] = value;
			DH_HT[i2] = value/4.;
			DS_HT[i2] = value/3.;
		}
		else if (j == 5)
		{
			Cp_HT[i1] = value;
			DH_HT[i2] = value/5.;
			DS_HT[i2] = value/4.;
		}
		else if (j == 6)
		{
			DH_HT[i2] = value;
		}
		else if (j == 7)
		{
			DS_HT[i2] = value;
		}
	}

	void ThermodynamicsMap_CHEMKIN::Change_a_LT(const unsigned int species, const unsigned int j, const double value)
	{

		const unsigned int i1 = species * 5 + (j - 1);
		const unsigned int i2 = species * 6 + (j - 1);

		if (j == 1)
		{
			Cp_LT[i1] = value;
			DH_LT[i2] = value;
			DS_LT[i2] = value;
		}
		else if (j == 2)
		{
			Cp_LT[i1] = value;
			DH_LT[i2] = value / 2.;
			DS_LT[i2] = value;
		}
		else if (j == 3)
		{
			Cp_LT[i1] = value;
			DH_LT[i2] = value / 3.;
			DS_LT[i2] = value / 2.;
		}
		else if (j == 4)
		{
			Cp_LT[i1] = value;
			DH_LT[i2] = value / 4.;
			DS_LT[i2] = value / 3.;
		}
		else if (j == 5)
		{
			Cp_LT[i1] = value;
			DH_LT[i2] = value / 5.;
			DS_LT[i2] = value / 4.;
		}
		else if (j == 6)
		{
			DH_LT[i2] = value;
		}
		else if (j == 7)
		{
			DS_LT[i2] = value;
		}
	}

	void ThermodynamicsMap_CHEMKIN::SetCoefficients(const unsigned k, const double* coefficients)
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

		this->MW__[k] = coefficients[17];
	}
 
	void ThermodynamicsMap_CHEMKIN::ImportCoefficientsFromASCIIFile(std::ifstream& fInput)
	{
                if(verbose_output_ == true)
		    std::cout << " * Reading thermodynamic coefficients of species..." << std::endl;
          
		double coefficients[18];
		for(unsigned int i=0;i<this->nspecies_;i++)
		{
			for(unsigned int j=0;j<18;j++)
				fInput >> coefficients[j];
			SetCoefficients(i, coefficients);
		}
	}
 
	void ThermodynamicsMap_CHEMKIN::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
		rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("NumberOfSpecies");
		rapidxml::xml_node<>* names_of_species_node = opensmoke_node->first_node("NamesOfSpecies");
		try
		{
			this->nspecies_ = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(number_of_species_node->value())));
			this->names_.resize(this->nspecies_);
			std::stringstream names_of_species_string;
			names_of_species_string << names_of_species_node->value();
			for(unsigned int i=0;i<this->nspecies_;i++)
				names_of_species_string >> this->names_[i]; 							
		}
		catch(...)
		{
			ErrorMessage("ThermodynamicsMap_CHEMKIN::ImportSpeciesFromXMLFile", "Error in reading the list of species.");
		}
	}
 
	void ThermodynamicsMap_CHEMKIN::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
	{
                if(verbose_output_ == true)
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
				ErrorMessage("ThermodynamicsMap_CHEMKIN::ImportCoefficientsFromXMLFile", "Thermodynamics type not supported!"); 
		}
		else
			ErrorMessage("ThermodynamicsMap_CHEMKIN::ImportCoefficientsFromXMLFile", "Thermodynamics tag was not found!");
	}
 
	inline double ThermodynamicsMap_CHEMKIN::MolecularWeight_From_MoleFractions(const double* x)
	{
		return Dot(this->nspecies_, x, this->MW__.data());
	}

	inline double ThermodynamicsMap_CHEMKIN::MolecularWeight_From_MassFractions(const double* y)
	{
		return 1./UDot(this->nspecies_, y, this->MW__.data());
	}

	inline void ThermodynamicsMap_CHEMKIN::MassFractions_From_MoleFractions(double* y, double& MW, const double* x)
	{
		MW = Dot(this->nspecies_, x ,this->MW__.data());
		ElementByElementProduct(this->nspecies_, x, this->MW__.data(), y);
		Prod(this->nspecies_, 1./MW, y);
	}
	
	inline void ThermodynamicsMap_CHEMKIN::MoleFractions_From_MassFractions(double* x, double& MW, const double* y)
	{
		MW = 1./UDot(this->nspecies_, y, this->MW__.data());
		Prod(this->nspecies_, MW, y, x);
		ElementByElementDivision(this->nspecies_, x, this->MW__.data(), x);
	}

	inline void ThermodynamicsMap_CHEMKIN::cp_over_R()
	{
		if (cp_must_be_recalculated_ == true)
		{
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;

			unsigned int j = 0;
			for (unsigned int k=0;k<this->nspecies_;k++)
			{
				species_cp_over_R__[k] = (this->T_>TM[k]) ?		Cp_HT[j] + this->T_*Cp_HT[j+1] + T2*Cp_HT[j+2] + T3*Cp_HT[j+3] + T4*Cp_HT[j+4] :
																Cp_LT[j] + this->T_*Cp_LT[j+1] + T2*Cp_LT[j+2] + T3*Cp_LT[j+3] + T4*Cp_LT[j+4] ;
				j += 5;
			}

			cp_must_be_recalculated_ = false;
		}
	}

	inline void ThermodynamicsMap_CHEMKIN::h_over_RT()
	{
		if (h_must_be_recalculated_ == true)
		{
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;
			const double uT = 1./this->T_;

			int j = 0;
			for (unsigned int k=0;k<this->nspecies_;k++)
			{
				species_h_over_RT__[k] = (this->T_>TM[k]) ?	DH_HT[j] + this->T_*DH_HT[j+1] + T2*DH_HT[j+2] + T3*DH_HT[j+3] + T4*DH_HT[j+4] + uT*DH_HT[j+5]:
															DH_LT[j] + this->T_*DH_LT[j+1] + T2*DH_LT[j+2] + T3*DH_LT[j+3] + T4*DH_LT[j+4] + uT*DH_LT[j+5];
				j += 6;
			}

			h_must_be_recalculated_ = false;
		}
	}
	
	inline void ThermodynamicsMap_CHEMKIN::s_over_R()
	{
		if (s_must_be_recalculated_ == true)
		{
			const double logT = std::log(this->T_);
			const double T2 = this->T_*this->T_;
			const double T3 = T2*this->T_;
			const double T4 = T3*this->T_;

			unsigned int j = 0;
			for (unsigned int k=0;k<this->nspecies_;k++)
			{
				species_s_over_R__[k] = (this->T_>TM[k]) ?	DS_HT[j]*logT + this->T_*DS_HT[j+1] + T2*DS_HT[j+2] + T3*DS_HT[j+3] + T4*DS_HT[j+4] + DS_HT[j+5] :
															DS_LT[j]*logT + this->T_*DS_LT[j+1] + T2*DS_LT[j+2] + T3*DS_LT[j+3] + T4*DS_LT[j+4] + DS_LT[j+5] ;
				j += 6;
			}

			s_must_be_recalculated_ = false;
		}
	}

	inline void ThermodynamicsMap_CHEMKIN::g_over_RT()
	{
		h_over_RT();
		s_over_R();

		Difference(this->nspecies_, species_h_over_RT__.data(), species_s_over_R__.data(), species_g_over_RT__.data());
	}
	
	double ThermodynamicsMap_CHEMKIN::cpMolar_Mixture_From_MoleFractions(const double* x)
	{
		cp_over_R();
		double cpmix = Dot(this->nspecies_, species_cp_over_R__.data(), x);
		cpmix *= PhysicalConstants::R_J_kmol;
		return cpmix;
	}
	
	double ThermodynamicsMap_CHEMKIN::hMolar_Mixture_From_MoleFractions(const double* x)
	{
		h_over_RT();
		double hmix = Dot(this->nspecies_, species_h_over_RT__.data(), x);
		hmix *= PhysicalConstants::R_J_kmol*this->T_;
		return hmix;
	}

	double ThermodynamicsMap_CHEMKIN::sMolar_Mixture_From_MoleFractions(const double* x)
	{
		s_over_R();

		const double eps = 1.e-32;
		double sum=0.;
		for(unsigned int i=1;i<=this->nspecies_;i++)
			sum += x[i]*std::log(x[i]+eps);
		
		double smix = Dot(this->nspecies_, species_s_over_R__.data(), x) - std::log(this->P_/101325.) - sum;
		smix *= PhysicalConstants::R_J_kmol;
		return smix;
	}
 
	double ThermodynamicsMap_CHEMKIN::uMolar_Mixture_From_MoleFractions(const double* x)
	{
		h_over_RT();
		double umix = Dot(this->nspecies_, species_h_over_RT__.data(), x) - 1.;
		umix *= PhysicalConstants::R_J_kmol*this->T_;
		return umix;
	}

	double ThermodynamicsMap_CHEMKIN::gMolar_Mixture_From_MoleFractions(const double* x)
	{
		const double hmix = hMolar_Mixture_From_MoleFractions(x);
		const double smix = sMolar_Mixture_From_MoleFractions(x);
		const double gmix = hmix-this->T_*smix;
		return gmix;
	}

	double ThermodynamicsMap_CHEMKIN::aMolar_Mixture_From_MoleFractions(const double* x)
	{
		const double umix = uMolar_Mixture_From_MoleFractions(x);
		const double smix = sMolar_Mixture_From_MoleFractions(x);
		const double amix = umix-this->T_*smix;
		return amix;
	}
 
	void ThermodynamicsMap_CHEMKIN::cpMolar_Species(double* cp_species)
	{
		cp_over_R();
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol, species_cp_over_R__.data(), cp_species);
	}

	void ThermodynamicsMap_CHEMKIN::hMolar_Species(double* h_species)
	{
		h_over_RT();
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, species_h_over_RT__.data(), h_species);
	}

	void ThermodynamicsMap_CHEMKIN::sMolar_Species(double* s_species)
	{
		s_over_R();
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol, species_s_over_R__.data(), s_species);
	}

	void ThermodynamicsMap_CHEMKIN::uMolar_Species(double* u_species)
	{
		h_over_RT();
		Sum(this->nspecies_, species_h_over_RT__.data(), -1., u_species);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, u_species);
	}

	void ThermodynamicsMap_CHEMKIN::gMolar_Species(double* g_species)
	{
		h_over_RT();
		s_over_R();
		Difference(this->nspecies_, species_h_over_RT__.data(), species_s_over_R__.data(), g_species);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, g_species);
	}

	void ThermodynamicsMap_CHEMKIN::aMolar_Species(double* a_species)
	{
		h_over_RT();
		s_over_R();
		Difference(this->nspecies_, species_h_over_RT__.data(), species_s_over_R__.data(), a_species);
		Sum(this->nspecies_, -1., a_species);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, a_species);
	}

	void ThermodynamicsMap_CHEMKIN::sMolar_Species_MixtureAveraged_From_MoleFractions(double* s_species, const double* x)
	{
		const double eps = 1.e-32;
		const double log_P_over_Patm = std::log(this->P_/101325.);
		s_over_R();
		for(unsigned int i=1;i<=this->nspecies_;i++)
			s_species[i] = (x[i]>eps) ? species_s_over_R__[i-1] - std::log(x[i]+eps) - log_P_over_Patm : species_s_over_R__[i-1] - log_P_over_Patm;
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol, s_species);
	}

	void ThermodynamicsMap_CHEMKIN::gMolar_Species_MixtureAveraged_From_MoleFractions(double* g_species, const double* x)
	{
		const double eps = 1.e-32;
		const double log_P_pver_Patm = std::log(this->P_/101325.);
		s_over_R();
		h_over_RT();
		for(unsigned int i=1;i<=this->nspecies_;i++)
			g_species[i] = (x[i]>eps) ? species_h_over_RT__[i-1] - (species_s_over_R__[i-1] - std::log(x[i]+eps) - log_P_pver_Patm) : 
			                            species_h_over_RT__[i-1] - (species_s_over_R__[i-1] - log_P_pver_Patm);
		Prod(this->nspecies_, PhysicalConstants::R_J_kmol*this->T_, g_species);
	}


	void ThermodynamicsMap_CHEMKIN::aMolar_Species_MixtureAveraged_From_MoleFractions(double* a_species, const double* x)
	{
		gMolar_Species_MixtureAveraged_From_MoleFractions(a_species, x);
		Sum(this->nspecies_, -PhysicalConstants::R_J_kmol*this->T_, a_species);
	}

	void ThermodynamicsMap_CHEMKIN::DerivativesOfConcentrationsWithRespectToMassFractions(const double cTot, const double MW, const double* omega, Eigen::MatrixXd* dc_over_omega)
	{
		for(unsigned int j=0;j<this->nspecies_;j++)
		{
			const double coefficient = cTot*MW/this->MW__[j];
			for(unsigned int i=0;i<this->nspecies_;i++)
				(*dc_over_omega)(j,i) = -coefficient*MW/this->MW__[i]*omega[j];
			(*dc_over_omega)(j, j) +=coefficient;
		}
	}

	void ThermodynamicsMap_CHEMKIN::Test(const int nLoops, const double& T, int* index)
	{
	}
	
	double ThermodynamicsMap_CHEMKIN::GetTemperatureFromEnthalpyAndMoleFractions(const double H, const double P_Pa, const double* x, const double TFirstGuess)
	{
		const double diff_temperature_to_switch_to_newtons = 15.;
		const unsigned int max_number_of_enlargements = 10;
		
		// Check how good is the first guess temperature
		double absolute_error;
		{
			SetTemperature(TFirstGuess);
			SetPressure(P_Pa);
			const double CpFirstGuess = cpMolar_Mixture_From_MoleFractions(x);
			const double HFirstGuess = hMolar_Mixture_From_MoleFractions(x);
			
			absolute_error = std::fabs(HFirstGuess-H)/CpFirstGuess;
		}

		if (absolute_error < 1e-16)
			return TFirstGuess;

		{
			double TA, TB;
			double HA, HB;
			const bool flag_find_interval = FindTheRightInterval(H, x, TFirstGuess, TA, TB, HA, HB, 100.);
			if (flag_find_interval == false)
				return 0.;

			const double tol  = 1e-14;
			const double xacc = 1e-14;

			const double TC = Brent(H, P_Pa, x, TA, TB, tol);

			return TC;
		}

		return 0;
	}
	
	bool ThermodynamicsMap_CHEMKIN::FindTheRightInterval(	const double H, const double* x, const double TFirstGuess,
															double& TA, double&TB, double& HA, double& HB,
															const double Diff_Temperature)
	{
		const double TMinimum = 250.;
		const double TMaximum = 6000.;
		const unsigned int max_number_of_enlargements = 50;

		TA = TB = TFirstGuess;

		SetTemperature(TA);
		HA = hMolar_Mixture_From_MoleFractions(x);
		
		SetTemperature(TB);
		HB = hMolar_Mixture_From_MoleFractions(x);

		unsigned int count = 0;
		for(;;)
		{
			if ((HA-H)*(HB-H)<=0.)
				return true;

			if ((HA-H)>0.)
			{	
				TA = std::max(TMinimum, TA-Diff_Temperature);
				SetTemperature(TA);
				HA = hMolar_Mixture_From_MoleFractions(x);
			}

			if ((HB-H)<0.)
			{	
				TB = std::min(TMaximum, TB+Diff_Temperature);
				SetTemperature(TB);
				HB = hMolar_Mixture_From_MoleFractions(x);
			}

			count++;
			if (count > max_number_of_enlargements)
				return false;
			//	ErrorMessage("FindTheRightInterval", "Maximum number of enlargements");
		}
		
		return false;
	}

	template<class T>
	inline const T SIGN(const T &a, const T &b)
		{return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

	double ThermodynamicsMap_CHEMKIN::Brent(const double H, const double P_Pa, const double* x, const double x1, const double x2, const double tol)
	{
		const unsigned int 	ITMAX=100;
		const double 		EPS=2.2204460492503131e-16;
	
		double a=x1;
		double b=x2;
		double c=x2;
		double d, e, min1, min2;
		double fa=Function(H, P_Pa, x, a);
		double fb=Function(H, P_Pa, x, b);
		double fc,p,q,r,s,tol1,xm;

		if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0))
			ErrorMessage("Brent", "Root must be bracketed in zbrent");
		fc=fb;
		
		for (unsigned int iter=0;iter<ITMAX;iter++) 
		{
			if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) 
			{
				c=a;
				fc=fa;
				e=d=b-a;
			}

			if (std::fabs(fc) < std::fabs(fb))
			{
				a=b;
				b=c;
				c=a;
				fa=fb;
				fb=fc;
				fc=fa;
			}
		
			tol1=2.0*EPS*std::fabs(b)+0.5*tol;
			xm=0.5*(c-b);
			
			if (std::fabs(xm) <= tol1 || fb == 0.0) 
				return b;
		
			if (std::fabs(e) >= tol1 && std::fabs(fa) > std::fabs(fb)) 
			{
				s=fb/fa;
				if (a == c) 
				{
					p=2.0*xm*s;
					q=1.0-s;
				} 
				else 
				{
					q=fa/fc;
					r=fb/fc;
					p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
					q=(q-1.0)*(r-1.0)*(s-1.0);
				}
			
				if (p > 0.0) 
					q = -q;
			
				p=std::fabs(p);
				min1=3.0*xm*q-std::fabs(tol1*q);
				min2=std::fabs(e*q);
				if (2.0*p < (min1 < min2 ? min1 : min2)) 
				{
					e=d;
					d=p/q;
				} 
				else 
				{
					d=xm;
					e=d;
				}
			} 
			else 
			{
				d=xm;
				e=d;
			}
			a=b;
			fa=fb;
			
			if (std::fabs(d) > tol1)
				b += d;
			else
				b += SIGN(tol1,xm);
			
			fb=Function(H, P_Pa, x, b);
		}
	
		ErrorMessage("Brent", "Maximum number of iterations exceeded in zbrent");
		return 0.;
	}
 
	double ThermodynamicsMap_CHEMKIN::Ridder(const double H, const double P_Pa, const double* x, const double x1, const double x2, const double xacc)
	{
		const int MAXIT=60;
		const double UNUSED=-1.11e30;
		int j;
		double ans,fh,fl,fm,fnew,s,xh,xl,xm,xnew;

		fl=Function(H, P_Pa, x, x1);
		fh=Function(H, P_Pa, x, x2);
		
		if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) 
		{
			xl=x1;
			xh=x2;
			ans=UNUSED;
			for (j=0;j<MAXIT;j++) 
			{
				xm=0.5*(xl+xh);
				fm=Function(H, P_Pa, x, xm);
				s=std::sqrt(fm*fm-fl*fh);
				if (s == 0.0) return ans;
				xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s);
				if (std::fabs(xnew-ans) <= xacc) return ans;
				ans=xnew;
				fnew=Function(H, P_Pa, x, ans);
				if (fnew == 0.0) return ans;
			
				if (SIGN(fm,fnew) != fm)
				{
					xl=xm;
					fl=fm;
					xh=ans;
					fh=fnew;
				} 
				else if (SIGN(fl,fnew) != fl) 
				{
					xh=ans;
					fh=fnew;
				} 
				else if (SIGN(fh,fnew) != fh) 
				{
					xl=ans;
					fl=fnew;
				} 
				else ErrorMessage("Ridder", "never get here.");
			
				if (std::fabs(xh-xl) <= xacc) 
					return ans;
			}
		
			ErrorMessage("Ridder", "zriddr exceed maximum iterations");
		}
		else 
		{
			if (fl == 0.) return x1;
			if (fh == 0.) return x2;
			
			ErrorMessage("Ridder", "root must be bracketed in zriddr.");
		}
	
		return 0.0;
	}
	
	double ThermodynamicsMap_CHEMKIN::NewtonRaphson(const double H, const double P_Pa, const double* x, const double x1, const double x2, const double xacc)
	{
		const int JMAX=20;
		int j;
		double df,dx,f,rtn;

		rtn=0.5*(x1+x2);
		for (j=0;j<JMAX;j++) 
		{
			Function(H, P_Pa, x, rtn, f, df);
			dx=f/df;
			rtn -= dx;
			if ((x1-rtn)*(rtn-x2) < 0.)
				ErrorMessage("NewtonRaphson", "Jumped out of brackets in rtnewt");
			if (std::fabs(dx) < xacc) 
			return rtn;
		}
	
		ErrorMessage("NewtonRaphson", "Maximum number of iterations exceeded in rtnewt");
		return 0.;
	}

	double ThermodynamicsMap_CHEMKIN::Function(const double H, const double P_Pa, const double* x, const double T)
	{

		SetTemperature(T);
		SetPressure(P_Pa);
		const double HCandidate = hMolar_Mixture_From_MoleFractions(x);
		return (H - HCandidate);
	}
	
	void ThermodynamicsMap_CHEMKIN::Function(const double H, const double P_Pa, const double* x, const double T, double& f, double& df)
	{
		SetTemperature(T);
		SetPressure(P_Pa);
		df = cpMolar_Mixture_From_MoleFractions(x);
		f = hMolar_Mixture_From_MoleFractions(x);
		df *= -1;
		f = H-f;
	}	

	void ThermodynamicsMap_CHEMKIN::NASA_LowT(const unsigned int i, double* coefficients) const
	{
		unsigned int k = i*5;
		unsigned int j = i*6;
		coefficients[0] = Cp_LT[k++];
		coefficients[1] = Cp_LT[k++];
		coefficients[2] = Cp_LT[k++];
		coefficients[3] = Cp_LT[k++];
		coefficients[4] = Cp_LT[k++];
		coefficients[5] = DH_LT[j+5];
		coefficients[6] = DS_LT[j+5];
	}

	void ThermodynamicsMap_CHEMKIN::NASA_HighT(const unsigned int i, double* coefficients) const
	{
		unsigned int k = i*5;
		unsigned int j = i*6;

		coefficients[0] = Cp_HT[k++];
		coefficients[1] = Cp_HT[k++];
		coefficients[2] = Cp_HT[k++];
		coefficients[3] = Cp_HT[k++];
		coefficients[4] = Cp_HT[k++];
		coefficients[5] = DH_HT[j+5];
		coefficients[6] = DS_HT[j+5];
	}
	
	void ThermodynamicsMap_CHEMKIN::NASA_Coefficients(const unsigned int i, double* coefficients) const
	{
		double sub_coefficients[7];

		// Intermediate temperature
		coefficients[0] = TM[i];

		// Low-T coefficients
		NASA_LowT(i, sub_coefficients);
		for (unsigned int j = 0; j < 7; j++)
			coefficients[1+j] = sub_coefficients[j];	

		// High-T coefficients
		NASA_HighT(i, sub_coefficients);
		for (unsigned int j = 0; j < 7; j++)
			coefficients[8+j] = sub_coefficients[j];	
		
	}

	void ThermodynamicsMap_CHEMKIN::NASA_Coefficients(double* coefficients) const
	{
		double sub_coefficients[15];
		
		for(unsigned int i=0;i<this->nspecies_;i++)
		{
			NASA_Coefficients(i, sub_coefficients);
			for (unsigned int j = 0; j < 15; j++)
				coefficients[j*this->nspecies_ + i] = sub_coefficients[j];
		}
	}
}

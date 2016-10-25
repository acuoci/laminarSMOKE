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
	ThermodynamicsMap_CHEMKIN<map>::ThermodynamicsMap_CHEMKIN(const unsigned int nSpecies, const unsigned int nPoints)
	{
		this->nspecies_ = nSpecies;
		this->npoints_ = nPoints;
                
                this->verbose_output_ = true;
                
		MemoryAllocation();
	}

	template<typename map> 
	ThermodynamicsMap_CHEMKIN<map>::ThermodynamicsMap_CHEMKIN(rapidxml::xml_document<>& doc, const unsigned int nPoints)
	{
		this->npoints_ = nPoints;
                
                this->verbose_output_ = true;
	
		ImportSpeciesFromXMLFile(doc);
		this->ImportElementsFromXMLFile(doc);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(doc);
	}
        
    template<typename map> 
	ThermodynamicsMap_CHEMKIN<map>::ThermodynamicsMap_CHEMKIN(rapidxml::xml_document<>& doc, bool verbose, const unsigned int nPoints)
	{
		this->npoints_ = nPoints;
                
                this->verbose_output_ = verbose;
	
		ImportSpeciesFromXMLFile(doc);
		this->ImportElementsFromXMLFile(doc);
		MemoryAllocation();
		ImportCoefficientsFromXMLFile(doc);
	}
    
        template<typename map> 
	ThermodynamicsMap_CHEMKIN<map>::ThermodynamicsMap_CHEMKIN( const ThermodynamicsMap_CHEMKIN<map>& rhs )
        {
            CopyFromMap(rhs);
        }
        
    
	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::CopyFromMap( const ThermodynamicsMap_CHEMKIN<map>& rhs )
        {
            this->nspecies_ = rhs.nspecies_;
            this->npoints_  = rhs.npoints_;
            
            MemoryAllocation();
             
	    this->names_ =  rhs.names_;
            this->MW_    =  rhs.MW_;
            this->uMW_   =  rhs.uMW_;

            this->atomic_composition_ =  rhs.atomic_composition_;
	    this->elements_ = rhs.elements_; 
	    this->verbose_output_          = rhs.verbose_output_;
            
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
	 
	template<typename map> 
	ThermodynamicsMap_CHEMKIN<map>::~ThermodynamicsMap_CHEMKIN(void)
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

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::MemoryAllocation()
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
	void ThermodynamicsMap_CHEMKIN<map>::SetTemperature(const map& T)
	{
		this->T_ = T;
		cp_must_be_recalculated_ = true;
		h_must_be_recalculated_ = true;
		s_must_be_recalculated_ = true;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::SetPressure(const map& P)
	{
		this->P_ = P;
	}

	
	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::SetCoefficients(const unsigned k, const double* coefficients)
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
	void ThermodynamicsMap_CHEMKIN<map>::ImportCoefficientsFromASCIIFile(std::ifstream& fInput)
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

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::ImportSpeciesFromXMLFile(rapidxml::xml_document<>& doc)
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
			ErrorMessage("ThermodynamicsMap_CHEMKIN<map>::ImportSpeciesFromXMLFile", "Error in reading the list of species.");
		}
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::ImportCoefficientsFromXMLFile(rapidxml::xml_document<>& doc)
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
				ErrorMessage("ThermodynamicsMap_CHEMKIN<map>::ImportCoefficientsFromXMLFile", "Thermodynamics type not supported!"); 
		}
		else
			ErrorMessage("ThermodynamicsMap_CHEMKIN<map>::ImportCoefficientsFromXMLFile", "Thermodynamics tag was not found!");
	}

	template<typename map> 
	inline void ThermodynamicsMap_CHEMKIN<map>::MolecularWeight_From_MoleFractions(map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		MW = Dot(x,this->MW_);
	}

	template<typename map> 
	inline void ThermodynamicsMap_CHEMKIN<map>::MolecularWeight_From_MassFractions(map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		MW = 1./Dot(y,this->uMW_);
	}

	template<typename map> 
	inline void ThermodynamicsMap_CHEMKIN<map>::MassFractions_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& y, map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		MW = Dot(x,this->MW_);
		ElementByElementProduct(x, this->MW_, &y);
		Product(1./MW, &y);
	}
	
	template<typename map> 
	inline void ThermodynamicsMap_CHEMKIN<map>::MoleFractions_From_MassFractions(OpenSMOKE::OpenSMOKEVectorDouble& x, map& MW, const OpenSMOKE::OpenSMOKEVectorDouble& y)
	{
		MW = 1./Dot(y,this->uMW_);
		Product(MW, y, &x);
		ElementByElementProduct(x,this->uMW_,&x);
	}

	template<> 
	inline void ThermodynamicsMap_CHEMKIN<double>::cp_over_R()
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
	inline void ThermodynamicsMap_CHEMKIN<map>::cp_over_R()
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
	inline void ThermodynamicsMap_CHEMKIN<double>::h_over_RT()
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
	inline void ThermodynamicsMap_CHEMKIN<map>::h_over_RT()
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
	inline void ThermodynamicsMap_CHEMKIN<double>::s_over_R()
	{
		if (s_must_be_recalculated_ == true)
		{
			const double logT = std::log(this->T_);
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
	inline void ThermodynamicsMap_CHEMKIN<double>::g_over_RT()
	{
		h_over_RT();
		s_over_R();

		Sub(species_h_over_RT_, species_s_over_R_, &species_g_over_RT_);
	}
	
	template<typename map> 
	inline void ThermodynamicsMap_CHEMKIN<map>::s_over_R()
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
				const double logT1 = std::log(T1);

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
	void ThermodynamicsMap_CHEMKIN<map>::cpMolar_Mixture_From_MoleFractions(map& cpmix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		cp_over_R();
		cpmix = Dot(species_cp_over_R_, x);
		cpmix *= PhysicalConstants::R_J_kmol;
	}
	
	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::hMolar_Mixture_From_MoleFractions(map& hmix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		h_over_RT();
		hmix = Dot(species_h_over_RT_, x);
		hmix *= PhysicalConstants::R_J_kmol*this->T_;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::sMolar_Mixture_From_MoleFractions(map& smix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		s_over_R();

		const double eps = 1.e-32;
		double sum=0.;
		for(unsigned int i=1;i<=this->nspecies_;i++)
			sum += x[i]*std::log(x[i]+eps);
		
		smix = Dot(species_s_over_R_, x) - std::log(this->P_/101325.) - sum;
		smix *= PhysicalConstants::R_J_kmol;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::uMolar_Mixture_From_MoleFractions(map& umix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		h_over_RT();
		umix = Dot(species_h_over_RT_, x) - 1.;
		umix *= PhysicalConstants::R_J_kmol*this->T_;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::gMolar_Mixture_From_MoleFractions(map& gmix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		double hmix, smix;
		hMolar_Mixture_From_MoleFractions(hmix,x);
		sMolar_Mixture_From_MoleFractions(smix,x);
		gmix = hmix-this->T_*smix;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::aMolar_Mixture_From_MoleFractions(map& amix, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		double umix, smix;
		uMolar_Mixture_From_MoleFractions(umix,x);
		sMolar_Mixture_From_MoleFractions(smix,x);
		amix = umix-this->T_*smix;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::cpMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& cp_species)
	{
		cp_over_R();
		Product(PhysicalConstants::R_J_kmol, species_cp_over_R_, &cp_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::hMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& h_species)
	{
		h_over_RT();
		Product(PhysicalConstants::R_J_kmol*this->T_, species_h_over_RT_, &h_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::sMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& s_species)
	{
		s_over_R();
		Product(PhysicalConstants::R_J_kmol, species_s_over_R_, &s_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::sMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& s_species, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		const double eps = 1.e-32;
		const double log_P_pver_Patm = std::log(this->P_/101325.);
		s_over_R();
		for(unsigned int i=1;i<=this->nspecies_;i++)
			s_species[i] = (x[i]>eps) ? species_s_over_R_[i] - std::log(x[i]+eps) - log_P_pver_Patm : species_s_over_R_[i] - log_P_pver_Patm;
		Product(PhysicalConstants::R_J_kmol, &s_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::uMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& u_species)
	{
		h_over_RT();
		Add(species_h_over_RT_, -1., &u_species);
		Product(PhysicalConstants::R_J_kmol*this->T_, &u_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::gMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& g_species)
	{
		h_over_RT();
		s_over_R();
		Sub(species_h_over_RT_, species_s_over_R_, &g_species);
		Product(PhysicalConstants::R_J_kmol*this->T_, &g_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::gMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& g_species, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		const double eps = 1.e-32;
		const double log_P_pver_Patm = std::log(this->P_/101325.);
		s_over_R();
		h_over_RT();
		for(unsigned int i=1;i<=this->nspecies_;i++)
			g_species[i] = (x[i]>eps) ? species_h_over_RT_[i] - (species_s_over_R_[i] - std::log(x[i]+eps) - log_P_pver_Patm) : 
			                            species_h_over_RT_[i] - (species_s_over_R_[i] - log_P_pver_Patm);
		Product(PhysicalConstants::R_J_kmol*this->T_, &g_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::aMolar_Species(OpenSMOKE::OpenSMOKEVectorDouble& a_species)
	{
		h_over_RT();
		s_over_R();
		Sub(species_h_over_RT_, species_s_over_R_, &a_species);
		a_species -= 1.;
		Product(PhysicalConstants::R_J_kmol*this->T_, &a_species);
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::aMolar_Species_MixtureAveraged_From_MoleFractions(OpenSMOKE::OpenSMOKEVectorDouble& a_species, const OpenSMOKE::OpenSMOKEVectorDouble& x)
	{
		gMolar_Species_MixtureAveraged_From_MoleFractions(a_species, x);
		a_species -= PhysicalConstants::R_J_kmol*this->T_;
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::DerivativesOfConcentrationsWithRespectToMassFractions(const double cTot, const double MW, const OpenSMOKE::OpenSMOKEVectorDouble& omega, OpenSMOKE::OpenSMOKEMatrixDouble* dc_over_omega)
	{
		for(unsigned int j=1;j<=this->nspecies_;j++)
		{
			const double coefficient = cTot*MW/this->MW_[j];
			for(unsigned int i=1;i<=this->nspecies_;i++)
				(*dc_over_omega)[j][i] = -coefficient*MW/this->MW_[i]*omega[j];
			(*dc_over_omega)[j][j] +=coefficient;
		}
	}

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::Test(const int nLoops, const map& T, int* index)
	{
	}
	
	template<typename map> 
	double ThermodynamicsMap_CHEMKIN<map>::GetTemperatureFromEnthalpyAndMoleFractions(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double TFirstGuess)
	{
		const double diff_temperature_to_switch_to_newtons = 15.;
		const unsigned int max_number_of_enlargements = 10;
		
		// Check how good is the first guess temperature
		double absolute_error;
		{
			double CpFirstGuess;
			double HFirstGuess;
			SetTemperature(TFirstGuess);
			SetPressure(P_Pa);
			cpMolar_Mixture_From_MoleFractions(CpFirstGuess, x);
			hMolar_Mixture_From_MoleFractions(HFirstGuess, x);
			
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
	
	template<typename map> 
	bool ThermodynamicsMap_CHEMKIN<map>::FindTheRightInterval(	const double H, const OpenSMOKEVectorDouble& x, const double TFirstGuess,
																double& TA, double&TB, double& HA, double& HB,
																const double Diff_Temperature)
	{
		const double TMinimum = 250.;
		const double TMaximum = 6000.;
		const unsigned int max_number_of_enlargements = 50;

		TA = TB = TFirstGuess;

		SetTemperature(TA);
		hMolar_Mixture_From_MoleFractions(HA, x);
		
		SetTemperature(TB);
		hMolar_Mixture_From_MoleFractions(HB, x);

		unsigned int count = 0;
		for(;;)
		{
			if ((HA-H)*(HB-H)<=0.)
				return true;

			if ((HA-H)>0.)
			{	
				TA = std::max(TMinimum, TA-Diff_Temperature);
				SetTemperature(TA);
				hMolar_Mixture_From_MoleFractions(HA, x);
			}

			if ((HB-H)<0.)
			{	
				TB = std::min(TMaximum, TB+Diff_Temperature);
				SetTemperature(TB);
				hMolar_Mixture_From_MoleFractions(HB, x);
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

	template<typename map> 
	double ThermodynamicsMap_CHEMKIN<map>::Brent(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double x1, const double x2, const double tol)
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
		return 0.0;
	}

	template<typename map> 
	double ThermodynamicsMap_CHEMKIN<map>::Ridder(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double x1, const double x2, const double xacc)
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
			if (fl == 0.0) return x1;
			if (fh == 0.0) return x2;
			
			ErrorMessage("Ridder", "root must be bracketed in zriddr.");
		}
	
		return 0.0;
	}
	
	template<typename map> 
	double ThermodynamicsMap_CHEMKIN<map>::NewtonRaphson(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double x1, const double x2, const double xacc)
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
			if ((x1-rtn)*(rtn-x2) < 0.0)
				ErrorMessage("NewtonRaphson", "Jumped out of brackets in rtnewt");
			if (std::fabs(dx) < xacc) 
			return rtn;
		}
	
		ErrorMessage("NewtonRaphson", "Maximum number of iterations exceeded in rtnewt");
		return 0.0;
	}

	template<typename map> 
	double ThermodynamicsMap_CHEMKIN<map>::Function(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double T)
	{
		double HCandidate;
		SetTemperature(T);
		SetPressure(P_Pa);
		hMolar_Mixture_From_MoleFractions(HCandidate, x);
		return H-HCandidate;
	}
	
	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::Function(const double H, const double P_Pa, const OpenSMOKEVectorDouble& x, const double T, double& f, double& df)
	{
		SetTemperature(T);
		SetPressure(P_Pa);
		cpMolar_Mixture_From_MoleFractions(df, x);
		hMolar_Mixture_From_MoleFractions(f, x);
		df *= -1;
		f = H-f;
	}	

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::NASA_LowT(const unsigned int i, double* coefficients) const
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

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::NASA_HighT(const unsigned int i, double* coefficients) const
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
	
	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::NASA_Coefficients(const unsigned int i, double* coefficients) const
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

	template<typename map> 
	void ThermodynamicsMap_CHEMKIN<map>::NASA_Coefficients(double* coefficients) const
	{
		double sub_coefficients[15];
		
		for(unsigned int i=0;i<this->nspecies_;i++)
		{
			NASA_Coefficients(i, sub_coefficients);
			for (unsigned int j = 0; j < 15; j++)
				coefficients[j*this->nspecies_ + i] = sub_coefficients[j];
		}
	}

	/*
	template<> 
	void ThermodynamicsMap_CHEMKIN<double>::Test(const int nLoops, const double& T, int* index)
	{
		// Loops
		int speciesLoopsThermo	= nLoops*100;
		int mixtureLoopsThermo	= nLoops*100;

		double cpmix, hmix, smix;

		// Composition (mass fractions)
		OpenSMOKEVectorDouble omega(this->nspecies_);
		for(int i=1;i<=this->nspecies_;i++)
			omega[i] = 1./double(this->nspecies_);


		// Times
		double speciesThermoTime, mixtureThermoTime;

		{
			double tStart = OpenSMOKEGetCpuTime();
			for(int k=1;k<=speciesLoopsThermo;k++)
			{
				cp(T);
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesThermoTime = tEnd-tStart;
			std::cout << "Species Cp " << "Time: " << tEnd - tStart << std::endl;
		}

		{
			double tStart = OpenSMOKEGetCpuTime();
			for(int k=1;k<=speciesLoopsThermo;k++)
			{
				h(T);
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesThermoTime += tEnd-tStart;
			std::cout << "Species H " << "Time: " << tEnd - tStart << std::endl;
		}
		
		{
			double tStart = OpenSMOKEGetCpuTime();
			for(int k=1;k<=speciesLoopsThermo;k++)
			{
				s(T);
			}
			double tEnd = OpenSMOKEGetCpuTime();
			speciesThermoTime += tEnd-tStart;
			std::cout << "Species S "<< "Time: " << tEnd - tStart << std::endl;
		}

		{
			double tStart = OpenSMOKEGetCpuTime();
			for(int k=1;k<=mixtureLoopsThermo;k++)
			{
				cpMix(cpmix, omega);
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureThermoTime = tEnd-tStart;
			std::cout << "Mixture Cp " << "Time: " << tEnd - tStart << std::endl;
		}

		{
			double tStart = OpenSMOKEGetCpuTime();
			for(int k=1;k<=mixtureLoopsThermo;k++)
			{
				hMix(hmix, omega);
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureThermoTime += tEnd-tStart;
			std::cout << "Mixture H " << "Time: " << tEnd - tStart << std::endl;
		}
		
		{
			double tStart = OpenSMOKEGetCpuTime();
			for(int k=1;k<=mixtureLoopsThermo;k++)
			{
				sMix(smix, omega);
			}
			double tEnd = OpenSMOKEGetCpuTime();
			mixtureThermoTime += tEnd-tStart;
			std::cout << "Mixture S "<< "Time: " << tEnd - tStart << std::endl;
		}

	std::ofstream fBenchmark;
	fBenchmark.open("Benchmark.plus", std::ios::out);
	fBenchmark.setf(std::ios::scientific);
	
	// Thermodynamics
	fBenchmark << cpSpecies_[index[1]] << " " << cpSpecies_[index[2]] << " " << cpSpecies_[index[3]] << std::endl; 
	fBenchmark << hSpecies_[index[1]] << " " << hSpecies_[index[2]] << " " << hSpecies_[index[3]] << std::endl; 
	fBenchmark << sSpecies_[index[1]] << " " << sSpecies_[index[2]] << " " << sSpecies_[index[3]] << std::endl; 
	fBenchmark << cpmix << std::endl;
	fBenchmark << hmix << std::endl;
	fBenchmark << smix << std::endl;
	fBenchmark << speciesLoopsThermo	<< " " << speciesThermoTime << " " << speciesThermoTime/double(speciesLoopsThermo)*1000.	<< std::endl;
	fBenchmark << mixtureLoopsThermo	<< " " << mixtureThermoTime << " " << mixtureThermoTime/double(mixtureLoopsThermo)*1000.	<< std::endl;
	
	fBenchmark.close();


	}
	*/

}

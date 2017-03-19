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
#include <numeric>

namespace OpenSMOKE
{
	unsigned int ThermodynamicsMap::IndexOfSpecies(const std::string name) const
	{
		for(unsigned int i=0;i<nspecies_;++i)
			if (name == names_[i])
				return i+1;

		ErrorMessage(	"const unsigned int ThermodynamicsMap::IndexOfSpecies(const std::string name) const", 
						"The requested species " + name + " is not available in the kinetic mechanism");
		
		return 0;
	}
 
	unsigned int ThermodynamicsMap::IndexOfSpeciesWithoutError(const std::string name) const
	{
		for(unsigned int i=0;i<nspecies_;++i)
			if (name == names_[i])
				return i+1;

		return 0;
	}

	unsigned int ThermodynamicsMap::IndexOfElement(const std::string name) const
	{
		for(unsigned int i=0;i<elements_.size();++i)
			if (name == elements_[i])
				return i+1;

		ErrorMessage(	"const unsigned int ThermodynamicsMap::IndexOfElement(const std::string name) const", 
						"The requested element " + name + "is not available in the kinetic mechanism");
		
		return 0;
	}

	unsigned int ThermodynamicsMap::IndexOfElementWithoutError(const std::string name) const
	{
		for(unsigned int i=0;i<elements_.size();++i)
			if (name == elements_[i])
				return i+1;

		return 0;
	}

	void ThermodynamicsMap::ImportElementsFromXMLFile(rapidxml::xml_document<>& doc)
	{
		rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

		rapidxml::xml_node<>* node_number = opensmoke_node->first_node("NumberOfElements");
		std::stringstream values_number(node_number->value());
		unsigned int number_of_elements = atoi(values_number.str().c_str());
			
		{
			this->elements_.resize(number_of_elements);
			rapidxml::xml_node<>* node_names = opensmoke_node->first_node("NamesOfElements");
			std::stringstream values_names(node_names->value());
			for(unsigned int j=0;j<number_of_elements;j++)
				values_names >> this->elements_[j];
				
			this->mw_elements__.resize(number_of_elements);
			for(unsigned int j=0;j<number_of_elements;j++)
				this->mw_elements__[j] = OpenSMOKE::AtomicWeights[this->elements_[j]];
		}

		{
			this->atomic_composition_.resize(this->nspecies_, number_of_elements);
			rapidxml::xml_node<>* node_composition = opensmoke_node->first_node("AtomicComposition");
			std::stringstream values_composition(node_composition->value());
			for(unsigned int i=0;i<this->nspecies_;i++)
				for(unsigned int j=0;j<number_of_elements;j++)
					values_composition >> this->atomic_composition_(i,j);
		}
	}
 
	std::vector<double> ThermodynamicsMap::GetMoleFractionsFromEquivalenceRatio(
		const double equivalence_ratio,
		const std::vector<std::string>& fuel_names, const std::vector<double>& moles_fuel,
		const std::vector<std::string>& oxidizer_names,	const std::vector<double>& moles_oxidizer)
	{
		const std::size_t number_of_fuels	= moles_fuel.size();
		const std::size_t number_of_oxidizers	= moles_oxidizer.size();

		std::vector<double> nC(number_of_fuels);
		std::vector<double> nH(number_of_fuels);
		std::vector<double> nO(number_of_fuels);
		std::vector<double> nOxidizers(number_of_oxidizers);

		const unsigned int jC = IndexOfElementWithoutError("C");
		const unsigned int jO = IndexOfElementWithoutError("O");
		const unsigned int jH = IndexOfElementWithoutError("H");
	
		if (jC>0)
			for(unsigned int j=0;j<number_of_fuels;j++)
				nC[j] = atomic_composition_(IndexOfSpecies(fuel_names[j])-1,jC-1);

		if (jH>0)
			for(unsigned int j=0;j<number_of_fuels;j++)
				nH[j]	= atomic_composition_(IndexOfSpecies(fuel_names[j])-1,jH-1);

		if (jO>0)
			for(unsigned int j=0;j<number_of_fuels;j++)
				nO[j]	= atomic_composition_(IndexOfSpecies(fuel_names[j])-1,jO-1);

		const double nFuel	= equivalence_ratio * std::accumulate(moles_fuel.begin(),moles_fuel.end(),0.);
		
		unsigned int jO2 = 0;
		for(unsigned int j=0;j<number_of_oxidizers;j++)
			if (oxidizer_names[j] == "O2")	
			{
				jO2 = j;
				break;
			}

		const double nO2		= (	2.*std::inner_product(nC.begin(), nC.end(), moles_fuel.begin(), 0.) +
									0.50*std::inner_product(nH.begin(), nH.end(), moles_fuel.begin(), 0.) - 
									std::inner_product(nO.begin(), nO.end(), moles_fuel.begin(), 0.) ) / 2.;
		
		for(unsigned int j=0;j<number_of_oxidizers;j++)
			nOxidizers[j] = moles_oxidizer[j]/moles_oxidizer[jO2]*nO2;
	
		const double n = nFuel + std::accumulate(nOxidizers.begin(),nOxidizers.end(),0.);

		std::vector<double> mole(NumberOfSpecies());	
		for(unsigned int j=0;j<number_of_fuels;j++)
			mole[IndexOfSpecies(fuel_names[j])-1] += equivalence_ratio*moles_fuel[j]/n;
		for(unsigned int j=0;j<number_of_oxidizers;j++)
			mole[IndexOfSpecies(oxidizer_names[j])-1] += nOxidizers[j]/n;

		return mole;
	}

	double ThermodynamicsMap::GetLocalEquivalenceRatio( 	const std::vector<double>& moles, const std::vector<double>& moles_st,
									const std::vector<std::string>& fuel_names)
	{
		double nFuelStoichiometric = 0.;
		double nFuel = 0.;
		for(unsigned int j=0;j<fuel_names.size();j++)
		{
			nFuelStoichiometric += moles_st[IndexOfSpecies(fuel_names[j])-1];
			nFuel               += moles[IndexOfSpecies(fuel_names[j])-1];
		}

		const double nOxygenStoichiometric = moles_st[IndexOfSpecies("O2")-1];
		const double nOxygen               = moles[IndexOfSpecies("O2")-1];

		const double phi = Min( nFuel*nOxygenStoichiometric/nFuelStoichiometric/(nOxygen+1e-12), 100.);

		return phi;
	}

	double 	ThermodynamicsMap::GetMixtureFractionFromMassFractions(
		const std::vector<double>& mass,
		const std::vector<std::string>& fuel_names, 	const std::vector<double>& mass_fuel,
		const std::vector<std::string>& oxidizer_names,	const std::vector<double>& mass_oxidizer)
	{
		const unsigned int number_of_fuels	= mass_fuel.size();
		const unsigned int number_of_oxidizers	= mass_oxidizer.size();

		double omegaFuelC = 0.;
		double omegaFuelO = 0.;
		double omegaFuelH = 0.;
		double omegaOxidizerC = 0.;
		double omegaOxidizerO = 0.;
		double omegaOxidizerH = 0.;
		double omegaC = 0.;
		double omegaO = 0.;
		double omegaH = 0.;

		const unsigned int jC = IndexOfElementWithoutError("C");
		const unsigned int jO = IndexOfElementWithoutError("O");
		const unsigned int jH = IndexOfElementWithoutError("H");
	
		if (jC>0)
		{
			for(unsigned int j=0;j<number_of_fuels;j++)
				omegaFuelC += atomic_composition_(IndexOfSpecies(fuel_names[j])-1,jC-1)*mass_fuel[j];
			for(unsigned int j=0;j<number_of_oxidizers;j++)
				omegaOxidizerC += atomic_composition_(IndexOfSpecies(oxidizer_names[j])-1,jC-1)*mass_oxidizer[j];
			for(unsigned int j=0;j<nspecies_;j++)
				omegaC += atomic_composition_(j,jC-1)*mass[j];
		}
		
		if (jH>0)
		{
			for(unsigned int j=0;j<number_of_fuels;j++)
				omegaFuelH += atomic_composition_(IndexOfSpecies(fuel_names[j])-1,jH-1)*mass_fuel[j];
			for(unsigned int j=0;j<number_of_oxidizers;j++)
				omegaOxidizerH += atomic_composition_(IndexOfSpecies(oxidizer_names[j])-1,jH-1)*mass_oxidizer[j];
			for(unsigned int j=0;j<nspecies_;j++)
				omegaH += atomic_composition_(j,jH-1)*mass[j];
		}

		if (jO>0)
		{
			for(unsigned int j=0;j<number_of_fuels;j++)
				omegaFuelO += atomic_composition_(IndexOfSpecies(fuel_names[j])-1,jO-1)*mass_fuel[j];
			for(unsigned int j=0;j<number_of_oxidizers;j++)
				omegaOxidizerO += atomic_composition_(IndexOfSpecies(oxidizer_names[j])-1,jO-1)*mass_oxidizer[j];
			for(unsigned int j=0;j<nspecies_;j++)
				omegaO += atomic_composition_(j,jO-1)*mass[j];
		}

		const double WC = OpenSMOKE::AtomicWeights["C"];
		const double WO = OpenSMOKE::AtomicWeights["O"];
		const double WH = OpenSMOKE::AtomicWeights["H"];

		const double z = ( 	2.*(omegaC-omegaOxidizerC)/WC 		+ (omegaH-omegaOxidizerH)/2./WH 	- (omegaO-omegaOxidizerO)/WO 		) /
				 ( 	2.*(omegaFuelC-omegaOxidizerC)/WC 	+ (omegaFuelH-omegaOxidizerH)/2./WH 	- (omegaFuelO-omegaOxidizerO)/WO 	) ;

		return z;
	}
}

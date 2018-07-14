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


#ifndef OpenSMOKE_Grammar_VirtualChemistry_H
#define OpenSMOKE_Grammar_VirtualChemistry_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_VirtualChemistry : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Table",
				OpenSMOKE::SINGLE_PATH,
				"Path to the look-up table",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Version",
				OpenSMOKE::SINGLE_INT,
				"Optimization table version: 170911 | 171013",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Reactions",
				OpenSMOKE::SINGLE_BOOL,
				"Reactions on/off (default: true)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelName",
				OpenSMOKE::SINGLE_STRING,
				"Fuel name",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerName",
				OpenSMOKE::SINGLE_STRING,
				"Oxidizer name",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InertName",
				OpenSMOKE::SINGLE_STRING,
				"Inert name",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FuelMW",
				OpenSMOKE::SINGLE_MEASURE,
				"Fuel molecular weight",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidizerMW",
				OpenSMOKE::SINGLE_MEASURE,
				"Oxidizer molecular weight",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InertMW",
				OpenSMOKE::SINGLE_MEASURE,
				"Inert molecular weight",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Viscosity_mu0",
				OpenSMOKE::SINGLE_MEASURE,
				"Viscosity at reference temperature T0",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Viscosity_T0",
				OpenSMOKE::SINGLE_MEASURE,
				"Reference temperature T0 for reference viscosity mu0",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Viscosity_Beta0",
				OpenSMOKE::SINGLE_DOUBLE,
				"Exponent for viscosity correlation",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Viscosity_Pr0",
				OpenSMOKE::SINGLE_DOUBLE,
				"Prandtl number of viscosity correlation",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SubMechanism_CO",
				OpenSMOKE::SINGLE_BOOL,
				"Sub-mechanism for CO formation (default: false)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SubMechanism_NO",
				OpenSMOKE::SINGLE_BOOL,
				"Sub-mechanism for NO formation (default: false)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Table_CO",
				OpenSMOKE::SINGLE_PATH,
				"Path to the look-up table for CO sub-mechanism",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Table_NO",
				OpenSMOKE::SINGLE_PATH,
				"Path to the look-up table for NO sub-mechanism",
				false));

			/*
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Regressions",
				OpenSMOKE::VECTOR_STRING,
				"Regressions on parameters (default: false)",
				false));
			*/
		}
	};
}

#endif /* OpenSMOKE_Grammar_VirtualChemistry_H */

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


#ifndef OpenSMOKE_Grammar_HMOM_H
#define OpenSMOKE_Grammar_HMOM_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_HMOM : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@HMOM",
				OpenSMOKE::SINGLE_BOOL,
				"Hybrid Method of Moments: on/off (default: true)",
				true));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NucleationModel",
				OpenSMOKE::SINGLE_INT,
				"Nucleation model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceGrowthModel",
				OpenSMOKE::SINGLE_INT,
				"Surface growth model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OxidationModel",
				OpenSMOKE::SINGLE_INT,
				"Oxidation model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CondensationModel",
				OpenSMOKE::SINGLE_INT,
				"Condensation model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CoagulationModel",
				OpenSMOKE::SINGLE_INT,
				"Coagulation model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ContinousCoagulationModel",
				OpenSMOKE::SINGLE_INT,
				"Continous coagulation model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThermophoreticModel",
				OpenSMOKE::SINGLE_INT,
				"Thermophoretic model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FractalDiameterModel",
				OpenSMOKE::SINGLE_INT,
				"Fractal diameter model: 0=off, 1=on (default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@CollisionDiameterModel",
				OpenSMOKE::SINGLE_INT,
				"Collision diameter model: 1 or 2 (default: 2)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NumberOfCarbonPAH",
				OpenSMOKE::SINGLE_INT,
				"Number of carbon atoms in PAHs (default 16)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PAH",
				OpenSMOKE::VECTOR_STRING,
				"Species to be assumed as PAH (example: @PAH C10H8 BIN1A;)",
				false));
				
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PAHConsumption",
				OpenSMOKE::SINGLE_BOOL,
				"Consumption of PAH is accounted for (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RadiativeHeatTransfer",
				OpenSMOKE::SINGLE_BOOL,
				"Radiative heat transfer (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PlanckCoefficient",
				OpenSMOKE::SINGLE_STRING,
				"Calculation of Planck Coefficient: Smooke (default) | Kent | Sazhin | none",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SchmidtNumber",
				OpenSMOKE::SINGLE_DOUBLE,
				"Schmidt number (default: 50)",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_HMOM_H */
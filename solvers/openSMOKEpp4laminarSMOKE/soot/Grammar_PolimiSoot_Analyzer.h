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


#ifndef OpenSMOKE_Grammar_PolimiSoot_Analyzer_H
#define OpenSMOKE_Grammar_PolimiSoot_Analyzer_H

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OpenSMOKE
{
	class Grammar_PolimiSoot_Analyzer : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SootLabel",
				OpenSMOKE::SINGLE_STRING,
				"Label indicating the soot particles (default: BIN)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@FractalDimension",
				OpenSMOKE::SINGLE_DOUBLE,
				"Fractal dimension of soot particles (default: 1.8)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SootMinimumSection",
				OpenSMOKE::SINGLE_INT,
				"Index of minimum discrete section which is considered soot (default: 5)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SootMinimumFractalDimension",
				OpenSMOKE::SINGLE_INT,
				"Index of minimum discrete section from which to apply the fractal dimension (default: 12)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Density",
				OpenSMOKE::VECTOR_DOUBLE,
				"Linear density: BINstart DENSITYstart BINend DENSITYend (default: 10 1500. 20 1700.)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PhysicalDiffusion",
				OpenSMOKE::SINGLE_BOOL,
				"Physical diffusion for soot particles (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PhysicalDiffusionReductionCoefficient",
				OpenSMOKE::SINGLE_DOUBLE,
				"Physical diffusion reduction coefficient for soot particles(default: 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThermophoreticEffect",
				OpenSMOKE::SINGLE_BOOL,
				"Thermophoretic effect for soot sections (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThermophoreticEffectAmplificationFactor",
				OpenSMOKE::SINGLE_DOUBLE,
				"Thermoforetic effect amplification factor",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThermophoreticEffectInCorrectionVelocity",
				OpenSMOKE::SINGLE_BOOL,
				"Thermophoretic effect for soot sections in calculation of correction velocity (default: false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThermophoreticEffectInEnthalpyFluxes",
				OpenSMOKE::SINGLE_STRING,
				"Thermophoretic effect for soot sections in calculation of enthalpy fluxes: DoNotExclude | Exclude | ExcludeOnlyThermophoreticEffect (default: DoNotExclude)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RadiativeHeatTransfer",
				OpenSMOKE::SINGLE_BOOL,
				"Radiative heat transfer (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThermophoreticEffectSmoothingTime",
				OpenSMOKE::SINGLE_MEASURE,
				"Smoothing time for thermophoretic effect (default: 0.)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@PlanckCoefficient",
				OpenSMOKE::SINGLE_STRING,
				"Calculation of Planck Coefficient: Smooke (default) | Kent | Sazhin | none",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@WritePSDF",
				OpenSMOKE::SINGLE_BOOL,
				"Write the soot particle size distribution function (PSDF) on file (default: true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@ThresholdForPSDF",
				OpenSMOKE::SINGLE_DOUBLE,
				"Threshold (soot volume fraction) for writing the PSDF on file (default: 1e-11)",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_PolimiSoot_Analyzer_H */

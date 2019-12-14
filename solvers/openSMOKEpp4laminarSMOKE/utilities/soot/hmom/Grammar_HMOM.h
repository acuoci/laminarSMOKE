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
|	License                                                           |
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

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@StickingCoefficient",
				OpenSMOKE::SINGLE_DOUBLE,
				"StickingCoefficient (default: 0.002)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SootDensity",
				OpenSMOKE::SINGLE_MEASURE,
				"Density of soot particles (default: 1800 kg/m3)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDensity",
				OpenSMOKE::SINGLE_MEASURE,
				"Surface density (default: 1.7 #/m2)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDensityCorrectionCoefficient",
				OpenSMOKE::SINGLE_BOOL,
				"SurfaceDensityCorrectionCoefficient (default: false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDensityCorrectionCoefficientA1",
				OpenSMOKE::SINGLE_DOUBLE,
				"SurfaceDensityCorrectionCoefficientA1 (default: 12.65)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDensityCorrectionCoefficientA2",
				OpenSMOKE::SINGLE_DOUBLE,
				"SurfaceDensityCorrectionCoefficientA2 (default: -0.00563)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDensityCorrectionCoefficientB1",
				OpenSMOKE::SINGLE_DOUBLE,
				"SurfaceDensityCorrectionCoefficientB1 (default: -1.38)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SurfaceDensityCorrectionCoefficientB2",
				OpenSMOKE::SINGLE_DOUBLE,
				"SurfaceDensityCorrectionCoefficientB2 (default: 0.00069)",
				false));

			// ----------------------------------------------------------------------------------------------------------- //
			// Reaction 1: Soot-H + OH = Soot* + H2O
			// ----------------------------------------------------------------------------------------------------------- //
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A1f",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 1 (forward): Soot-H + OH = Soot* + H2O (default: 6.72e1 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A1b",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 1 (backward): Soot-H + OH = Soot* + H2O (default: 6.44e-1 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n1f",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 1 (forward): Soot-H + OH = Soot* + H2O (default: 3.33)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n1b",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 1 (backward): Soot-H + OH = Soot* + H2O (default: 3.79)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E1f",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 1 (forward): Soot-H + OH = Soot* + H2O (default: 6.09 kJ/mol)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E1b",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 1 (backward): Soot-H + OH = Soot* + H2O (default: 27.96 kJ/mol)",
				false));

			// ----------------------------------------------------------------------------------------------------------- //
			// Reaction 2: Soot-H + H = Soot* + H2
			// ----------------------------------------------------------------------------------------------------------- //
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A2f",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 2 (forward): Soot-H + H = Soot* + H2 (default: 1e8 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A2b",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 2 (backward): Soot-H + H = Soot* + H2 (default: 8.68e4 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n2f",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 2 (forward): Soot-H + H = Soot* + H2 (default: 1.80)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n2b",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 2 (backward): Soot-H + H = Soot* + H2 (default: 2.36)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E2f",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 2 (forward): Soot-H + H = Soot* + H2 (default: 68.42 kJ/mol)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E2b",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 2 (backward): Soot-H + H = Soot* + H2 (default: 25.46 kJ/mol)",
				false));

			// ----------------------------------------------------------------------------------------------------------- //
			// Reaction 3: Soot + H = Soot* + H
			// ----------------------------------------------------------------------------------------------------------- //
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A3f",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 3 (forward): Soot + H = Soot* + H (default: 1.13e16 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A3b",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 3 (backward): Soot + H = Soot* + H (default: 4.17e13 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n3f",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 3 (forward): Soot + H = Soot* + H (default: -0.06)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n3b",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 3 (backward): Soot + H = Soot* + H (default: 0.15)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E3f",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 3 (forward): Soot + H = Soot* + H (default: 476.05 kJ/mol)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E3b",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 3 (backward): Soot + H = Soot* + H (default: 0 kJ/mol)",
				false));

			// ----------------------------------------------------------------------------------------------------------- //
			// Reaction 4: Soot* + C2H2 => Soot-H
			// ----------------------------------------------------------------------------------------------------------- //
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A4",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 4: Soot* + C2H2 => Soot-H (default: 2.52e9 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n4",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 4: Soot* + C2H2 => Soot-H (default: 1.10)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E4",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 4: Soot* + C2H2 => Soot-H (default: 17.13 kJ/mol)",
				false));

			// ----------------------------------------------------------------------------------------------------------- //
			// Reaction 5: Soot* + O2 => Soot-H + 2CO
			// ----------------------------------------------------------------------------------------------------------- //
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@A5",
				OpenSMOKE::SINGLE_MEASURE,
				"Frequency factor reaction 5: Soot* + O2 => Soot-H + 2CO (default: 2.20e12 cm3,mol,s)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@n5",
				OpenSMOKE::SINGLE_DOUBLE,
				"Temperature exponent reaction 5: Soot* + O2 => Soot-H + 2CO (default: 0)",
				false));
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@E5",
				OpenSMOKE::SINGLE_MEASURE,
				"Activation energy reaction 5: Soot* + O2 => Soot-H + 2CO (default: 31.38 kJ/mol)",
				false));

			// ----------------------------------------------------------------------------------------------------------- //
			// Reaction 6: Soot-H + OH => Soot-H + CO
			// ----------------------------------------------------------------------------------------------------------- //
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Efficiency6",
				OpenSMOKE::SINGLE_DOUBLE,
				"Efficiency reaction 6: Soot-H + OH => Soot-H + CO (default: 0.13)",
				false));
		}
	};
}

#endif /* OpenSMOKE_Grammar_HMOM_H */

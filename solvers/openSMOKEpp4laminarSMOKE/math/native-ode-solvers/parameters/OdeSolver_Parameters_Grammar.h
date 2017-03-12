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

#ifndef OdeSolver_Parameters_Grammar_H
#define	OdeSolver_Parameters_Grammar_H

#include <string>
#include <vector>
#include <iostream>
#include <iomanip>
#include <math/OpenSMOKEFunctions.h>
#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

namespace OdeSMOKE
{
	class OdeSolver_Parameters_Grammar : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
	{
	protected:

		virtual void DefineRules()
		{
			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@OdeSolver",
				OpenSMOKE::SINGLE_STRING,
				"ODE Solver: OpenSMOKE++ (default) | BzzOde | Cvode | DASPK",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Jacobian",
				OpenSMOKE::SINGLE_STRING,
				"Jacobian pattern: Band (default) | TridiagonalBlock | Sparse",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@SparseSolver",
				OpenSMOKE::SINGLE_STRING,
				"Jacobian pattern: EigenSparseLU (default) | EigenBiCGSTAB | EigenGMRES | EigenDGMRES | Pardiso | SuperLUSerial | UMFPack",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Preconditioner",
				OpenSMOKE::SINGLE_STRING,
				"Preconditioner to be used by iterative solvers: ILUT (default) | diagonal",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@RelativeTolerance",
				OpenSMOKE::SINGLE_DOUBLE,
				"Releative tolerance (default 1e-6)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@AbsoluteTolerance",
				OpenSMOKE::SINGLE_DOUBLE,
				"Releative tolerance (default 1e-10)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumNumberOfSteps",
				OpenSMOKE::SINGLE_INT,
				"Maximum number of steps (default 500000)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumStep",
				OpenSMOKE::SINGLE_DOUBLE,
				"Maximum step (default: automatically chosen by the solver)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumStep",
				OpenSMOKE::SINGLE_DOUBLE,
				"Minimum step (default: automatically chosen by the solver)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@InitialStep",
				OpenSMOKE::SINGLE_DOUBLE,
				"Initial step (default: automatically chosen by the solver)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumOrder",
				OpenSMOKE::SINGLE_INT,
				"Maximum order to be used during the ODE integration",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MeanResidualThreshold",
				OpenSMOKE::SINGLE_DOUBLE,
				"Stop the integration when the mean residual is below this threshold (default: 0.)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@VerbosityLevel",
				OpenSMOKE::SINGLE_INT,
				"Verbosity Level (default 1)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MinimumConstraints",
				OpenSMOKE::SINGLE_BOOL,
				"Constraints on minimum values (available only for OpenSMOKE++ and BzzOde solvers, default true)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@MaximumConstraints",
				OpenSMOKE::SINGLE_BOOL,
				"Constraints on maximum values (available only for OpenSMOKE++ and BzzOde solvers, default false)",
				false));

			AddKeyWord(OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@NonNegativeVariables",
				OpenSMOKE::SINGLE_BOOL,
				"Constraints on minimum values (available only for CVODE and DASPK solvers, default false)",
				false));
		}
	};
}

#endif	/* OdeSolver_Parameters_Grammar_H */


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

#ifndef OdeSolverUtilities_H
#define OdeSolverUtilities_H

namespace OdeSMOKE
{

	enum OdeHStatus
	{
		H_STATUS_DECREASED,
		H_STATUS_CONST,
		H_STATUS_INCREASED
	};

	enum OdeOrderStatus
	{
		ORDER_STATUS_DECREASED,
		ORDER_STATUS_CONST,
		ORDER_STATUS_INCREASED
	};

	enum JacobianType
	{
		JACOBIAN_TYPE_CONST,
		JACOBIAN_TYPE_USERDEFINED,
		JACOBIAN_TYPE_NUMERICAL
	};

	enum JacobianStatus
	{
		JACOBIAN_STATUS_HAS_TO_BE_CHANGED,
		JACOBIAN_STATUS_MODIFIED,
		JACOBIAN_STATUS_OK
	};

	enum FactorizationStatus
	{
		MATRIX_HAS_TO_BE_FACTORIZED,
		MATRIX_FACTORIZED,
	};

	enum OdeConvergence
	{
		CONVERGENCE_STATUS_FAILURE,
		CONVERGENCE_STATUS_OK
	};

	enum OdeStatus
	{
		ODE_STATUS_TO_BE_INITIALIZED = 1,
		ODE_STATUS_CONTINUATION = 2,
		ODE_STATUS_STOP_INTEGRATION_FOR_SMALL_YPRIME_NORM1 = 3,

		ODE_STATUS_MAX_NUMBER_OF_STEPS_REACHED = -1,
		ODE_STATUS_TOO_STRICT_TOLERANCES = -2,
		ODE_STATUS_ILLEGAL_CONTINUATION_REQUEST = -3,
		ODE_STATUS_MAX_NUMBER_ERRORTEST_FAILURES = -4,
		ODE_STATUS_MAX_NUMBER_CONVERGENCETEST_FAILURES = -5,
		ODE_STATUS_TOO_SMALL_STEP_SIZE = -6,
		ODE_STATUS_ILLEGAL_MAX_INDEPENDENT_VARIABLE = -7,
		ODE_STATUS_ILLEGAL_CONSTRAINTS = -9,
		ODE_STATUS_EXCEPTION_HANDLING_STOP = -10,
		ODE_STATUS_YOU_CANNOT_OVERSHOOT_TCRITIC = -12
	};

	unsigned int Factorial(unsigned int n);
}

#include "OdeSolverUtilities.hpp"

#endif // ODESolverUtilities_H
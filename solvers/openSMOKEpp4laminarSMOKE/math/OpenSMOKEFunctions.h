/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef OpenSMOKE_OpenSMOKEFunctions_Hpp
#define OpenSMOKE_OpenSMOKEFunctions_Hpp

#include "OpenSMOKEStdInclude.h"

namespace OpenSMOKE
{
	/**
	* Utility to print fatal error messages
	*/
	void ErrorMessage(const std::string functionName, const std::string errorMessage);

	/**
	* Utility to print fatal error messages
	*/
	int FatalErrorMessage(const std::string errorMessage);

	/**
	* Returns the MachEps (float)
	*/
	float MachEpsFloat();
	
	/**
	* Returns the MachEps (double)
	*/
	double MachEps();

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	float SqrtSumSqr(const int n, float *x);

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	double SqrtSumSqr(const int n, double *x);

	/**
	* Returns the n-th power of x
	* Example: z=PowInt(x,n) means \f$ z=x^n \f$
	*/
	double PowInt(const double x, const int n);

	/**
	* Returns the cubic root of x
	* Example: z=CubicRoot(x) means \f$ z=x^(1/3) \f$
	*/
	double CubicRoot(const double x);
        
    /**
	* Returns the elapsed time in seconds since the process started
	*/
    double OpenSMOKEClock(void);

    /**
	* Returns the elapsed time in seconds since the OpenSMOKEDiffClock() started
	*/
	double OpenSMOKEGetCpuTime(void);

    /**
	* Returns the number of cpu clocks
	*/
	#if defined(_WIN32) || defined(_WIN64) 
		unsigned __int64 OpenSMOKEGetCpuClocks(void);
	#else
		unsigned long int OpenSMOKEGetCpuClocks(void);
	#endif

	/**
	* Returns the cpu frequency
	*/
    double OpenSMOKEGetCpuFrequency(void);

    /**
	* Returns the max cpu frequency
	*/
    double OpenSMOKEGetMaxCpuFrequency(void);

    /**
	* Returns the max cpu clock frequency
	*/
    double OpenSMOKEGetCpuClocksFrequency(void);
}

#include "OpenSMOKEFunctions.hpp"

#endif	// OpenSMOKE_OpenSMOKEFunctions_Hpp


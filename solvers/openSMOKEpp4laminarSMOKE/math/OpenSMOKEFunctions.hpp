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

#include <stdlib.h>
#include <time.h>

#if OPENSMOKE_USE_MKL == 1
	#include <mkl.h>
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "lapacke.h"
#endif

namespace OpenSMOKE
{
	static const float  OPENSMOKE_MACH_EPS_FLOAT = MachEpsFloat();
	static const double OPENSMOKE_MACH_EPS_DOUBLE = MachEps();

	void ErrorMessage(const std::string functionName, const std::string errorMessage)
	{
		std::cout << "Function:     " << functionName << std::endl;
		std::cout << "Fatal error:  " << errorMessage << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	int FatalErrorMessage(const std::string errorMessage)
	{
		std::cout << "Fatal error:  " << errorMessage << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
		return OPENSMOKE_FATAL_ERROR_EXIT;
	}

	float MachEpsFloat()
	{
		float macheps = 1.F;
		float eps = 2.F;
	
		while(eps != 1.F)
		{
			macheps /= 2.F;
			eps = float(1.F + macheps); // in simple rounded form
		}
	
		if(macheps < 1.e-7)
			return 1.192092896e-07F;
		else
			return macheps*2.F;
	}

	double MachEps()
	{
		double macheps = 1.;
		double eps = 2.;
		while(eps != 1.)
		{
			macheps /= 2.;
			eps = 1. + macheps; // in double rounded form
		}
		return macheps*2.;
	}

	float SqrtSumSqr(const int n, float *x)
	{	
		double norm = 0.;
		for(int i=0;i<n;i++)
		{
			double aux = x[i];
			norm += aux*aux;
		}
		norm = std::sqrt(norm);
		if(norm > OPENSMOKE_BIG_FLOAT)
			norm = OPENSMOKE_BIG_FLOAT;
		return float(norm);
	}

	double SqrtSumSqr(const int n, double *x)
	{
		double aux;
		double xmax = 0.;
		double xmin = OPENSMOKE_BIG_DOUBLE;
		for(int j=0;j<n;j++)
		{
			aux = std::fabs(x[j]);
			if(xmax<aux)	xmax = aux;
			if(xmin>aux)	xmin = aux;
		}
		
		if(xmax == 0.)
			return xmax;
		
		if (xmin == 0.)
			xmin = OPENSMOKE_TINY_DOUBLE;
		
		aux = std::sqrt(OPENSMOKE_BIG_DOUBLE/((double)n));
		
		// to avoid the problems of
		#if OPENSMOKE_LONG_DOUBLE == OPENSMOKE_DOUBLE
			if	( xmax < aux &&												// overflow
				  xmax > OPENSMOKE_TINY_DOUBLE/OPENSMOKE_MACH_EPS_FLOAT)	 // small numbers
		#else
			long double longaux = (long double)xmax/(long double)xmin;
			if(	xmax < aux &&									// overflow
				xmax > OPENSMOKE_TINY_DOUBLE/OPENSMOKE_MACH_EPS_FLOAT &&	// small numbers
				longaux < 1./OPENSMOKE_MACH_EPS_FLOAT)				// sort
		#endif
			{
				double norm = 0.;	// without problems: double
				for(int i = 0;i < n;i++)
				{
					aux = x[i];
					norm += aux*aux;
				}
				return std::sqrt(norm);
			}
			else	// if there are problems it works in long double
			{
				#if OPENSMOKE_LONG_DOUBLE == OPENSMOKE_DOUBLE
		
					double *y = new double[n];
					
					if(!y)
					{
						std::cout << "Not enough memory in double SqrtSumSqr(int n,double *x) function" << std::endl;
						std::cout << "Press enter to exit..." << std::endl;
						getchar();
						exit(OPENSMOKE_FATAL_ERROR_EXIT);
					}
					
					for(int i=0;i<n;i++)
						y[i] = x[i]/xmax;
					double norm = 0.;
					for(int i=0;i<n;i++)
						norm += y[i]*y[i];
					delete y;

					return xmax*std::sqrt(norm);
		
				#else
					long double norm = 0.;
					for(int i = 0;i < n;i++)
					{
						longaux = x[i];
						norm += longaux*longaux;
					}
					if(norm < OPENSMOKE_BIG_DOUBLE && norm > OPENSMOKE_TINY_DOUBLE)
						return std::sqrt(double(norm));
					longaux = (long double)xmax*(long double)n;
					norm /= longaux; // avoids overflow
					norm /= longaux;
					norm = longaux*std::sqrt(double(norm)); // renormalises
					
					// avoids overflow
					if(norm > OPENSMOKE_BIG_DOUBLE) norm = OPENSMOKE_BIG_DOUBLE;
						return double(norm);
				#endif
		}
	}

	double PowInt(const double x, const int n)
	{
		if(n > 0)
		{
			double y;
			double x2,x4,x5,x6,x8,x9,x10;
			switch(n)
			{
				case 0:
					return 1.;
				case 1:
					return x;
				case 2:
					return x*x;
				case 3:
					return x*x*x;
				case 4:
					x2 = x*x;
					return x2*x2;
				case 5:
					x2 = x*x;
					return x2*x2*x;
				case 6:
					x2 = x*x;
					return x2*x2*x2;
				case 7:
					x2 = x*x;
					return x2*x2*x2*x;
				case 8:
					x2 = x*x;
					x4 = x2*x2;
					return x4*x4;
				case 9:
					x2 = x*x;
					x4 = x2*x2;
					return x4*x4*x;
				case 10:
					x2 = x*x;
					x4 = x2*x2;
					x5 = x4*x;
					return x5*x5;
				case 11:
					x2 = x*x;
					x4 = x2*x2;
					x5 = x4*x;
					return x*x5*x5;
				case 12:
					x2 = x*x;
					x4 = x2*x2;
					x6 = x4*x2;
					return x6*x6;
				case 13:
					x2 = x*x;
					x4 = x2*x2;
					x6 = x4*x2;
					return x*x6*x6;
				case 14:
					x2 = x*x;
					x4 = x2*x2;
					x6 = x4*x2;
					return x2*x6*x6;
				case 15:
					x2 = x*x;
					x4 = x2*x2;
					x6 = x4*x2;
					return x*x2*x6*x6;
				case 16:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					return x8*x8;
				case 17:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					return x*x8*x8;
				case 18:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					x9 = x*x8;
					return x9*x9;
				case 19:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					x9 = x*x8;
					return x*x9*x9;
				case 20:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					x10 = x2*x8;
					return x10*x10;
				case 21:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					x10 = x2*x8;
					return x*x10*x10;
				case 22:
					x2 = x*x;
					x4 = x2*x2;
					x8 = x4*x4;
					x10 = x2*x8;
					return x2*x10*x10;

				default:

					int n_ = n;
					double x_ = x;

					y = 1.;
					while(n_ != 0)
					{
						if(n_ & 1)
							y *= x_;
						x_ *= x_;
						n_ >>= 1;
					}
			}
		
			return y;
		}
	
		if(x == 0.)
		{
			ErrorMessage("PowInt(double x, int n)", "x^n with x=0. and n<=0");
			return 0.;
		}
		else if(n < 0)
			return 1./PowInt(x,-n);
		else
			return 1.;
	}

	double CubicRoot(const double x)
	{
		if(x<0.)
			return -std::pow(-x, 1./3.);
		else
			return  std::pow(x, 1./3.);
	}
        
    double OpenSMOKEClock(void)
	{
		return (double)(std::clock())/CLOCKS_PER_SEC;
	}

    double OpenSMOKEGetCpuTime(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return dsecnd();
		#else
            return OpenSMOKEClock();
		#endif
    }

	#if defined(_WIN32) || defined(_WIN64) 
		unsigned __int64 OpenSMOKEGetCpuClocks(void)
	#else
		unsigned long int OpenSMOKEGetCpuClocks(void)
	#endif
    {
		#if OPENSMOKE_USE_MKL == 1
			unsigned MKL_INT64 clocks;
			mkl_get_cpu_clocks( &clocks );
			return clocks;
		#else
            return std::clock();
		#endif	
	}

    double OpenSMOKEGetCpuFrequency(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return mkl_get_cpu_frequency();
		#else
            return 0;
		#endif
	}

    double OpenSMOKEGetMaxCpuFrequency(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return 0;
			//return mkl_get_max_cpu_frequency();
		#else
            return 0;
		#endif
	}

    double OpenSMOKEGetCpuClocksFrequency(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return 0;
			//return mkl_get_clocks_frequency();
		#else
            return 0;
		#endif
	}
}

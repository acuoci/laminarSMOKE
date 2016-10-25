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

// LIS libraries: must be included before any other inclusion
#if (OPENSMOKE_USE_LIS == 1)
	#ifdef HAVE_CONFIG_H
		#include "lis_config.h"
	#else
		#ifdef HAVE_CONFIG_WIN_H
			#include "lis_config_win.h"
		#endif
	#endif
	#include "lis.h"
#endif

// OPENSMOKE_USE_MKL and OPENSMOKE_USE_OPENBLAS cannot be used together
#if (OPENSMOKE_USE_MKL == 1 && OPENSMOKE_USE_OPENBLAS == 1)
	BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "Only 1 between OPENSMOKE_USE_MKL and OPENSMOKE_USE_OPENBLAS can be available");
#endif

// If OPENSMOKE_USE_MKL is available, Eigen can exploit it
#if OPENSMOKE_USE_MKL == 1
	#define EIGEN_USE_MKL_ALL
#endif

// Bug with Microsoft Visual Studio and g++
#if OPENSMOKE_USE_OPENBLAS == 1
	#include <complex>
	#define lapack_complex_float std::complex<float>
	#define lapack_complex_double std::complex<double>
#endif

// Sundials package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_SUNDIALS == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
		BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_SUNDIALS is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// DVODE package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_DVODE == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
		BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_DVODE is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// ODEPACK package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_ODEPACK == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
		BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_ODEPACK is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// DASPK package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_DASPK == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
		BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_DASPK is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// RADAU package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_RADAU == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
		BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_RADAU is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif


// MEBDF package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_MEBDF == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
		BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_MEBDF is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// SuperLUSerial package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_SUPERLU_SERIAL == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
			BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_SUPERLU_SERIAL is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// SuperLUSerial package requires OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS
#if (OPENSMOKE_USE_UMFPACK == 1)
	#if (OPENSMOKE_USE_MKL == 0 && OPENSMOKE_USE_OPENBLAS == 0)
				BOOST_STATIC_ASSERT_MSG(sizeof(T) == 0, "If OPENSMOKE_USE_UMFPACK is enabled, OPENSMOKE_USE_MKL or OPENSMOKE_USE_OPENBLAS is required");
	#endif
#endif

// ARMADILLO C++ Library
#define OPENSMOKE_USE_ARMADILLO	0

// PLASMA C++ Library
#define OPENSMOKE_USE_PLASMA	0

// FLAME C Library
#define OPENSMOKE_USE_FLAME	0


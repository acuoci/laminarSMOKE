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

// Intel MKL Libraries (fast numerical functions)
#ifndef OPENSMOKE_USE_MKL
#define OPENSMOKE_USE_MKL         0
#endif

// CLAPACK (e.g. OpenBLAS)
#ifndef OPENSMOKE_USE_OPENBLAS
#define OPENSMOKE_USE_OPENBLAS    0
#endif

// Sundials C Library
#ifndef OPENSMOKE_USE_SUNDIALS
#define OPENSMOKE_USE_SUNDIALS    0
#endif

// DVODE Fortran Library
#ifndef OPENSMOKE_USE_DVODE
#define OPENSMOKE_USE_DVODE       0
#endif

// ODEPACK Fortran Library
#ifndef OPENSMOKE_USE_ODEPACK
#define OPENSMOKE_USE_ODEPACK     0
#endif

// DASPK Fortran Library
#ifndef OPENSMOKE_USE_DASPK
#define OPENSMOKE_USE_DASPK       0
#endif

// RADAU Fortran Library
#ifndef OPENSMOKE_USE_RADAU
#define OPENSMOKE_USE_RADAU       0
#endif

// MEBDF Fortran Library
#ifndef OPENSMOKE_USE_MEBDF
#define OPENSMOKE_USE_MEBDF       0
#endif

// SuperLU Serial Library
#ifndef OPENSMOKE_USE_SUPERLU_SERIAL
#define OPENSMOKE_USE_SUPERLU_SERIAL       0
#endif

// UMFPACK Library
#ifndef OPENSMOKE_USE_UMFPACK
#define OPENSMOKE_USE_UMFPACK       0
#endif

// LIS Libraries
#ifndef OPENSMOKE_USE_LIS
#define OPENSMOKE_USE_LIS         0
#endif

// BzzMath C++ Library
#ifndef OPENSMOKE_USE_BZZMATH
#define OPENSMOKE_USE_BZZMATH     0
#endif

// PLASMA Library (To be added in the future)
#ifndef OPENSMOKE_USE_PLASMA
#define OPENSMOKE_USE_PLASMA     0
#endif

// FLAME Library (To be added in the future)
#ifndef OPENSMOKE_USE_FLAME
#define OPENSMOKE_USE_FLAME     0
#endif

// Check the definitions reported above
#include "math/OpenSMOKE_CheckDefinitions.h"

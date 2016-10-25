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

#ifndef OpenSMOKE_OpenSMOKELoad_Hpp
#define OpenSMOKE_OpenSMOKELoad_Hpp

#include "OpenSMOKEStdInclude.h"
#include "OpenSMOKEBaseClass.h"
#include "OpenSMOKEUtilities.h"
#include "OpenSMOKEFunctions.h"

namespace OpenSMOKE
{
	/**
	*  OpenSMOKE Load Class
	*/

	class OpenSMOKELoad : public OpenSMOKEBaseClass
	{
	
		public:

			// Default constructor
			//inline OpenSMOKELoad(void);

			// Copy-initializer
			inline OpenSMOKELoad(OpenSMOKELoad &rval);
	
			// File constructor;
			inline OpenSMOKELoad(const std::string file, const OpenSMOKE_File_Format fileFormat);

			// Destructor
			inline ~OpenSMOKELoad(void);

			// Operator ()
			inline void operator()(const std::string file, const OpenSMOKE_File_Format fileFormat);
			

			// Utilities
			inline std::string FileName() const;
			inline void End(void);

			// Load functions
			template<typename T_>
			inline OpenSMOKELoad &operator >> (T_ &x);

			template<typename T_, typename IndexPolicy_>
			inline OpenSMOKELoad &operator >> (OpenSMOKEVector<T_, IndexPolicy_> &v);

			template<typename T_, typename IndexPolicy_>
			inline OpenSMOKELoad &operator >> (OpenSMOKEMatrix<T_, IndexPolicy_> &A);
	
			//template<typename T_, typename IndexPolicy_>
			//void readVectorFromBinaryFile(std::ifstream &fInput, int N, OpenSMOKEVector<T_, IndexPolicy_> &v);
std::ifstream fileLoad_;
		private:
	
			OpenSMOKE_File_Format fileFormat_;
			std::string fileName_;
			
	};
}

#include "OpenSMOKELoad.hpp"

#endif	// OpenSMOKE_OpenSMOKEVector_Hpp


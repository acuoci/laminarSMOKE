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


//#include "math/OpenSMOKELoad.h"
#include "math/OpenSMOKEVector.h"
#include "math/OpenSMOKEMatrix.h"
#include <typeinfo>
#include <stdio.h>
#include <string.h>

namespace OpenSMOKE
{

	/*inline OpenSMOKELoad::OpenSMOKELoad(void)
	{
	}*/

	inline OpenSMOKELoad::OpenSMOKELoad(OpenSMOKELoad &rval)
	{	
		ErrorMessage("OpenSMOKELoad::OpenSMOKELoad(OpenSMOKELoad &rval): No copy constructor available");
	}

	inline OpenSMOKELoad::~OpenSMOKELoad(void)
	{
		fileLoad_.close();
	}

	inline OpenSMOKELoad::OpenSMOKELoad(const std::string file, const OpenSMOKE_File_Format fileFormat)
	{
		fileName_ = file;
		fileFormat_ = fileFormat;

		if (fileFormat_ == OPENSMOKE_FORMATTED_FILE)
			fileLoad_.open(file.c_str(), std::ios::in);
		else if (fileFormat_ == OPENSMOKE_BINARY_FILE)
			fileLoad_.open(file.c_str(), std::ios::in | std::ios::binary);

		if (!fileLoad_.is_open()) 
			ErrorMessage("OpenSMOKELoad::OpenSMOKELoad(const std::string file, const OpenSMOKE_File_Format fileFormat): File already open");
	}

	inline void OpenSMOKELoad::operator()(const std::string file, const OpenSMOKE_File_Format fileFormat)
	{
		fileName_ = file;
		fileFormat_ = fileFormat;

		if (fileFormat_ == OPENSMOKE_FORMATTED_FILE)
			fileLoad_.open(file.c_str(), std::ios::in);
		else if (fileFormat_ == OPENSMOKE_BINARY_FILE)
			fileLoad_.open(file.c_str(), std::ios::in | std::ios::binary);

		if (!fileLoad_.is_open()) 
			ErrorMessage("void OpenSMOKELoad::operator()(const std::string file, const OpenSMOKE_File_Format fileFormat): File already open");
	}

	template<typename T>
	inline OpenSMOKELoad &OpenSMOKELoad::operator >> (T &x)
	{
		if(fileFormat_ == OPENSMOKE_BINARY_FILE)
		{
			if(!fileLoad_.read(reinterpret_cast < char * > (&x),sizeof(T)))
				ErrorMessage("OpenSMOKELoad &OpenSMOKELoad::operator >> (T &x): Error in reading");
		}
		else 
		{
			fileLoad_ >> x;
		}

		return *this;
	}

	template<typename T, typename IndexPolicy>
	inline OpenSMOKELoad &OpenSMOKELoad::operator >> (OpenSMOKEVector<T, IndexPolicy> &v)
	{
		v.Load(fileLoad_, fileFormat_);
		return *this;
	}

	template<typename T, typename IndexPolicy>
	inline OpenSMOKELoad &OpenSMOKELoad::operator >> (OpenSMOKEMatrix<T, IndexPolicy> &A)
	{
		A.Load(fileLoad_, fileFormat_);
		return *this;
	}


	inline std::string OpenSMOKELoad::FileName() const
	{
		return fileName_;
	}

	inline void OpenSMOKELoad::End(void)
	{
		fileLoad_.close();
	}
}

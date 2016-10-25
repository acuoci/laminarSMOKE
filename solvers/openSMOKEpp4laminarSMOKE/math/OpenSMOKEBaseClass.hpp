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

#include <typeinfo>
#include <stdio.h>
#include <stdlib.h>
#include "math/OpenSMOKEBaseClass.h"

namespace OpenSMOKE
{
	
	void OpenSMOKEBaseClass::ErrorMessage(const std::string message) const
	{
		std::cout << "Class:   " << typeid(OpenSMOKEBaseClass).name() << std::endl;
		std::cout << "Index:   " << whoAmI_ << "/" << count_ << std::endl;
		std::cout << "Message: " << message << std::endl;
		std::cout << "Press enter to continue..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	unsigned int OpenSMOKEBaseClass::count_ = 0;
	
	unsigned int OpenSMOKEBaseClass::countInScope_ = 0;

}

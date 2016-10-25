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

#ifndef OpenSMOKE_OpenSMOKEBaseClass_Hpp
#define OpenSMOKE_OpenSMOKEBaseClass_Hpp

#include "OpenSMOKEStdInclude.h"

namespace OpenSMOKE
{
	//!  A base class for the OpenSMOKE library
	/*!
		 This is the OpenSMOKE Base Class, which is common to every OpenSMOKE class
	*/

	class OpenSMOKEBaseClass
	{
		protected:

			static unsigned int count_;				/**< total number of available object of this type */ 
			static unsigned int countInScope_;		/**< total number of available object of this type ???  */  

			unsigned int whoAmI_;					/**< index of current object  */  

			/**
			* Return a fatal error message
			*/
			inline void ErrorMessage(const std::string message) const;

			
			/**
			* Initialize counters
			*/
			inline void SetCounters()
			{
				count_++;
				countInScope_++;	
				whoAmI_ = count_;
			}

		public:
	
			/**
			* Returns the index of current object type
			*/
			inline unsigned int WhoAmI() const				{	return whoAmI_;	}

			/**
			* Returns the total number of currently available objects
			*/
			inline static unsigned int ObjectCount()		{	return count_;	}

			/**
			* Returns the total number of currently available objects (???)
			*/
			inline static unsigned int ObjectCountInScope()	{	return countInScope_;	}

	};
}

#include "OpenSMOKEBaseClass.hpp"

#endif	// OpenSMOKE_OpenSMOKEBaseClass_Hpp

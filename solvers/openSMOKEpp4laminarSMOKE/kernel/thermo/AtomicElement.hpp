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

#include <boost/algorithm/string.hpp>

namespace OpenSMOKE
{
	AtomicElement::AtomicElement()
	{
		description_ = "[not defined]";
		name_ = "[not defined]";
		mw_ = 0;
	}

	AtomicElement::AtomicElement(const std::string description__, const std::string name__, const double mw__) : 
		description_(description__), 
		name_(name__),
		mw_(mw__)
	{
		boost::to_upper(name_);
	}

	AtomicElement::AtomicElement(const std::string name__, const double mw__) :
		name_(name__),
		mw_(mw__)
	{
		name_ = "[not defined]";
	}

	AtomicElement::AtomicElement( const AtomicElement& other ) :
		description_(other.description_), 
		name_(other.name_),
		mw_(other.mw_)
	{
		boost::to_upper(name_);
	}

	AtomicElement& AtomicElement::operator=( const AtomicElement& rhs )
	{
		description_ = rhs.description_;
		name_ = rhs.name_;
		mw_ = rhs.mw_;
		boost::to_upper(name_);

		return *this;
	}
	
}
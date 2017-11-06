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

#ifndef OpenSMOKE_PreProcessorKinetics_H
#define	OpenSMOKE_PreProcessorKinetics_H

namespace OpenSMOKE
{

	//!  Abstract class to preprocess a kinetic mechanism 
	/*!
		 This class provides a common interface to different user-defined
		 kinetic pre-processors.
	*/

	template<typename KineticsPolicy>
	class PreProcessorKinetics : public KineticsPolicy 
	{
	public:

		/**
		* Default constructor
		*/
		PreProcessorKinetics(std::ostream& flog);

		/**
		* Copy constructor
		*/
		PreProcessorKinetics(const PreProcessorKinetics& orig);
    
		/**
		* Default destructor
		*/
		virtual ~PreProcessorKinetics();

		/**
		* Reads the kinetic scheme from ASCII file
		*/
		bool ReadKineticsFromASCIIFile(AtomicCompositionTable& compositionTable);

		/**
		* Reads the kinetic scheme from ASCII file
		*/
		template<typename PreProcessor>
		bool ReadKineticsFromASCIIFile(AtomicCompositionTable& compositionTable, const PreProcessor& preprocessor_kinetics);

    
	private:

		std::ostream& flog_;	//!< log file where useful information is written during the pre-processing of a kinetic mechanism
  
	};

}

#include "PreProcessorKinetics.hpp"

#endif	/* OpenSMOKE_PreProcessorKinetics_H */


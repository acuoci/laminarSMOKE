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

namespace OpenSMOKE
{

	template<typename KineticsPolicy>
	PreProcessorKinetics<KineticsPolicy>::PreProcessorKinetics(std::ostream& flog) :
		flog_(flog)
	{
	}

	template<typename KineticsPolicy>
	PreProcessorKinetics<KineticsPolicy>::PreProcessorKinetics(const PreProcessorKinetics& orig) {
	}

	template<typename KineticsPolicy>
	PreProcessorKinetics<KineticsPolicy>::~PreProcessorKinetics() {
	}

	template<typename KineticsPolicy>
	bool PreProcessorKinetics<KineticsPolicy>::ReadKineticsFromASCIIFile(AtomicCompositionTable& compositionTable)
	{
		return this->KineticsFromASCIIFile(compositionTable, flog_);
	}

	template<typename KineticsPolicy>
	template<typename PreProcessor>
	bool PreProcessorKinetics<KineticsPolicy>::ReadKineticsFromASCIIFile(AtomicCompositionTable& compositionTable, const PreProcessor& preprocessor_kinetics)
	{
		return this->KineticsFromASCIIFile(compositionTable, preprocessor_kinetics, flog_);
	}
}

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

#ifndef OpenSMOKE_FluxAnalysisMap_H
#define OpenSMOKE_FluxAnalysisMap_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	//!  A class to perform flux analysis
	/*!
			A class to perform flux analysis
	*/

	class FluxAnalysisMap
	{
	public:

		FluxAnalysisMap(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMapXML,
						OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML);

		void SetDestructionAnalysis(const bool flag) { destruction_analysis_ = flag; }

		void SetLogarithmicThickness(const bool flag) { logarithmic_thickness_ = flag; }

		void SetNormalThickness(const bool flag) { normal_thickness_ = flag; }
	
		void SetNormalTags(const bool flag) { normal_tags_ = flag; }

		void SetMaxDepth(const int max_depth) { max_depth_ = max_depth; }
	
		void SetMinPercentageThreshold(const double min_percentage_threshold) { min_percentage_threshold_ = min_percentage_threshold; }
	
		void SetMaxWidth(const int max_width) { max_width_ = max_width; }
	
		void SetAtom(const unsigned int index_atom) { index_atom_ = index_atom; }
	
		void SetReactionRates(unsigned int n, const double* r);

		void GloballyAnalyze(const std::vector<unsigned int>& important_indices, const int current_depth );

		void CalculateThickness();

		void Plot(const std::string file_name);

		void OpenGraphFile(const std::string file_name);

		void CloseGraphFile();

		void WriteStoichiometricMatrix( const std::string file_name );


	private:

		OpenSMOKE::ThermodynamicsMap_CHEMKIN&  thermodynamicsMapXML_;
		OpenSMOKE::KineticsMap_CHEMKIN& kineticsMapXML_;

		unsigned int NC;
		unsigned int NR;
		unsigned int index_atom_;
		std::vector<double> r__;
		int max_width_;
		int max_depth_;
		double min_percentage_threshold_;
		std::vector<unsigned int> list_of_analyzed_species_;
		std::ofstream fOut;
		bool normal_thickness_;
		bool normal_tags_;
		bool destruction_analysis_;
		bool logarithmic_thickness_;

		void AnalyzeNetFluxes(	const unsigned int index_j, std::vector<unsigned int>& important_indices,
								std::vector<double>& important_normal_fluxes, std::vector<double>& important_fluxes);
	
		void AddSpeciesToGraphFile(	const unsigned int index_j, std::vector<unsigned int>& local_indices,
									std::vector<double>& local_thickness,
									std::vector<double>& local_normal_fluxes,
									std::vector<double>& local_fluxes);

		std::vector<std::vector<unsigned int> > global_important_indices_;
		std::vector<std::vector<double>	>		global_important_normal_fluxes_;
		std::vector<std::vector<double>	>		global_important_fluxes_;
		std::vector<std::vector<double>	>		global_relative_thickness_;
	};
}

#include "FluxAnalysisMap.hpp"

#endif // OpenSMOKE_FluxAnalysisMap_H

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

#ifndef OpenSMOKE_JacobianSparsityPatternMap_H
#define OpenSMOKE_JacobianSparsityPatternMap_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include <Eigen/Sparse>

namespace OpenSMOKE
{
	//!  A class containing the data about the stoichiometry and the reaction orders
	/*!
		 This class provides the tools to manage the stoichiometry and the reaction orders of all the 
		 reactions included in the kinetic scheme. This class is specifically conceived in order to perform 
		 the calculations as fast as possible. To reach this goal the kinetic data are stored in a
		 format which is not so easy to manage and to understand.
	*/

	template<typename map>
	class JacobianSparsityPatternMap
	{

	public:

		/**
		*@brief Default constructor
		*@param nspecies number of species 
		*@param nreactions number of reactions
		*/
		JacobianSparsityPatternMap(map& kinetics_map);

		/**
		*@brief Default destructor
		*/
		~JacobianSparsityPatternMap();

		void SetEpsilon(const double epsilon);

		void epsilon() const { return epsilon_; }
                
		void RecognizeJacobianSparsityPattern(std::vector<unsigned int>& row, std::vector<unsigned int>& col);
		
		void Jacobian(const OpenSMOKE::OpenSMOKEVectorDouble& omega, const double T, const double P_Pa, Eigen::SparseMatrix<double> &J);

		void Jacobian(const OpenSMOKE::OpenSMOKEVectorDouble& omega, const double T, const double P_Pa, Eigen::VectorXd &Jdiagonal);


		Eigen::SparseMatrix<double>* jacobian_matrix() { return jacobian_matrix_; }
		Eigen::SparseMatrix<double>* drf_over_domega() { return drf_over_domega_; }
		Eigen::SparseMatrix<double>* drb_over_domega() { return drb_over_domega_; }

	private:

		map& kinetics_map_;	//!< reference to the kinetic map 

		Eigen::SparseMatrix<double>* drf_over_domega_;
		Eigen::SparseMatrix<double>* drb_over_domega_;
		Eigen::SparseMatrix<double>* dthirdbody_over_domega_;
		Eigen::SparseMatrix<double>* dfalloff_over_domega_;
		Eigen::SparseMatrix<double>* dcabr_over_domega_;
		Eigen::SparseMatrix<double>* jacobian_matrix_;

		Eigen::VectorXi analytical_thirdbody_reactions_;
		Eigen::VectorXi analytical_thirdbody_species_;
		Eigen::VectorXi analytical_falloff_reactions_;

		OpenSMOKE::OpenSMOKEVectorDouble analytical_omegaStar_;
		OpenSMOKE::OpenSMOKEVectorDouble analytical_cStar_;
		OpenSMOKE::OpenSMOKEVectorDouble analytical_xStar_;
		OpenSMOKE::OpenSMOKEVectorDouble analytical_RStar_;
		OpenSMOKE::OpenSMOKEVectorDouble analytical_rf_;
		OpenSMOKE::OpenSMOKEVectorDouble analytical_rb_;

		unsigned int nr;
		unsigned int nc;

		double epsilon_;
		double sum_mass_fractions_;
	};
}

#include "JacobianSparsityPatternMap.hpp"

#endif /* JacobianSparsityPatternMap_H */

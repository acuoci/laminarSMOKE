/*-----------------------------------------------------------------------*\
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
|	License                                                           |
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

#include <Eigen/Sparse>

namespace OpenSMOKE
{
	//!  A class for managing the Jacobian matrix in sparse formulation
	/*!
	This class provides the tools to manage the Jacobian matrix in sparse formulation
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

		/**
		*@brief Returns the sparsity pattern as a couple o vectors containing
		*       the indices (rows and colums) of non-zero elements
		*/
		void RecognizeJacobianSparsityPattern(std::vector<unsigned int>& row, std::vector<unsigned int>& col);

		/**
		*@brief Returns the Jacobian matrix (sparse), given the current operating conditions
		*@param omega current mass fractions
		*@param T current temperature (in K)
		*@param P_Pa current pressure (in Pa)
		*@returns the Jacobian matrix
		*/
		void Jacobian(const double* omega, const double T, const double P_Pa, Eigen::SparseMatrix<double> &J);

		/**
		*@brief Returns only the diagonal elements of the Jacobian matrix (sparse), given the current operating conditions
		*@param omega current mass fractions
		*@param T current temperature (in K)
		*@param P_Pa current pressure (in Pa)
		*@returns the diagonal elements of the Jacobian matrix
		*/
		void Jacobian(const double* omega, const double T, const double P_Pa, Eigen::VectorXd &Jdiagonal);

		/**
		*@brief Returns the Jacobian matrix
		*/
		Eigen::SparseMatrix<double>* jacobian_matrix() { return jacobian_matrix_; }

		/**
		*@brief Returns the matrix of derivatives of forward reaction rates with respect to mass fractions
		*/
		Eigen::SparseMatrix<double>* drf_over_domega() { return drf_over_domega_; }

		/**
		*@brief Returns the matrix of derivatives of backward reaction rates with respect to mass fractions
		*/
		Eigen::SparseMatrix<double>* drb_over_domega() { return drb_over_domega_; }

		/**
		*@brief Sets epsilon, i.e. the minimum mass fraction of species
		*       needed to numerically calculate the Jacobian (default 1e-15)
		*/
		void SetEpsilon(const double epsilon);

		/**
		*@brief Returns epsilon, i.e. the minimum mass fraction of species
		*       needed to numerically calculate the Jacobian (default 1e-15)
		*/
		double epsilon() const { return epsilon_; }

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

		Eigen::VectorXd analytical_omegaStar_;
		Eigen::VectorXd analytical_cStar_;
		Eigen::VectorXd analytical_xStar_;
		Eigen::VectorXd analytical_RStar_;
		Eigen::VectorXd analytical_rf_;
		Eigen::VectorXd analytical_rb_;

		unsigned int nr;
		unsigned int nc;

		double epsilon_;
		double sum_mass_fractions_;
	};
}

#include "JacobianSparsityPatternMap.hpp"

#endif /* JacobianSparsityPatternMap_H */
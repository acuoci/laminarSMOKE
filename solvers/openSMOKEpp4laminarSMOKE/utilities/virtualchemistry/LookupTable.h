/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Alberto Cuoci, Giampaolo Maio, Benoit Fiorina                |
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
|   Copyright(C) 2018 Alberto Cuoci                                       |
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

#include "Eigen/Dense"
#include "math.h"

#ifndef OpenSMOKEpp_LookupTable
#define OpenSMOKEpp_LookupTable

namespace OpenSMOKE
{
	//!  A class to manage lookup tables for virtual chemistry applications
	/*!
	This class provides the tools to manage lookup tables for virtual chemistry applications
	*/

	class LookupTable
	{
	public:

		void Setup(const boost::filesystem::path path_table);
		void Interpolation(const double x);
		double Interpolation(const double x, const unsigned int k);

		inline const Eigen::VectorXd& interpolated() const { return interpolated_; }
		unsigned int nv() const { return nv_; }
		unsigned int np() const { return np_; }

	private:

		unsigned int nv_;
		unsigned int np_;
		double min_x_;
		double max_x_;
		Eigen::VectorXd	x_;
		Eigen::MatrixXd v_;
		Eigen::MatrixXd ratios_;
		std::vector<bool> iRegressions_;
		Eigen::VectorXd	interpolated_;

	};
}

#include "LookupTable.hpp"

#endif /* OpenSMOKEpp_LookupTable */

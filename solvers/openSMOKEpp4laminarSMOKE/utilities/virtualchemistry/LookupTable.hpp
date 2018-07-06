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

namespace OpenSMOKE
{
	void LookupTable::Setup(const boost::filesystem::path path_table)
	{
		std::cout << "Importing virtual chemistry table..." << std::endl;

		// Open the table file
		std::ifstream fInput(path_table.c_str(), std::ios::in);

		// Read the comment line
		std::string comment_line;
		std::getline(fInput, comment_line);

		// Read table size
		fInput >> nv_;
		fInput >> np_;

		// Memory allocation
		v_.resize(np_, nv_);
		x_.resize(np_);
		ratios_.resize(np_ - 1, nv_);
		interpolated_.resize(nv_);
		iRegressions_.resize(nv_);
		for (unsigned int i = 0; i < nv_; i++)
			iRegressions_[i] = false;

		// Fill the table
		for (unsigned int i = 0; i<np_; i++)
		{
			fInput >> x_(i);
			for (unsigned int k = 0; k<nv_; k++)
				fInput >> v_(i, k);
		}

		// Check closure
		std::string dummy;
		fInput >> dummy;
		if (dummy != "END")
		{
			std::cout << "Virtual Chemistry: error in reading the table. Expected END, Found " << dummy << std::endl;
			abort();
		}

		// Close input file
		fInput.close();

		// Min and max values
		min_x_ = x_.minCoeff();
		max_x_ = x_.maxCoeff();

		std::cout << "Virtual chemistry: min(x)=" << min_x_ << " - max(x)=" << max_x_ << std::endl;

		// Precalculation of ratios to be used in interpolation
		for (unsigned int i = 1; i<np_; i++)
			for (unsigned int k = 0; k<nv_; k++)
				ratios_(i - 1, k) = (v_(i, k) - v_(i - 1, k)) / (x_(i) - x_(i - 1));

		std::cout << "Lookup table" << std::endl;
		for (unsigned int i = 0; i < np_; i++)
		{
			std::cout << std::setw(12) << std::left << x_(i);
			for (unsigned int k = 0; k < nv_; k++)
				std::cout << std::setw(12) << std::left << v_(i, k);
			std::cout << std::endl;
		}
		std::cout << std::endl;

		std::cout << "Lookup table (ratios)" << std::endl;
		for (unsigned int i = 0; i < np_ - 1; i++)
		{
			std::cout << std::setw(12) << std::left << x_(i);
			std::cout << std::setw(12) << std::left << x_(i + 1);
			for (unsigned int k = 0; k < nv_; k++)
				std::cout << std::setw(12) << std::left << ratios_(i, k);
			std::cout << std::endl;
		}
		std::cout << std::endl;
	}

	void LookupTable::Interpolation(const double x)
	{
		if (x <= min_x_)
		{
			for (unsigned int k = 0; k < nv_; k++)
				interpolated_(k) = v_(0, k);
			return;
		}
		else if (x >= max_x_)
		{
			for (unsigned int k = 0; k < nv_; k++)
				interpolated_(k) = v_(np_ - 1, k);
			return;
		}
		else
		{
			for (unsigned int i = 1; i < np_; i++)
				if (x_(i) >= x)
				{
					for (unsigned int k = 0; k < nv_; k++)
						interpolated_(k) = v_(i - 1, k) + ratios_(i - 1, k)*(x - x_(i - 1));
					return;
				}
		}
	}

	double LookupTable::Interpolation(const double x, const unsigned int k)
	{
		if (x <= min_x_)
		{
			return v_(0, k);
		}
		else if (x >= max_x_)
		{
			return v_(np_ - 1, k);
		}
		else
		{
			for (unsigned int i = 1; i < np_; i++)
				if (x_(i) >= x)
				{
					return (v_(i - 1, k) + ratios_(i - 1, k)*(x - x_(i - 1)));
				}
		}
	}
}

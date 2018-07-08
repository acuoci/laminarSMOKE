/*-----------------------------------------------------------------------*\
|                                                                         |
|   Authors: G. D'Alessio, A. Cuoci, A. Parente                           |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE solver.                             |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2018 A. Cuoci                                            |
|   laminarSMOKE is free software: you can redistribute it and/or modify  |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   laminarSMOKE is distributed in the hope that it will be useful,       |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE. If not, see <http://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

#ifndef ClassifierVQ2_HPP
#define ClassifierVQ2_HPP

#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <math.h>
#include <numeric>
#include <armadillo>

arma::Mat<double> GetMatrix(const std::string file_name);

class ClassifierVQ2
{
	public:
	
		void Setup(const std::string folder_name);
		std::vector<int> Classify(arma::Mat<double> Y);

		unsigned int nc() const { return nc_; }
		unsigned int nv() const { return nv_; }
	
	private:

		void FatalErrorMessage(const std::string message);
	
		unsigned int nv_;
		unsigned int nc_;
		
		arma::Mat<double> C_;
		arma::Mat<double> D_;
		arma::Mat<double> mu_;
		arma::Mat<double> sigma_;
		
		std::vector< arma::Mat<double> > eig_; 
		std::vector< arma::Mat<double> > eig_tran_; 
};

void ClassifierVQ2::Setup(const std::string folder_name)
{
	std::string file_name_C = folder_name + "/" + "centroids.txt";
	std::string file_name_D = folder_name + "/" + "D.txt";
	std::string file_name_mu = folder_name + "/" + "mu.txt";
	std::string file_name_sigma = folder_name + "/" + "sigma.txt";	
	
	C_ = GetMatrix(file_name_C.c_str());
	nv_ = C_.n_cols;
	nc_ = C_.n_rows;

	std::cout << "Number of variables (nv): " << nv_ << std::endl;
	std::cout << "Number of clusters (nc):  " << nc_ << std::endl;

    	D_ = GetMatrix(file_name_D.c_str());
	if (D_.n_rows != nv_ || D_.n_cols != nv_)
		FatalErrorMessage("Wrong D matrix size");

    	mu_ = GetMatrix(file_name_mu.c_str());
    	sigma_ = GetMatrix(file_name_sigma.c_str());
	
	eig_.resize(nc_);
	eig_tran_.resize(nc_);	
	for(int j = 0; j < nc_; j++)
	{
		std::cout << "Reading eig matrices: " << j+1 << std::endl;

		const std::string filename_eig = folder_name + "/" + "eigv" + std::to_string(j+1) + ".txt";
		eig_[j] = GetMatrix(filename_eig);
		if (eig_[j].n_rows != nv_ || eig_[j].n_cols != nv_)
			FatalErrorMessage("Wrong eig matrix size");

		eig_tran_[j].resize(nv_, nv_);
		eig_tran_[j] = eig_[j].t();
	}
}

std::vector<int> ClassifierVQ2::Classify(arma::Mat<double> Y)
{
	const unsigned int n = Y.n_rows;

	for(int i=0; i<n; i++)
	{
		Y.row(i) -= mu_.row(0);

		for(int j =0; j<Y.n_cols; j++)
			Y(i,j) /= sigma_(0,j);
	}
	
	arma::Mat<double> rec_err_os(n, nc_, arma::fill::zeros);
	arma::Mat<double> Y0(n, nv_, arma::fill::zeros);
	arma::Mat<double> DIST(n, nc_, arma::fill::zeros);

	for(int j = 0; j < nc_; j++)
	{
		for(int i = 0; i < n; i++)
			Y0.row(i) = Y.row(i) - C_.row(j);
	
		rec_err_os = (Y0 - Y0 * D_*eig_[j]*eig_tran_[j]*D_);

		DIST.col(j) = sqrt(sum(pow(rec_err_os, 2), 1));
	}

	std::vector<int> idx(n);
	for(int i=0; i<n; i++)
		idx[i] = DIST.row(i).index_min();

	return idx;
}

void ClassifierVQ2::FatalErrorMessage(const std::string message)
{
	std::cout << "ClassifierVQ2 Fatal error message" << std::endl;
	std::cout << message << std::endl;
	exit(-1);
}

arma::Mat<double> GetMatrix(const std::string file_name)
{
	std::string s;
	std::vector<std::vector<double>> matrix;

	std::ifstream fInput(file_name.c_str(), std::ios::in);
	while ( std::getline(fInput, s) )
	{
		matrix.push_back(std::vector<double>());
		std::vector<double>& vBack = matrix.back();
		std::istringstream ss(s);
		std::copy(std::istream_iterator<double>(ss), std::istream_iterator<double>(), std::back_inserter(vBack));
	}
	
	fInput.close();

	arma::Mat<double> matrix_mat(matrix.size(), matrix[0].size(), arma::fill::zeros);
	for(int i=0; i<matrix.size(); i++)
	{
		for(int j=0; j<matrix[0].size(); j++)
		{
			matrix_mat(i,j) = matrix[i][j];
		}
	}
	return matrix_mat;
}

#endif // ClassifierVQ2_HPP

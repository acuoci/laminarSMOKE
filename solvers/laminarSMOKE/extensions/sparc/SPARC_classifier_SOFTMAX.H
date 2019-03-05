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

#ifndef ClassifierSOFTMAX_HPP
#define ClassifierSOFTMAX_HPP

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

arma::Mat<double> CenterScale(arma::Mat<double> X);

class ClassifierSOFTMAX
{
	public:
	
		void Setup(const std::string folder_name);
		void Classify(arma::Mat<double> Y, std::vector<int>& idx);

		unsigned int nc() const { return nc_; }
		unsigned int nv() const { return nv_; }
	
	private:

		void FatalErrorMessage(const std::string message);
	
		unsigned int nv_;
		unsigned int nc_;
	
		arma::Mat<double> weightsMatrix_;
		
};

void ClassifierSOFTMAX::Setup(const std::string folder_name)
{
	std::string filename = folder_name + "/" + "weightsMatrix.txt";        //txt containing the idx from the training

	std::cout << "Reading weights matrix.." << std::endl;
	weightsMatrix_ = GetMatrix<double>(filename); 
	std::cout << "Weight matrix loaded successfully" << std::endl;

	nv_ = weightsMatrix_.n_cols;
	nc_ = weightsMatrix_.n_rows-1;

	std::cout << "Number of variables (nv): " << nv_ << std::endl;
	std::cout << "Number of clusters (nc):  " << nc_ << std::endl;
}

void ClassifierSOFTMAX::Classify(arma::Mat<double> Y, std::vector<int>& idx)
{    
	int n_classes = weightsMatrix_.n_rows;
	if(Y.n_rows>Y.n_cols)
		Y = Y.t();

	if(arma::max(Y.row(0.)) != 1.)
		Y = CenterScale(Y);

	// Classification
	{	    
		arma::Mat<double> scores = weightsMatrix_ * Y;
		arma::Mat<double> temp = arma::repmat(arma::max(scores, 0.), n_classes, 1);
		
		scores = scores - temp;
		scores = arma::exp(scores);
		temp = arma::repmat(arma::sum(scores, 0.), n_classes, 1);
		scores = scores / temp;
		scores = arma::log(scores);
		arma::rowvec result = arma::zeros<arma::rowvec>(scores.n_cols);

		for(int i=0; i<scores.n_cols; i++)
		{
			// WARNING: -1 BECAUSE THE LPCA TRAINING STARTS FROM 1, WHILE THE MECHANISMS START FROM 0.
			idx[i] = scores.col(i).index_max() -1; 
		}
	}
}

void ClassifierSOFTMAX::FatalErrorMessage(const std::string message)
{
	std::cout << "ClassifierSOFTMAX Fatal error message" << std::endl;
	std::cout << message << std::endl;
	exit(-1);
}

arma::Mat<double> CenterScale(arma::Mat<double> X)
{
    arma::Mat<double> X0(X.n_rows, X.n_cols);
    arma::Mat<double> rowave(X.n_rows, 1);
    arma::Mat<double> d(X.n_rows,1);

    for(int i =0; i<X.n_rows; i++)
        rowave.row(i) = arma::min(X.row(i));

    for(int i =0; i<X.n_rows; i++)
        d.row(i) = arma::max(X.row(i)) - arma::min(X.row(i)) + 1e-16;

    for(int i =0; i<X.n_cols; i++)
        X0.col(i) = X.col(i) - rowave;

    for(int i =0; i<X.n_cols; i++)
        for(int j=0; j<X.n_rows; j++)
            X0(j,i) = X0(j,i)/d(j,0);

    return X0;
}

#endif // ClassifierSOFTMAX_HPP

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

arma::Mat<double> CenterScale(const arma::Mat<double>& X, const int scaling);

arma::field<arma::Mat<double>> GetMuSigma(const arma::Mat<double>& matrix, const int scaling);

arma::field <arma::Mat<double>> GetCluster(const arma::Mat<double>& matrix, const arma::Mat<int>& idx);

arma::Mat<double> GetCentroids(const arma::field<arma::Mat<double>>& clusters);

arma::field<arma::Mat<double>> GetModes(const arma::field<arma::Mat<double>>& clusters, const int eig);

template<typename T>
arma::Mat<T> GetMatrix(const std::string file_name);

class ClassifierVQ2
{
	public:
	
		void Setup(const std::string folder_name, const int eig, const int scaling);
		void Classify(arma::Mat<double> Y, std::vector<int>& idx);

		unsigned int nc() const { return nc_; }
		unsigned int nv() const { return nv_; }
	
	private:

		void FatalErrorMessage(const std::string message);
	
		unsigned int nv_;
		unsigned int nc_;
		
		arma::Mat<double> C_;
		arma::field<arma::Mat<double>> MuSigma_;
		arma::field<arma::Mat<double>> Modes_; 

	        arma::Mat<double> Y0_;
	        arma::Mat<double> DIST_;
	        arma::Mat<double> rec_err_os_;
};

void ClassifierVQ2::Setup(const std::string folder_name, const int eig, const int scaling)
{
    // eig: number of eigens used for the training 
    // scaling: 0 = auto, 1 = vast, 2 = pareto

    std::string file_name_idx = folder_name + "/" + "idx.txt";
    std::string file_name_X = folder_name + "/" + "X.txt";

    std::cout << "Reading idx..." << std::endl;                                                  
    arma::Mat<int> idx = GetMatrix<int>(file_name_idx);
    std::cout << " * number of elements: " << idx.size() << std::endl;
    std::cout << " * minimum index: " << idx.min() << std::endl;
    std::cout << " * maximum index: " << idx.max() << std::endl;

    std::cout << "Reading X matrix..." << std::endl;                                                  
    arma::Mat<double> X = GetMatrix<double>(file_name_X); 
    std::cout << " * number of elements:  " << X.n_rows << std::endl;
    std::cout << " * number of variables: " << X.n_cols << std::endl;
    
    std::cout << "Calculating mu and sigma from training matrix..." << std::endl;
    MuSigma_ = GetMuSigma(X, scaling);
	
    std::cout << "Centering and scaling the training matrix..." << std::endl;
    arma::Mat<double> centered_scaled_X = CenterScale(X, scaling);
    
    std::cout << "Calculating clusters.." << std::endl;
    arma::field<arma::Mat<double>> clusters = GetCluster(centered_scaled_X, idx);
    
    std::cout << "Calculating centroids..." << std::endl;
    C_ = GetCentroids(clusters);                                                
    std::cout << "Centroids " << std::endl;

    std::cout<<"Calculating modes.." << std::endl;
    Modes_ = GetModes(clusters, eig);

    nv_ = C_.n_cols;
    nc_ = C_.n_rows;

    std::cout << "Number of variables (nv): " << nv_ << std::endl;
    std::cout << "Number of clusters (nc):  " << nc_ << std::endl;

    // Check consistency of data
   if (idx.size() != X.n_rows)
	FatalErrorMessage("The number of elements in idx vector does not match the number of elements in the X matrix");
    if (nv_ != X.n_cols)
	FatalErrorMessage("The number of variables in X matrix does not match the number of variables in C matrix");
    if (nc_ != idx.max())
	FatalErrorMessage("The number of clusters from the idx vector does not match the number of cluster from C matrix");
}

void ClassifierVQ2::Classify(arma::Mat<double> Y, std::vector<int>& idx)
{    
    const unsigned int n = Y.n_rows;

    // Memory allocation carried out only in case n changes
    if (Y0_.n_rows != n)
    {
	    Y0_.resize(n, nv_);
	    DIST_.resize(n, nc_);
	    rec_err_os_.resize(n, nc_);
    }

    arma::Row<double> mu = MuSigma_(0).row(0);
    arma::Mat<double> sigma = MuSigma_(1);
    
   for(int i=0; i<n; i++)
        Y.row(i) = Y.row(i) - mu;

    for(int i =0; i<n; i++)
        for(int j =0; j<nv_; j++)
            Y(i,j) /= sigma(0,j);

    for(int j = 0; j < C_.n_rows; j++)
    {
		for(int i = 0; i < n; i++)
			Y0_.row(i) = Y.row(i)-C_.row(j);
			 
		rec_err_os_ = (Y0_ - Y0_ * Modes_(j) * Modes_(j).t());
		DIST_.col(j) = arma::sum(arma::pow(rec_err_os_, 2), 1);
    }

    idx.resize(n);
    for(int i=0; i<n; i++)
        idx[i] = DIST_.row(i).index_min();
}

void ClassifierVQ2::FatalErrorMessage(const std::string message)
{
	std::cout << "ClassifierVQ2 Fatal error message" << std::endl;
	std::cout << message << std::endl;
	exit(-1);
}

arma::field<arma::Mat<double>> GetMuSigma(const arma::Mat<double>& matrix, const int scaling)
{
    const int n_oss = matrix.n_rows;
    const int n_variab = matrix.n_cols;
    
    arma::Mat<double> mediarow(1, n_variab);
    for(int j=0; j<n_variab; j++)
		mediarow.col(j) = arma::mean(arma::mean(matrix.col(j)));
    
    arma::Mat<double> d(1, n_variab);
    if(scaling == 0)
    {
        for(int j=0; j<n_variab; j++)
            d.col(j) = arma::stddev(matrix.col(j)) + 1e-16; 
    } 
    else if(scaling == 1)
    {
        for(int j=0; j<n_variab; j++)
            d.col(j) = std::pow(arma::stddev(matrix.col(j)),2.)/(mediarow.col(j)) + 1e-16; 
    } 
    else if(scaling == 2)
    {
        for(int j=0; j<n_variab; j++)
            d.col(j) = arma::stddev(matrix.col(j)) + 1e-16; 
    }

    arma::field<arma::Mat<double>> MuSigma(2,1);
    MuSigma(0,0) = mediarow;
    MuSigma(1,0) = d;

    return MuSigma;
}

arma::Mat<double> CenterScale(const arma::Mat<double>& X, const int scaling)
{
    // Declarations
    arma::Mat<double> X0(X.n_rows, X.n_cols);
    arma::Mat<double> rowave(1, X.n_cols);
    arma::Mat<double> d(1, X.n_cols);
	
    // Get variables' mean values
    for(int j=0; j<X.n_cols; j++)
        rowave.col(j) = arma::mean(X.col(j));
    
    if(scaling == 0)
    {
        for(int j=0; j<X.n_cols; j++)
            d.col(j) = arma::stddev(X.col(j)) + 1e-16; 
    } 
    else if(scaling == 1)
    {
        for(int j=0; j<X.n_cols; j++)
            d.col(j) = std::pow(arma::stddev(X.col(j)),2)/(rowave.col(j)) + 1e-16; 
    } 
    else if(scaling == 2)
    {
        for(int j=0; j<X.n_cols; j++)
            d.col(j) = arma::stddev(X.col(j)) + 1e-16; 
    }
    
    // Subtract the mean and divide by the scaling factor
    for(int i=0; i<X.n_rows; i++)
        X0.row(i) = X.row(i) - rowave;
    
    for(int i=0; i<X.n_rows; i++)
        for(int j=0; j<X.n_cols; j++)
            X0(i,j) /= d(0,j);
        
    return X0;
}

arma::field<arma::Mat<double>> GetCluster(const arma::Mat<double>& matrix, const arma::Mat<int>& idx)
{      
    const int n_var = matrix.n_cols;
    const int max_idx = idx.max();

    arma::field <arma::Mat<double> > clusters(max_idx, 1);

    for(int i = 1; i<=max_idx; i++)
    {
        arma::uvec q = arma::find(idx == i);
        clusters(i-1) = matrix.rows(q);
    }
    return clusters;
}

arma::Mat<double> GetCentroids(const arma::field<arma::Mat<double>>& clusters)
{                    
    const int n_clust = clusters.n_rows;
    const int n_var = clusters(1).n_cols;
    arma::Mat<double> C(n_clust, n_var);

    for(int i=0; i<clusters.n_rows; i++)
    {
        arma::Mat<double> X = clusters(i);
        for(int j=0; j<X.n_cols; j++)
            C(i,j) = mean(X.col(j));
    }
    return C;
}

arma::field<arma::Mat<double>> GetModes(const arma::field<arma::Mat<double>>& clusters, const int eig)
{
    arma::field<arma::Mat<double>> Modes(clusters.n_rows, 1);

/*
    for(int i=0; i<clusters.n_rows; i++)
    {
	arma::mat Cov = (clusters(i).t() * clusters(i)) / static_cast<double>(clusters(i).n_rows);
		
        arma::mat U;
        arma::vec s;
        arma::mat V;

        svd(U,s,V,Cov);
        Modes(i,0) = U.cols(0,eig);
    }
*/
    for(int i=0; i<clusters.n_rows; i++)
    {
	arma::mat Cov = (clusters(i).t() * clusters(i)) / static_cast<double>(clusters(i).n_rows-1);
		
        arma::vec eigval;
        arma::mat eigvec;

        eig_sym(eigval, eigvec, Cov);
	const int max_mod = eigvec.n_cols;
        Modes(i,0) = eigvec.cols(max_mod-eig-1, max_mod-1);
    }

    return Modes;
}

template<typename T>
arma::Mat<T> GetMatrix(const std::string file_name)
{
	std::vector<std::vector<T>> matrix;

	std::ifstream pFile(file_name);
	std::string s;
	
	while ( std::getline(pFile, s) )
	{
		matrix.push_back(std::vector<T>());
		std::vector<T>& vBack = matrix.back();
		std::istringstream ss(s);
		std::copy(std::istream_iterator<T>(ss), std::istream_iterator<T>(), std::back_inserter(vBack));
	}

       pFile.close();
	
	arma::Mat<T> matrix_mat(matrix.size(), matrix[0].size(), arma::fill::zeros);
	for(int i=0; i<matrix.size(); i++)
		for(int j=0; j<matrix[0].size(); j++)
			matrix_mat(i,j) = matrix[i][j];

	return matrix_mat;
}

#endif // ClassifierVQ2_HPP

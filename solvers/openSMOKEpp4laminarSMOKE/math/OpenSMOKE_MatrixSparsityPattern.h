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
|   Copyright(C) 2016, 2015, 2014, 2013, 2012  Alberto Cuoci              |
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

#ifndef MatrixSparsityPattern_H
#define MatrixSparsityPattern_H

#include "math/OpenSMOKEVector.h"

namespace OpenSMOKE
{
	struct ElementSparsityPattern
	{
		int column;
		ElementSparsityPattern *next;
	};

	//!  A class to manage the sparsity pattern of a matrix
	/*!
	The purpose of this class is to manage the sparsity pattern of a matrix
	In particular, the class is able to find the dependencies among the variables. This information can be 
	exploited for smart calculation of sparse Jacobian matrices, by minimizing the number of 
	evaluations of the functions
	*/

	class OpenSMOKE_MatrixSparsityPattern
	{

	public:

		/**
		*@brief Default constructor (no memory allocation)
		*/
		OpenSMOKE_MatrixSparsityPattern();

		/**
		*@brief Constructor with memory allocation
		*@param rows the number of rows
		*@param cols the number of cols
		*/
		OpenSMOKE_MatrixSparsityPattern(const int rows, const int cols);

		/**
		*@brief Destructor
		*/
		~OpenSMOKE_MatrixSparsityPattern();

		/**
		*@brief Removes the internal variables
		*/
		void Delete();

		/**
		*@brief Resizes the matrix
		*@param rows the number of rows
		*@param cols the number of cols
		*/
		void resize(const int rows, const int columns);

		/**
		*@brief Sets the non-zero element
		*@param row index of row of non-zero element
		*@param col index of column of non-zero element
		*/
		void operator ()(const int row, const int col);

		/**
		*@brief Analyzes the pattern to find the dependencies among the variables
		*/
		void FindDependence();

		/**
		*@brief Calculates and returns the total number of non-zero elements
		*/
		int CountElements();

		/**
		*@brief Reset the scanning operation
		*/
		void ResetScanning();

		/**
		*@brief Scanning procedure: returns the indices (row and column) of the current non-zero element
		*/
		int Scanning(int *i, int *j);

		/**
		*@brief Assignement operator: A=B
		*@param rhs the matrix to be copied
		*/
		OpenSMOKE_MatrixSparsityPattern& operator = (const OpenSMOKE_MatrixSparsityPattern &rval);

		/**
		*@brief Returns the number of groups which can be recognized
		*/
		inline int number_groups() const { return number_groups_; }

		/**
		*@brief Returns the total number of non-zero elements
		*/
		inline int number_elements() const { return number_elements_; }

		/**
		*@brief Returns the number of variables available in the requested group
		*@param j index of group
		*/
		inline int number_variables_in_group(const int j) const { return number_variables_in_group_[j]; }

		/**
		*@brief Number of equationsin which the specified variable is involved
		*@param j index of variable
		*/
		inline int number_equations_per_variable(const int j) const { return number_equations_per_variable_[j]; }

		/**
		*@brief Returns the index of k-th variable in the j-th group
		*@param group index of group
		*@param variable index of variable
		*/
		inline int variable_in_group(const int group, const int variable) const { return variable_in_group_[group][variable]; }

		/**
		*@brief Returns the dependencies of the specified variable
		*@param j index of variable
		*@param equation index of the equation
		*/
		inline int dependencies(const int j, const int equation) const { return dependencies_[j][equation]; }

		/**
		*@brief Returns the number of rows
		*/
		inline int rows() const { return rows_; }

		/**
		*@brief Returns the number of columns
		*/
		inline int cols() const { return rows_; }

		/**
		*@brief Returns the upper diagonal dimension
		*/
		inline int up() const { return up_; }

		/**
		*@brief Returns the lower diagonal dimension
		*/
		inline int low() const { return low_; }

		/**
		*@brief Returns the density with respect to the corresponding band matrix
		*/
		inline double density() const { return density_; }

		/**
		*@brief Returns the maximum number of non-zero elements which are present in a single column
		*/
		inline int max_elements_in_cols() const { return max_elements_in_cols_; }

		/**
		*@brief Returns the minimum number of non-zero elements which are present in a single column
		*/
		inline int min_elements_in_cols() const { return min_elements_in_cols_; }

		/**
		*@brief Returns the maximum number of non-zero elements which are present in a single row
		*/
		inline int max_elements_in_rows() const { return max_elements_in_rows_; }

		/**
		*@brief Returns the minimum number of non-zero elements which are present in a single column
		*/
		inline int min_elements_in_rows() const { return min_elements_in_rows_; }

		/**
		*@brief Returns the number of zeros along the main diagonal
		*/
		inline int number_zeros_on_diagonal() const { return number_zeros_on_diagonal_; }

		/**
		*@brief Returns the index of the (first) column in which there is the maximum number of non-zero elements
		*/
		inline int index_max_col() const { return index_max_col_; }

		/**
		*@brief Returns the index of the (first) row in which there is the maximum number of non-zero elements
		*/
		inline int index_max_row() const { return index_max_row_; }

		/**
		*@brief Returns the index of the (first) column in which there is the minimum number of non-zero elements
		*/
		inline int index_min_col() const { return index_min_col_; }

		/**
		*@brief Returns the index of the (first) row in which there is the minimum number of non-zero elements
		*/
		inline int index_min_row() const { return index_min_row_; }

	private:

		/**
		*@brief Allocates the memory
		*/
		void Initialize(const int rows, const int cols);

		/**
		*@brief Inserts a new non-zero element
		*/
		void InsertElement(ElementSparsityPattern *elem, const int row, const int col, const int first);

		/**
		*@brief Copy a given pattern into the current object
		*/
		void Copy(const OpenSMOKE_MatrixSparsityPattern &rval);

		/**
		*@brief Analyzes the pattern to find the lower and upper bands
		*/
		void AnalyzeBand();

		/**
		*@brief Calculated the density of the sparsity pattern compared with a band matrix
		*/
		void CalculateDensity();

		/**
		*@brief Remove the previous dependencies analysis (if any)
		*/
		void DeleteDependence();

		/**
		*@brief Remove all the elements in the matrix
		*/
		void RemoveAllElementsInMatrix();

	private:

		int rows_;							//!< number of rows
		int columns_;						//!< number of columns
		int low_;							//!< lower band
		int up_;							//!< upper band
		double density_;					//!< density with respect to the band matrix

		int	number_groups_;					//!< number of groups
		int	number_elements_;				//!< number of elements
		int number_zeros_on_diagonal_;		//!< number of non-zero elements along the main diagonal
		int max_number_of_elements_;		//!< maximum number of elements (rows*columns)

		int max_elements_in_cols_;			//!< maximum number of non-zero elements which can be found in a single column
		int min_elements_in_cols_;			//!< minimum number of non-zero elements which can be found in a single column
		int max_elements_in_rows_;			//!< maximum number of non-zero elements which can be found in a single row
		int min_elements_in_rows_;			//!< minimum number of non-zero elements which can be found in a single row
		int index_max_col_;					//!< index of column with the maximum number of non-zero elements
		int index_max_row_;					//!< index of row with the maximum number of non-zero elements
		int index_min_col_;					//!< index of column with the minimum number of non-zero elements
		int index_min_row_;					//!< index of row with the minimum number of non-zero elements

		bool is_dependence_available_;		//!< internal flag: true means that the dependencies analysis has been performed

		ElementSparsityPattern **element_row_;
		ElementSparsityPattern	*element_scanning_;

		int	**dependencies_;					//!< list of dependencies among the variables
		int	row_scanning_;						//!< current row under scanning analysis
		int *ptr_derivatives_;

		int	*number_equations_per_variable_;	//!< number of equations per variable
		int	**variable_in_group_;				//!< list of indices of variables per group
		int	*number_variables_in_group_;		//!< number of variables per group
	};

}

#include "OpenSMOKE_MatrixSparsityPattern.hpp"

#endif // OpenSMOKE_MatrixSparsityPattern_H

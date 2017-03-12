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
|   Copyright(C) 2016  Alberto Cuoci                                      |
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
	OpenSMOKE_MatrixSparsityPattern::OpenSMOKE_MatrixSparsityPattern()
	{
		Initialize(0, 0);
	}

	OpenSMOKE_MatrixSparsityPattern::OpenSMOKE_MatrixSparsityPattern(const int rows, const int cols)
	{
		Initialize(rows, cols);
	}

	void OpenSMOKE_MatrixSparsityPattern::resize(const int rows, const int columns)
	{
		RemoveAllElementsInMatrix();

		if (columns != columns_)
			columns_ = columns;
		
		if (rows != rows_)
		{
			delete element_row_;
			Initialize(rows, columns);
		}
	}

	void OpenSMOKE_MatrixSparsityPattern::Initialize(const int nrows, const int ncols)
	{
		dependencies_ = 0;
		variable_in_group_ = 0;
		number_equations_per_variable_ = 0;
		number_variables_in_group_ = 0;
		ptr_derivatives_ = 0;
		is_dependence_available_ = false;
		number_groups_ = 0;
		row_scanning_ = 1;
		element_scanning_ = 0;

		if (nrows < 1 || ncols < 1)
		{
			rows_ = columns_ = 0;
			element_row_ = 0;
			return;
		}

		if (nrows < 0 || ncols < 0)
			OpenSMOKE::FatalErrorMessage("MatrixSparsityPattern: the number of rows and columns must be positive or equal to 0");

		rows_ = nrows;
		columns_ = ncols;
		element_row_ = new ElementSparsityPattern *[rows_ + 1];

		if (!element_row_)
			OpenSMOKE::FatalErrorMessage("MatrixSparsityPattern: it was not possible to allocate enough memory to initialize the matrix");

		for (int row = 0; row <= rows_; row++)
			element_row_[row] = 0;

		low_ = 0;
		up_ = 0;
	}

	void OpenSMOKE_MatrixSparsityPattern::operator ()(const int row, const int col)
	{
		if (row < 1 || row > rows_ || col < 1 || col > columns_)
		{
			std::cout << "Matrix size: " << rows_ << " x " << columns_ << std::endl;
			std::cout << "Requested element: " << row << " " << col << std::endl;
			OpenSMOKE::FatalErrorMessage("MatrixSparsityPattern: the row and column indices of inserted element must fit the matrix dimension");
		}
		ElementSparsityPattern *elem = element_row_[row];
		if (elem == 0 || col < elem->column)
		{
			InsertElement(elem, row, col, 1);
			return;
		}

		if (col == elem->column)
			return;

		for (; elem->next != 0; elem = elem->next)
		{
			if (col == elem->next->column)
				return;
			else if (col < elem->next->column)
				break;
		}

		InsertElement(elem, row, col, 0);
	}

	void OpenSMOKE_MatrixSparsityPattern::InsertElement(ElementSparsityPattern *elem, const int row, const int col, const int first)
	{
		if (is_dependence_available_ == true)
			DeleteDependence();

		ElementSparsityPattern *newElement = new ElementSparsityPattern;

		if (!newElement)
			OpenSMOKE::FatalErrorMessage("MatrixSparsityPattern: it was not possible to allocate enough memory to insert a new element");

		newElement->column = col;
		if (first == 1)
		{
			newElement->next = elem;
			element_row_[row] = newElement;
		}
		else
		{
			newElement->next = elem->next;	// inserted next
			elem->next = newElement;
		}
	}

	void OpenSMOKE_MatrixSparsityPattern::DeleteDependence()
	{
		if (is_dependence_available_ == false)
			return;

		is_dependence_available_ = false;

		delete number_equations_per_variable_;
		delete dependencies_[0];
		delete dependencies_;
		delete ptr_derivatives_;
		delete variable_in_group_;
		delete number_variables_in_group_;
	}

	void OpenSMOKE_MatrixSparsityPattern::FindDependence()
	{
		const int max_length_int_vector = static_cast<int> (std::pow(2., 32.) - 1);

		if (is_dependence_available_ == true)
			return;

		is_dependence_available_ = true;

		ElementSparsityPattern *elem;
		number_elements_ = 0;
		max_elements_in_rows_ = 0;
		min_elements_in_rows_ = columns_;
		max_elements_in_cols_ = 0;
		min_elements_in_cols_ = rows_;

		number_equations_per_variable_ = new int[columns_ + 1];
		memset(number_equations_per_variable_, 0, (columns_ + 1)*sizeof(int));

		int maxRow, minRow;
		OpenSMOKE::OpenSMOKEVectorInt elementsInColumn(columns_);
		number_zeros_on_diagonal_ = 0;
		for (int row = 1; row <= rows_; row++)
		{
			maxRow = 0;
			minRow = 0;
			elem = element_row_[row];
			while (elem != 0)
			{
				elementsInColumn(elem->column)++;
				if (elem->column == row)
					number_zeros_on_diagonal_++;
				maxRow++;
				minRow++;
				number_elements_++;
				number_equations_per_variable_[elem->column]++;
				if (number_elements_ == max_length_int_vector)
					OpenSMOKE::FatalErrorMessage("MatrixSparsityPattern::FindDependence(): Too many element dependencies");

				elem = elem->next;
			}
			if (maxRow > max_elements_in_rows_)
			{
				max_elements_in_rows_ = maxRow;
				index_max_row_ = row;
			}
			if (minRow < min_elements_in_rows_)
			{
				min_elements_in_rows_ = minRow;
				index_min_row_ = row;
			}
		}

		max_elements_in_cols_ = elementsInColumn.Max(&index_max_col_);
		min_elements_in_cols_ = elementsInColumn.Min(&index_min_col_);
		number_zeros_on_diagonal_ = Min(rows_, columns_) - number_zeros_on_diagonal_;

		dependencies_ = new int *[columns_ + 1];
		dependencies_[0] = new int[number_elements_ + 1];
		dependencies_[1] = dependencies_[0];
		int *ptrInt = dependencies_[0] + 1;
		int *ivar = new int[columns_ + 1];
		memset(ivar, 0, (columns_ + 1)*sizeof(int));

		for (int var = 1; var < columns_; var++)
			dependencies_[var + 1] = dependencies_[var] + number_equations_per_variable_[var];
		for (int row = 1; row <= rows_; row++)
		{
			elem = element_row_[row];
			while (elem != 0)
			{
				ivar[elem->column]++;
				dependencies_[elem->column][ivar[elem->column]] = row;
				elem = elem->next;
			}
		}
		delete ivar;
		ptrInt = dependencies_[0] + 1;
		char *tempFunz = new char[rows_ + 1];
		char *tempVar = new char[columns_ + 1];
		memset(tempVar, 0, (columns_ + 1)*sizeof(char));
		ptr_derivatives_ = new int[columns_ + 1];
		int *tempDer = new int[columns_ + 1];

		int k;
		int ider = 1;
		int jder = 1;
		int kder;
		number_groups_ = 0;

		int var;
		while (1)
		{
			memset(tempFunz, 0, (rows_ + 1)*sizeof(char));

			for (var = 1; var <= columns_; var++)
			if (tempVar[var] == 0)
				break;

			if (var > columns_)
				break;

			tempVar[var] = 1;
			ptr_derivatives_[ider++] = var;
			number_groups_++;
			kder = 1;
			for (k = 1; k <= number_equations_per_variable_[var]; k++)
				tempFunz[dependencies_[var][k]] = 1;
			while (1)
			{
				for (var += 1; var <= columns_; var++)
				if (tempVar[var] == 0)break;
				if (var > columns_)
				{
					tempDer[jder++] = kder;
					break;
				}
				char no = 0;
				for (k = 1; k <= number_equations_per_variable_[var]; k++)
				{
					if (tempFunz[dependencies_[var][k]] == 1)
					{
						no = 1; break;
					}
				}
				if (no == 0)
				{
					for (k = 1; k <= number_equations_per_variable_[var]; k++)
						tempFunz[dependencies_[var][k]] = 1;
					tempVar[var] = 1;
					ptr_derivatives_[ider++] = var;
					kder++;
				}
			}
		}
		variable_in_group_ = new int *[number_groups_ + 1];
		number_variables_in_group_ = new int[number_groups_ + 1];
		var = 0;
		for (k = 1; k <= number_groups_; k++)
		{
			number_variables_in_group_[k] = tempDer[k];
			variable_in_group_[k] = &ptr_derivatives_[var];
			var += tempDer[k];
		}
		max_number_of_elements_ = rows_*columns_;

		delete tempFunz;
		delete tempVar;
		delete tempDer;

		// Additional analyses
		AnalyzeBand();
		CalculateDensity();
	}

	void OpenSMOKE_MatrixSparsityPattern::CalculateDensity()
	{
		density_ = 1.;
		density_ = double((rows_ - low_)*low_ + (rows_ - up_)*up_);
		if (density_ != 0.)
			density_ = double(this->number_elements()) / density_;
	}

	int OpenSMOKE_MatrixSparsityPattern::CountElements()
	{
		int countElements = 0;

		if (element_row_ == 0)
			return 0;

		ElementSparsityPattern *elem;
		for (int row = 1; row <= rows_; row++)
		{
			elem = element_row_[row];
			while (elem != 0)
			{
				countElements++;
				elem = elem->next;
			}
		}

		return countElements;
	}

	void OpenSMOKE_MatrixSparsityPattern::AnalyzeBand()
	{
		low_ = 0;
		up_ = 0;

		if (element_row_ == 0)
			return;

		ElementSparsityPattern *elem;

		int i;
		for (int row = 1; row <= rows_; row++)
		{
			elem = element_row_[row];
			i = 0;
			while (elem != 0)
			{
				i = row - elem->column;
				if (i > low_)
					low_ = i;
				i = elem->column - row;
				if (i > up_)
					up_ = i;
				elem = elem->next;
			}
		}
	}

	void OpenSMOKE_MatrixSparsityPattern::ResetScanning()
	{
		row_scanning_ = 1;
		element_scanning_ = 0;
	}

	int OpenSMOKE_MatrixSparsityPattern::Scanning(int *i, int *j)
	{
		if (element_row_ == 0)
			return 0;

		if (row_scanning_ == 0)
		{
			*i = *j = 0;
			row_scanning_ = 1;
			element_scanning_ = element_row_[1];
			return 0;
		}

		if (row_scanning_ == 1 && element_scanning_ == 0)
		{
			for (; row_scanning_ <= rows_; row_scanning_++)
			{
				element_scanning_ = element_row_[row_scanning_];
				if (element_scanning_ != 0)break;
			}
		}

		if (row_scanning_ > rows_)
			return 0;

		*i = row_scanning_;
		*j = element_scanning_->column;
		element_scanning_ = element_scanning_->next;
		do
		{
			if (element_scanning_ == 0)
			{
				row_scanning_++;
				if (row_scanning_ > rows_)
				{
					row_scanning_ = 0;
					element_scanning_ = 0;
					break;
				}
				else
					element_scanning_ = element_row_[row_scanning_];
			}
		} while (element_scanning_ == 0);

		return 1;
	}

	void OpenSMOKE_MatrixSparsityPattern::RemoveAllElementsInMatrix()
	{
		if (element_row_ == 0)
			return;

		if (is_dependence_available_ == true)
			DeleteDependence();

		ElementSparsityPattern *elem;
		ElementSparsityPattern *temp;
		for (int row = 1; row <= rows_; row++)
		{
			elem = element_row_[row];
			while (elem != 0)
			{
				temp = elem;
				elem = elem->next;
				delete temp;
			}
			element_row_[row] = 0;
		}
	}

	void OpenSMOKE_MatrixSparsityPattern::Copy(const OpenSMOKE_MatrixSparsityPattern &rval)
	{
		row_scanning_ = 1;
		element_scanning_ = 0;
		for (int row = 1; row <= rows_; row++)
		{
			ElementSparsityPattern *elem;
			ElementSparsityPattern *elemRval;

			int first = 1;
			for (elemRval = rval.element_row_[row]; elemRval != 0; elemRval = elemRval->next)
			{
				ElementSparsityPattern *newElement = new ElementSparsityPattern;

				if (!newElement)
					OpenSMOKE::FatalErrorMessage("MatrixSparsityPattern::FindDependence(): Too many element dependencies");

				newElement->column = elemRval->column;
				newElement->next = 0;
				if (first == 1)
				{
					elem = element_row_[row] = newElement;
					first = 0;
				}
				else
				{
					elem->next = newElement;
					elem = elem->next;
				}
			}
		}
	}

	OpenSMOKE_MatrixSparsityPattern &OpenSMOKE_MatrixSparsityPattern::operator = (const OpenSMOKE_MatrixSparsityPattern &rval)
	{
		RemoveAllElementsInMatrix();

		if (columns_ != rval.columns_)
			columns_ = rval.columns_;

		if (rows_ != rval.rows_)
		{
			delete element_row_;
			Initialize(rval.rows_, rval.columns_);
		}

		Copy(rval);

		return *this;
	}

	void OpenSMOKE_MatrixSparsityPattern::Delete()
	{
		if (element_row_ == 0)
		{
			return;
		}
		
		RemoveAllElementsInMatrix();
		delete element_row_;
		element_row_ = 0;
		rows_ = 0;
		columns_ = 0;
	}


	OpenSMOKE_MatrixSparsityPattern::~OpenSMOKE_MatrixSparsityPattern()
	{
		Delete();
		if (is_dependence_available_ == true)
			DeleteDependence();
	}
}


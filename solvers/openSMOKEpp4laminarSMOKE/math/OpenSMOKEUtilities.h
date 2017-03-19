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
|   License                                                               |
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

#ifndef OpenSMOKE_OpenSMOKEUtilities_Hpp
#define OpenSMOKE_OpenSMOKEUtilities_Hpp

#include "OpenSMOKEStdInclude.h"
#include "OpenSMOKEBaseClass.h"
#include <boost/filesystem.hpp>
#include "rapidxml.hpp"

namespace OpenSMOKE
{
	/**
	* Policy class for OpenSMOKEVector with zero index basis
	*/
	class ZeroIndexPolicy //: public OpenSMOKEBaseClass
	{
		protected:
			static const int index_ = 0;	/**< index basis */
	};

	/**
	* Policy class for OpenSMOKEVector with one index basis
	*/
	class OneIndexPolicy //: public OpenSMOKEBaseClass
	{
		protected:
			static const int index_ = 1;	/**< index basis */
	};

	template <typename T>
	void Swap(T* x, T* y);

	template <typename T>
	inline T Abs(T const& a);

	template <typename T>
	inline T const& Max(T const& a, T const& b);

	template <typename T>
	inline T const& Min(T const& a, T const& b);

	template <typename T>
	inline T MaxAbs(T const& a, T const& b);

	template <typename T>
	inline T MinAbs(T const& a, T const& b);

	template <typename T>
	inline T Max(const int n, T *x, int *im);

	template <class T>
	inline T Max(const int n, T *x);

	template <class T>
	inline T MaxAbs(const int n, T *x,int *im);

	template <class T>
	inline T MaxAbs(const int n, T *x);

	template <class T>
	inline T Min(const int n, T *x, int *im);

	template <class T>
	inline T Min(const int n, T *x);

	template <class T>
	inline T MinAbs(const int n, T *x, int *im);

	template <class T>
	inline T MinAbs(const int n, T *x);

	/**
	* The sum of two vectors with control
	* Example: Sum(n,x,y,z) means \f$ z=x+y \f$;
	*/	
	template <class T>
	void Sum(const int n, const T* lval, const T* rval, T* result);

	/**
	* The sum of a constant to a vector
	* Example: Sum(n,x,c,z) means \f$ z=x+c \f$;
	*/	
	template <class T>
	void Sum(const int n, const T* lval, const T rval, T* result);

	/**
	* Adds a constant value to a vector
	* Example: Sum(n,2,&x) means \f$ x_i=x_i+2 \f$;
	*/	
	template <class T>
	void Sum(const int n, const T rval, T* lvalAndResult);

	/**
	* The sum of two vectors in the first vector
	* Example: Sum(n,x,y)  means \f$ x=x+y \f$;
	*/	
	template <class T>
	void Sum(const int n, T* lvalAndResult, const T* rval);

	/**
	* The sum of two equal vectors
	* Example: Sum(n,x)  means \f$ x=x+x \f$;
	*/	
	template <class T>
	void Sum(const int n, T* lvalRvalAndResult);

	/**
	* The difference between two vectors in a third vector with control
	* Example: Difference(n,x,y,z)  means \f$ z=x-y \f$;
	*/	
	template <class T>
	void Difference(const int n, T* lval, T* rval, T* result);

	/**
	* The difference between two vectors in the first vector
	* Example: Difference(n,x,y)  means \f$ x=x-y \f$;
	*/	
	template <class T>
	void Difference(const int n, T* lvalAndResult, T* rval);

	/**
	* The difference between two vectors in the first vector
	* Example: DifferenceBis(n,x,y)  means \f$ y=x-y \f$;
	*/	
	template <class T>
	void DifferenceBis(const int n, T* lval, T* rvalAndResult);

	/**
	* The scalar (dot) product of two vectors	
	* Meaning: \f$ todo \f$;
	*/	
	template <class T>
	T Dot(const int n, const T* lval, const T* rval);

	/**
	* The scalar (dot) product of two vectors	
	* Meaning: \f$ todo \f$;
	*/
	template <class T>
	T UDot(const int n, const T* lval, const T* rval);

	/**
	* The product of a value with a vector	
	* Example: Product(n,c,v,r)	multiplies n values of the vector v by c, resulting in r
	*/	
	template <class T>
	void Prod(const int n, const T lval, const T* rval, T* result);

	/**
	* The product of a value with a vector	
	* Example: Product(n,c,v) multiplies n values of the vector v by c, substituting in v
	*/	
	template <class T>
	void Prod(const int n, const T lval, T* rvalAndResult);

	/**
	* Dividing a vector by a value
	* Example: Div(n,v,c,r) divides n elements of the vector v by c, resulting in r
	*/	
	template <class T>
	void Div(const int n, T* lval, T rval, T* result);

	/**
	* Dividing a vector by a value
	* Example: Div(n,v,c) divides n terms of the vector v by c, substituting in v
	*/	
	template <class T>
	void Div(const int n, T* lvalAndResult, T rval);

	template <class T>
	void ElementByElementProduct(const int n, const T* lval, const T* rval, T* result);

	template <class T>
	void ElementByElementDivision(const int n, const T* lval, const T* rval, T* result);

	template<typename T>
	void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, T& value);

	template<typename T>
	void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, T& value);

	template<typename T>
	void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, const int n, T* values);

	template<typename T>
	void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, const int n, T* values);

	template <typename T>
	inline unsigned char Compare_LE(const T a, const T b);

	template <typename T>
	inline unsigned char Compare_LT(const T a, const T b);

	template <typename T>
	inline unsigned char Compare_GE(const T a, const T b);

	template <typename T>
	inline unsigned char Compare_GT(const T a, const T b);

	template <typename T>
	void Sort(const int n, T *x, int *iS);

	template <typename T>
	std::vector<size_t> SortAndTrackIndicesIncreasing(std::vector<T> const& values);

	template <typename T>
	std::vector<size_t> SortAndTrackIndicesDecreasing(std::vector<T> const& values);

	template <typename T>
	void Load(T* value, std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat);

	template <typename T>
	void Save(const double value, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat);

	template<typename T>
	void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat);

	template<typename T>
	void CheckIfFileIsOpen(const T& fOut, const std::string file_name);

	/**
	*@brief Returns the estimation of the error (see Buzzi-Ferraris, formula 29.190): E = ||w*e|| / sqrt(ne)
	*@param w vector for whcih the error has to be estimated (typically z)
	*@param e the current vector of error weights: e(j) = 1. / ( tolAbs + tolRel*y(j) )
	*/
	template<typename T>
	double ErrorControl(const T& w, const T& e);

	/**
	*@brief Checks whether a value is present in a vector
	*@param T value to check
	*@param x vector to analyse
	*/
	template <typename T>
	bool IsValuePresent (const T value, const std::vector<T>& x);

	/**
	*@brief Returns the index of an element within a vector
	*@param T value to check
	*@param x vector to analyse
	*/
	template <typename T>
	int Index(const T value, const std::vector<T>& value_list);

	/**
	*@brief Transforms a variable into a string
	*@param T value to add
	*/
	template <typename T>
	std::string ToString(const T value);

	template<typename T>
	void Load(std::vector<T>& v, std::istream& fInput, const OpenSMOKE_File_Format fileFormat);

	template<typename T>
	void Exp(const std::vector<T>& lval, std::vector<T>* rval);
}

#include "OpenSMOKEUtilities.hpp"

#endif	// OpenSMOKE_OpenSMOKEUtilities_Hpp


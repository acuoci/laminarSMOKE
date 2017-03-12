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

#ifndef OpenSMOKE_OpenSMOKEUtilities_Cpp
#define OpenSMOKE_OpenSMOKEUtilities_Cpp

#include "OpenSMOKEStdInclude.h"
#include <numeric>
#include "OpenSMOKEFunctions.h"
#include <boost/math/constants/constants.hpp> 
#include<boost/tokenizer.hpp>

#if __cplusplus <= 199711L
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

//#include "Crypto-SHA256.h"

namespace OpenSMOKE
{
	template <typename T>
	void Swap(T* x, T *y)
	{
		T temp = *x;
		*x = *y;
		*y = temp;
	}


	template <typename T>
	inline T Abs(T const& a)
	{
	   return a < 0 ? -a : a;
	}


	template <typename T>
	inline T const& Max(T const& a, T const& b)
	{
		return (a > b ? a : b);
	}


	template <typename T>
	inline T const& Min(T const& a, T const& b)
	{
		return (a < b ? a : b);
	}


	template <typename T>
	inline T MaxAbs(T const& a, T const& b)
	{
		return ( Abs(a) > Abs(b) ? Abs(a) : Abs(b));
	}

	template <typename T>
	inline T MinAbs(T const& a, T const& b)
	{
		return ( Abs(a) < Abs(b) ? Abs(a) : Abs(b) );
	}

	template <typename T>
	inline T Max(const int n, T *x, int *im)
	{
		if(n < 0) 
			return x[0];
		
		T temp;
		temp = x[0];
		if(im != 0) 
			*im = 0;
		
		for(int i=1;i<n;i++)
			if(temp<x[i])
			{
				temp = x[i];
				if(im != 0)
					*im = i;
			}
		
		return temp;
	}

	template <class T>
	inline T Max(const int n, T *x)
	{
		if(n < 0) 
			return x[0];
		
		T temp;
		temp = x[0];
		for(int i=1;i<n;i++)
			if(temp < x[i])
				temp = x[i];
		
		return temp;
	}

	template <class T>
	inline T MaxAbs(const int n, T *x,int *im)
	{
		T temp = Abs(x[0]);
		if(n < 0) 
			return temp;
		
		if(im != 0) 
			*im = 0;
		
		for(int i=1;i<n;i++)
			if(temp < Abs(x[i]))
			{
				temp = Abs(x[i]);
				if(im != 0)
				*im = i;
			}

		return temp;
	}

	template <class T>
	inline T MaxAbs(const int n, T *x)
	{
		T temp = Abs(x[0]);
		if(n < 0) 
			return temp;
		
		for(int i=1;i<n;i++)
			if(temp < Abs(x[i]))
				temp = Abs(x[i]);

		return temp;
	}

	template <class T>
	inline T Min(const int n, T *x, int *im)
	{
		if(n < 0) 
			return x[0];
		
		T temp = x[0];
		if(im != 0) 
			*im = 0;
		
		for(int i = 1;i < n;i++)
			if(temp > x[i])
			{
				temp = x[i];
				if(im != 0)
					*im = i;
			}

		return temp;
	}

	template <class T>
	inline T Min(const int n, T *x)
	{
		if(n < 0)
			return x[0];
	
		T temp = x[0];
		for(int i=1;i<n;i++)
			if(temp > x[i])
				temp = x[i];

		return temp;
	}

	template <class T>
	inline T MinAbs(const int n, T *x, int *im)
	{
		T temp = Abs(x[0]);
		if(n < 0)
			return temp;
		if(im != 0) 
			*im = 0;
		for(int i = 1;i < n;i++)
			if(temp > Abs(x[i]))
			{
				temp = Abs(x[i]);
				if(im != 0)
				*im = i;
			}
		return temp;
	}

	template <class T>
	inline T MinAbs(const int n, T *x)
	{
		T temp = Abs(x[0]);
		
		if(n < 0) 
			return temp;
		
		for(int i = 1;i < n;i++)
			if(temp > Abs(x[i]))
				temp = Abs(x[i]);
		
		return temp;
	}

	template <class T>
	void Sum(const int n, T* lval, T* rval, T* result)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lval++) + (*rval++);
			*result++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, T* lval, const T rval, T* result)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lval++) + rval;
			*result++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, const T rval, T* lvalAndResult)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lvalAndResult) + rval;
			*lvalAndResult++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, T* lvalAndResult, T* rval)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lvalAndResult) + (*rval++);
			*lvalAndResult++ = sum;
		}
	}

	template <class T>
	void Sum(const int n, T* lvalRvalAndResult)
	{
		T sum;
		for(int i=0;i<n;i++)
		{
			sum = (*lvalRvalAndResult) + (*lvalRvalAndResult);
			*lvalRvalAndResult++ = sum;
		}
	}

	template <class T>
	void Difference(const int n, T* lval, T* rval, T* result)
	{
		T diff;
		for(int i=0;i<n;i++)
		{
			diff = (*lval++) - (*rval++);
			*result++ = diff;
		}
	}

	template <class T>
	void Difference(const int n, T* lvalAndResult, T* rval)
	{
		T diff;
		for(int i=0;i<n;i++)
		{
			diff = (*lvalAndResult) - (*rval++);
			*lvalAndResult++ = diff;
		}
	}

	template <class T>
	void DifferenceBis(const int n, T* lval, T* rvalAndResult)
	{
		T diff;
		for(int i=0;i<n;i++)
		{
			diff = (*lval++) - (*rvalAndResult);
			*rvalAndResult++ = diff;
		}
	}

	template <class T>
	T Dot(const int n, T* lval, T* rval)
	{
		T result = 0;
		for(int i=0;i<n;i++)
			result += (*lval++) * (*rval++);
		return result;
	}

	template <class T>
	T UDot(const int n, T* lval, T* rval)
	{
		T result = 0;
		for(int i=0;i<n;i++)
			result += (*lval++) / (*rval++);
		return result;
	}

	template <class T>
	void Prod(const int n, T lval, T* rval, T* result)
	{
		for(int i=0;i<n;i++)
			*result++ = lval * (*rval++);
	}

	template <class T>
	void Prod(const int n, T lval, T* rvalAndResult)
	{
		for(int i=0;i<n;i++)
		{
			*rvalAndResult = (lval) * (*rvalAndResult);
			 rvalAndResult++;
		}
	}

	template <class T>
	void Div(const int n, T* lval, T rval, T* result)
	{
		if(rval == 0)
			ErrorMessage("Div(const int n, T* lval, T rval, T* result)", "Division with 0");
		
		for(int i=0;i<n;i++)
			*result++ = (*lval++)/double(rval);
	}

	template <class T>
	void Div(const int n, T* lvalAndResult, T rval)
	{
		if(rval == 0)
			ErrorMessage("Div(const int n, T* lvalAndResult, T rval)", "Division with 0");
		for(int i=0;i<n;i++)
		{
			*lvalAndResult = (*lvalAndResult)/double(rval);
			 lvalAndResult++;
		}
	}

	template<typename T>
	void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, T& value)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fInput.read(reinterpret_cast < char * > (&value), sizeof(T)))
				ErrorMessage(	"Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, T& value)", 
								"I was unable to read an std::string coefficient from binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fInput >> value;
		}
	}

	template<typename T>
	void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, T& value)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fOutput.write( reinterpret_cast < char * > (&value), sizeof(T)))
				ErrorMessage(	"Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, T const& value)", 
								"I was unable to write on a binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << value << std::endl;
		}
	}

	template<typename T>
	void Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, const int n, T* values)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			for(int i=0;i<n;i++)
				if(!fInput.read(reinterpret_cast < char * > (&values[i]), sizeof(T)))
					ErrorMessage(	"Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat, const int n, T* values)", 
									"I was unable to read from binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			for(int i=0;i<n;i++)
				fInput >> values[i];
		}
	}

	template<typename T>
	void Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, const int n, T* values)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			for(int i=0;i<n;i++)
				if(!fOutput.write( reinterpret_cast < char * > (&values[i]), sizeof(T)))
					ErrorMessage(	"Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat, const int n, const T* values)", 
									"I was unable to write on a binary file");
		}
		else if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			for(int i=0;i<n;i++)
				fOutput << values[i] << std::endl;
		}
	}

	template <typename T>
	inline unsigned char Compare_LE(const T a, const T b)
	{
		return (a <= b);
	}

	template <typename T>
	inline unsigned char Compare_LT(const T a, const T b)
	{
		return (a < b);
	} 

	template <typename T>
	inline unsigned char Compare_GE(const T a, const T b)
	{
		return (a >= b);
	}

	template <typename T>
	inline unsigned char Compare_GT(const T a, const T b)
	{
		return (a > b);
	}

	template <typename T>
	void Sort(const int n, T *x, int *iS)
	{
		if(n <= 1)
			return;
	   
		int j,k,ik,jk;
		
		for(int node=1;node<n;node++)
		{
			int i = node;
			j = ((i + 1)/2) - 1;
			while(i != 0 && Compare_LE(x[j],x[i]))
			{
				Swap(x + j,x + i);
				Swap(iS + j,iS + i);
				i = j;
				j = ((i + 1)/2) - 1;
			}
		}

		for(int i=n-1;i>=1;i--)
		{
			Swap(x + i,x);
			Swap(iS + i,iS);
			k = i - 1;
			ik = 0;
			jk = 1;

			if(k >= 2 && Compare_GT(x[2],x[1]))
				jk = 2;

			while(jk <= k && Compare_GT(x[jk],x[ik]))
			{
				Swap(iS + jk,iS + ik);
				Swap(x + jk,x + ik);
				ik = jk;
				jk = (2*(ik + 1)) - 1;
				if(jk + 1 <= k)
					if(Compare_GT(x[jk + 1],x[jk]))
						jk++;
			}
		}
	}

	template <typename T>
	std::vector<size_t> SortAndTrackIndicesDecreasing(std::vector<T> const& values)
	{
		std::vector<size_t> indices(values.size());

		#if defined(_WIN32) || defined(_WIN64) 

		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort(begin(indices), end(indices), [&](size_t a, size_t b) { return values[a] > values[b]; });

        #elif __APPLE__

			ErrorMessage("sort_and_track_indices_decreasing", "sort_and_track_indices_decreasing not yet available for MacOSX");

        #else

			#if __cplusplus > 199711L

				struct IdxCompare
				{
					const std::vector<T>& target;
					IdxCompare(const std::vector<T>& target) : target(target) {}
					bool operator()(size_t a, size_t b) const { return target[a] > target[b]; }
				};

				for (size_t i = 0; i < indices.size(); ++i)
					indices[i] = i;

				std::sort(indices.begin(), indices.end(), IdxCompare(values));

			#else

				using namespace boost::phoenix;
				using namespace boost::phoenix::arg_names;

				int i = 0;
				std::transform(values.begin(), values.end(), indices.begin(), ref(i)++);
				std::sort(indices.begin(), indices.end(), ref(values)[arg1] > ref(values)[arg2]);

			#endif

		#endif

		return indices;
	}

	template <typename T>
	std::vector<size_t> SortAndTrackIndicesIncreasing(std::vector<T> const& values)
	{
		std::vector<size_t> indices(values.size());

		#if defined(_WIN32) || defined(_WIN64) 

		std::iota(begin(indices), end(indices), static_cast<size_t>(0));
		std::sort(begin(indices), end(indices), [&](size_t a, size_t b) { return values[a] < values[b]; });

		#elif __APPLE__
		
			ErrorMessage("sort_and_track_indices_increasing", "sort_and_track_indices_increasing not yet available for MacOSX");
	
		#else

			#if __cplusplus > 199711L

				struct IdxCompare
				{
					const std::vector<T>& target;
					IdxCompare(const std::vector<T>& target) : target(target) {}
					bool operator()(size_t a, size_t b) const { return target[a] < target[b]; }
				};

				for (size_t i = 0; i < indices.size(); ++i)
					indices[i] = i;

				std::sort(indices.begin(), indices.end(), IdxCompare(values));

			#else

				using namespace boost::phoenix;
				using namespace boost::phoenix::arg_names;

				int i = 0;
				std::transform(values.begin(), values.end(), indices.begin(), ref(i)++);
				std::sort(indices.begin(), indices.end(), ref(values)[arg1] < ref(values)[arg2]);

			#endif

		#endif

		return indices;
	}

	template <typename T>
	void Save(T value, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
				fOutput << value << std::endl;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fOutput.write( reinterpret_cast < char * > (&value), sizeof(T)))
				ErrorMessage("Save(T value, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)", "I was unable to write on binary file");
		}
	}

	template <typename T>
	void Load(T* value, std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
				fInput >> *value;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			if(!fInput.read(reinterpret_cast<char *>(value), sizeof(T)))
				ErrorMessage("Load(T* value, std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat)", "I was unable to read from binary file");
		}
	}

	template<typename T>
	void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_FORMATTED_FILE)
		{
			fOutput << vector.size() << std::endl;
			for(unsigned int i=0;i<vector.size();i++)
				fOutput << vector[i] << std::endl;
		}
		else if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			// Writing vector size
			if(!fOutput.write( reinterpret_cast < char * > (&vector.size()), sizeof(vector.size())))
				ErrorMessage("void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)", "I was unable to write on binary file");

			// Writing vector elements
			for(unsigned int i=0;i<vector.size();i++)
			if(!fOutput.write( reinterpret_cast < char * > (&vector[i]), sizeof(T)))
				ErrorMessage("void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)", "I was unable to write on binary file");
		}
	}

	template<typename T>
	void CheckIfFileIsOpen(const T& fOut, const std::string file_name)
	{
		if (!fOut.is_open())
		{
			std::cout << " ! Fatal error: I was unable to open this file: " << std::endl;
			std::cout << "   " << file_name << std::endl;
			std::cout << "   Press enter to exit..." << std::endl;
			getchar();
			exit(OPENSMOKE_FATAL_ERROR_EXIT);
		}
	}
        
	template<typename T>
	double ErrorControl(const T& w, const T& e)
	{
		const std::size_t ne = w.size();

		double sum = 0.;
		for (unsigned int i = 0; i < ne; i++)
		{
			const double coeff = w(i)*e(i);
			sum += coeff*coeff;
		}

		return std::sqrt( sum/double(ne) );
	}

	template <typename T>
	bool IsValuePresent (const T value, const std::vector<T>& x)
	{
		bool isPresent = false;

		for (unsigned int i = 0; i < x.size(); i++)
		if (x[i] == value)
		{
			isPresent = true;
			break;
		}

		return isPresent;
	}

	template <typename T>
	int Index(const T value, const std::vector<T>& value_list)
	{
		int index = -1;
		bool iFound = false;
		for (unsigned int i = 0; i < value_list.size (); i++)
		if (value == value_list[i])
		{
			index = i;
			iFound = true;
			break;
		}

		if (iFound == false)
		{
			std::stringstream err_str;
			err_str << "Value/String " << value << " not found in the corresponding list!" << std::endl;
			OpenSMOKE::FatalErrorMessage (err_str.str ());
		}

		return index;
	}

	template <typename T>
	std::string ToString(const T value)
	{
		std::ostringstream strs;
		strs << value;
		std::string str = strs.str ();
		return str;
	}
}

#endif	// OpenSMOKE_OpenSMOKEUtilities_Cpp

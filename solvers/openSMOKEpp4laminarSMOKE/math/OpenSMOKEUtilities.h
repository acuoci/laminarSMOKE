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
	void Sum(const int n, T* lval, T* rval, T* result);

	/**
	* The sum of a constant to a vector
	* Example: Sum(n,x,c,z) means \f$ z=x+c \f$;
	*/	
	template <class T>
	void Sum(const int n, T* lval, const T rval, T* result);

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
	void Sum(const int n, T* lvalAndResult, T* rval);

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
	T Dot(const int n, T* lval, T* rval);

	/**
	* The scalar (dot) product of two vectors	
	* Meaning: \f$ todo \f$;
	*/
	template <class T>
	T UDot(const int n, T* lval, T* rval);

	/**
	* The product of a value with a vector	
	* Example: Product(n,c,v,r)	multiplies n values of the vector v by c, resulting in r
	*/	
	template <class T>
	void Prod(const int n, T lval, T* rval, T* result);

	/**
	* The product of a value with a vector	
	* Example: Product(n,c,v) multiplies n values of the vector v by c, substituting in v
	*/	
	template <class T>
	void Prod(const int n, T lval, T* rvalAndResult);

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
	std::vector<size_t> sort_and_track_indices_increasing(std::vector<T> const& values);

	template <typename T>
	std::vector<size_t> sort_and_track_indices_decreasing(std::vector<T> const& values);

	template <typename T>
	void Load(T* value, std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat);

	template <typename T>
	void Save(const double value, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat);

	template<typename T>
	void Save(const std::vector<T>& vector, std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat);

	template<typename T>
	void CheckIfFileIsOpen(const T& fOut, const std::string file_name);

	std::string GetCurrentTime();
	std::string GetCurrentDate();
	void OpenSMOKE_logo(const std::string application_name, const std::string author_name);
	std::string SplitStringIntoSeveralLines(std::string source, const std::size_t width = 60, const std::string whitespace = " \t\r");
	void PrintTagOnASCIILabel(const unsigned int width, std::ostream &fOut, const std::string tag, unsigned int &counter);
	void SetXMLFile(std::ostream& xml_string);
	void CreateDirectory(const boost::filesystem::path& path_output);
	void OpenSMOKE_CheckLicense(const boost::filesystem::path executable_folder, const std::string compilation_date_os, const long maximum_license_days);
	boost::filesystem::path GetExecutableFileName(char** argv);
        
        void OpenInputFileXML(rapidxml::xml_document<>& doc, std::vector<char>& xml_copy, const boost::filesystem::path& file_name);
        
	void OpenInputFileASCII(std::ifstream &fASCII, const boost::filesystem::path input_file_ascii);

	

        /**
        *@brief Open the output files and checks any error
        *@param fXML the output stream used
        *@param output_file_xml XML file where the output will be written 
        */
        void OpenOutputFileXML(std::ofstream &fXML, const boost::filesystem::path output_file_xml);
    
        /**
        *@brief Opens the output files and checks any error
        *@param fASCII the output stream used
        *@param output_file_ascii ASCII file where the output will be written
        */
        void OpenOutputFileASCII(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii);
        
		/**
        *@brief Opens an existing output files, checks any error and continue writing on it
        *@param fASCII the output stream used
        *@param output_file_ascii ASCII file where the output will be written
        */
	void OpenOutputFileASCII_Append(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii);


        void CheckKineticsFolder(const boost::filesystem::path& path_output);

	/**
	*@brief Returns the estimation of the error (see Buzzi-Ferraris, formula 29.190): E = ||w*e|| / sqrt(ne)
	*@param w vector for whcih the error has to be estimated (typically z)
	*@param e the current vector of error weights: e(j) = 1. / ( tolAbs + tolRel*y(j) )
	*/
	template<typename T>
	double ErrorControl(const T& w, const T& e);

	/**
	*@brief Returns the sparsity pattern of a tridiagonal block matrix
	*@param number_blocks number of blocks
	*@param block_size size of a single block
	*@param rows row indices of non zero elements (1 index based)
	*@param cols column indices of non zero elements (1 index based)
	*/

	void SparsityPatternTridiagonal(const unsigned int number_equations, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);
	void SparsityPatternPentadiagonal(const unsigned int number_equations, const unsigned int width, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

	void SparsityPatternBlock(const unsigned int number_blocks, const unsigned int block_size,
		const std::vector<unsigned int>& rows_single, const std::vector<unsigned int>& cols_single,
		std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

	/**
	*@brief Returns the number of lines of a text file
	*@param path path of the file
	*/
	unsigned int filelines(const boost::filesystem::path path);

	/**
	*@brief Returns the number of words in a stringstream
	*@param str stringstream to analyze
	*/
	int countwords (std::stringstream& str);

	/**
	*@brief Returns the real solution of a cubic equation
	*@param a third order coefficient
	*@param b second order coefficient
	*@param c first order coefficient
	*@param d zero order coefficient
	*/
	std::vector<double> cubicrootsreal(const double a, const double b, const double c, const double d);
  
	/**
	*@brief Returns a normalized vector (sum of elements = 1)
	*@param x vector to normalize
	*/
	std::vector<double> normalize(std::vector<double>& x);

	/**
	*@brief Checks whether a value is present in a vector
	*@param T value to check
	*@param x vector to analyse
	*/
	template <typename T>
	bool isValuePresent (const T value, const std::vector<T>& x);

	/**
	*@brief Evaluates the median value of a vector<double>
	*@param vec vector to analyse
	*/
	double median (const std::vector<double> vec);

	/**
	*@brief Evaluates the median absolute deviation of a vector<double>
	*@param vec vector to analyse
	*/
	double mad (const std::vector<double> vec);

	/**
	*@brief Returns the index of an element within a vector
	*@param T value to check
	*@param x vector to analyse
	*/
	template <typename T>
	int index (const T value, const std::vector<T>& value_list);

	/**
	*@brief Transforms a variable into a string
	*@param T value to add
	*/
	template <typename T>
	std::string ToString(const T value);

	/**
	*@brief Checks whether sum of elements of a vector<double> is equal to 1, without elements < 0 or > 1
	*@param x vector to analyse
	*/
	void CheckSumOfFractions(std::vector<double>& x);

	/**
	*@brief Converts temperature into [K]
	*@param temperature temperature to convert
	*@param unit original unit of measure
	*/
	double SI_temperature (const double temperature, const std::string unit);

	/**
	*@brief Converts pressure into [Pa]
	*@param pressure pressure to convert
	*@param unit original unit of measure
	*/
	double SI_pressure (const double pressure, const std::string unit);

	/**
	*@brief Converts length into [m]
	*@param length length to convert
	*@param unit original unit of measure
	*/
	double SI_length (const double length, const std::string unit);

	/**
	*@brief Converts volume into [m3]
	*@param volume volume to convert
	*@param unit original unit of measure
	*/
	double SI_volume (const double volume, const std::string unit);

	/**
	*@brief Converts velocity into [m/s]
	*@param velocity velocity to convert
	*@param unit original unit of measure
	*/
	double SI_velocity (const double velocity, const std::string unit);

	/**
	*@brief Converts time into [s]
	*@param time time to convert
	*@param unit original unit of measure
	*/
	double SI_time (const double time, const std::string unit);
}

#include "OpenSMOKEUtilities.hpp"

#endif	// OpenSMOKE_OpenSMOKEUtilities_Hpp


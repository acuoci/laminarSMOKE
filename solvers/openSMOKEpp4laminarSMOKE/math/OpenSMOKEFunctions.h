/***************************************************************************
 *   Copyright (C) 2011 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef OpenSMOKE_OpenSMOKEFunctions_Hpp
#define OpenSMOKE_OpenSMOKEFunctions_Hpp

#include "OpenSMOKEStdInclude.h"

namespace OpenSMOKE
{
	/**
	* Utility to print fatal error messages
	*/
	void ErrorMessage(const std::string functionName, const std::string errorMessage);

	/**
	* Utility to print fatal error messages
	*/
	int FatalErrorMessage(const std::string errorMessage);

	/**
	* Returns the MachEps (float)
	*/
	float MachEpsFloat();
	
	/**
	* Returns the MachEps (double)
	*/
	double MachEps();

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	float SqrtSumSqr(const int n, float *x);

	/**
	* Returns the square root of the sum of the squares of x
	* Example: z=SqrtSumSqr(n,&x) means \f$ z=todo \f$
	*/
	double SqrtSumSqr(const int n, double *x);
    
    /**
	* Returns the elapsed time in seconds since the process started
	*/
    double OpenSMOKEClock(void);

    /**
	* Returns the elapsed time in seconds since the OpenSMOKEDiffClock() started
	*/
	double OpenSMOKEGetCpuTime(void);

    /**
	* Returns the number of cpu clocks
	*/
	#if defined(_WIN32) || defined(_WIN64) 
		unsigned __int64 OpenSMOKEGetCpuClocks(void);
	#else
		unsigned long int OpenSMOKEGetCpuClocks(void);
	#endif

	/**
	* Returns the cpu frequency
	*/
    double OpenSMOKEGetCpuFrequency(void);

    /**
	* Returns the max cpu frequency
	*/
    double OpenSMOKEGetMaxCpuFrequency(void);

    /**
	* Returns the max cpu clock frequency
	*/
    double OpenSMOKEGetCpuClocksFrequency(void);

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

	/**
	*@brief Checks if the folder containing the preprocessed kinetic mechanism exists
	*@param path_output the path to the kinetic folder
	*/
	void CheckKineticsFolder(const boost::filesystem::path& path_output);

	/**
	*@brief Open a XML file
	*@param doc the XML class
	*@param xml_copy the content of XML file
	*@param file_name the path to the XML file
	*/
	void OpenInputFileXML(rapidxml::xml_document<>& doc, std::vector<char>& xml_copy, const boost::filesystem::path& file_name);

	/**
	*@brief Open a ASCII file
	*@param fASCII the file stream
	*@param input_file_ascii the path to the ASCII file
	*/
	void OpenInputFileASCII(std::ifstream &fASCII, const boost::filesystem::path input_file_ascii);

	/**
	*@brief Writes the header lines on a XML file
	*@param xml_streamthe XML stream
	*/
	void SetXMLFile(std::ostream& xml_stream);

	/**
	*@brief Writes the tag (label+index of column) in ASCII files
	*@param width width of the column
	*@param fOut output file stream
	*@param tag label to be written
	*@param counter column index to be written
	*/
	void PrintTagOnASCIILabel(const unsigned int width, std::ostream &fOut, const std::string tag, unsigned int &counter);

	/**
	*@brief Creates a new directory
	*@param path_output path to the new directory to be created
	*/
	void CreateDirectory(const boost::filesystem::path& path_output);

	/**
	*@brief Returns the current time in readable format
	*/
	std::string GetCurrentTime();

	/**
	*@brief Returns the current time in raw format
	*/
	std::string GetRawCurrentTime();
	
	/**
	*@brief Returns the current date in readable format
	*/
	std::string GetCurrentDate();

	/**
	*@brief Splits a given string into several lines
	*@param source string to be split
	*@param width maximum length of each line
	*@param whitespace chars to be considered targets for splitting the string
	*/
	std::string SplitStringIntoSeveralLines(std::string source, const std::size_t width = 60, const std::string whitespace = " \t\r");

	/**
	*@brief Get the complete name of an exectutable file
	*@param argv argv[0] is the name (only) of the executable 
	*/
	boost::filesystem::path GetExecutableFileName(char** argv);

	/**
	*@brief Writes the OpenSMOKE++ logo on the screen
	*@param application_name name of the OpenSMOKE++ solver
	*@param author_name author's name
	*/
	void OpenSMOKE_logo(const std::string application_name, const std::string author_name);

	/**
	*@brief Writes the OpenSMOKE++ logo on a file
	*@param fOut output stream
	*@param application_name name of the OpenSMOKE++ solver
	*/
	void OpenSMOKE_logo(std::ofstream& fOut, const std::string application_name);

	/**
	*@brief Calculates the width of a column related to a chemical species
	*@param species_name species' name
	*@param max_number_species maximum number of species
	*/
	unsigned int CalculateSpeciesFieldWidth(const std::string species_name, const int max_number_species);

	/**
	*@brief Imports the reaction strings from a XML file
	*@param file_name XML file name
	*@param n total number of reactions
	*@param reaction_strings vector containing the reaction strings
	*/
	void ImportReactionNames(const boost::filesystem::path file_name, const unsigned int n, std::vector<std::string>& reaction_strings);

	/**
	*@brief Calculates the real roots of a cubic function a*x^3+b*x^2+c*x+d=0
	*@param a a coefficient
	*@param b b coefficient
	*@param c c coefficient
	*@param d d coefficient
	*/
	std::vector<double> CubicRootsReal(const double a, const double b, const double c, const double d);

	/**
	*@brief Given an input vector, returns a normalized vector
	*@param x original vector
	*/
	std::vector<double> Normalize(const std::vector<double>& x);

	/**
	*@brief Returns the median of a given vector
	*@param v the given vector
	*/
	double Median(const std::vector<double>& v);

	/**
	*@brief Returns the median absolute deviation of a given vector
	*@param v the given vector
	*/
	double MedianAbsoluteDeviation(const std::vector<double>& v);

	/**
	*@brief Checks and corrects the elements of a vector whose sum of elements is supposed to be 1.
	*@param x the given vector
	*@param boundary_eps maximum error with respect to 0. or 1. for each element
	*@param eps maximum error on the sum
	*/
	void CheckAndCorrectSumOfFractions(std::vector<double>& x, const double boundary_eps = 1.e-8, const double eps = 1.e-4);
	
	/**
	*@brief Returns the number of lines of a text file
	*@param path path of the file
	*/
	unsigned int NumberOfLinesInFile(const boost::filesystem::path& path);

	/**
	*@brief Returns the number of words in a stringstream
	*@param str stringstream to analyze
	*/
	unsigned int CountWordsInString(std::stringstream& str);

	/**
	*@brief Returns the sparsity pattern of a tridiagonal matrix
	*@param number_equations total number of equations
	*@param rows row indices of non zero elements (1 index based)
	*@param cols column indices of non zero elements (1 index based)
	*/
	void SparsityPatternTridiagonal(const unsigned int number_equations, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);
	
	/**
	*@brief Returns the sparsity pattern of a pentadiagonal matrix
	*@param number_equations total number of equations
	*@param width distance between the main diagonal and the external diagonals
	*@param rows row indices of non zero elements (1 index based)
	*@param cols column indices of non zero elements (1 index based)
	*/
	void SparsityPatternPentadiagonal(const unsigned int number_equations, const unsigned int width, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

	/**
	*@brief Returns the sparsity pattern of a block matrix
	*@param number_blocks number of blocks
	*@param block_size size of a single block
	*@param rows_single rows indices of non-zero blocks (as if they were 1x1)
	*@param cols_single columns indices of non-zero blocks (as if they were 1x1)
	*@param rows row indices of non zero elements (1 index based)
	*@param cols column indices of non zero elements (1 index based)
	*/
	void SparsityPatternBlock(	const unsigned int number_blocks, const unsigned int block_size,
								const std::vector<unsigned int>& rows_single, const std::vector<unsigned int>& cols_single,
								std::vector<unsigned int>& rows, std::vector<unsigned int>& cols);

}

#include "OpenSMOKEFunctions.hpp"

#endif	// OpenSMOKE_OpenSMOKEFunctions_Hpp


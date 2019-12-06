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

#ifndef OpenSMOKE_KineticsUtilityFunctions_H
#define	OpenSMOKE_KineticsUtilityFunctions_H

namespace OpenSMOKE_Utilities
{
	#if defined(_WIN32) || defined(_WIN64) 
		#define snprintf _snprintf 
		#define vsnprintf _vsnprintf 
		#define strcasecmp _stricmp 
		#define strncasecmp _strnicmp 
	#endif

	/**
	*@brief This is a struct for reordering a couple of vectors
	*/	
	template<class T1, class T2, class Pred = std::less<T1> >
	struct sort_pair_second 
	{
		bool operator()(const std::pair<T1,T2>&left, const std::pair<T1,T2>&right) 
		{
			Pred p;
			return p(left.first, right.first);
		}
	};

	typedef std::vector<double>::const_iterator myiterdouble;

	/**
	*@brief This is a struct for reordering a vector in ascending order
	*/	
	struct ordering 
	{
		bool operator ()(std::pair<size_t, myiterdouble> const& a, std::pair<size_t, myiterdouble> const& b) 
		{
			return *(a.second) < *(b.second);
		}
	};

	/**
	*@brief Modification of Manuel's answer
	*/
	/*
	struct ciLessBoost : std::binary_function<std::string, std::string, bool>
	{
		bool operator() (const std::string & s1, const std::string & s2) const 
		{
			return boost::algorithm::lexicographical_compare(s1, s2, boost::is_iless());
		}
	};
	*/

	struct ciLessBoost
	{
		bool operator() (const std::string & s1, const std::string & s2) const
		{
			return boost::algorithm::lexicographical_compare(s1, s2, boost::is_iless());
		}
	};

	/**
	* @brief recommended in Meyers, Effective STL when internationalization and embedded. 
	*        NULLs aren't an issue.  Much faster than the STL or Boost lex versions.
	*/
	/*
	struct ciLessLibC : public std::binary_function<std::string, std::string, bool> 
	{
		bool operator()(const std::string &lhs, const std::string &rhs) const {
			return strcasecmp(lhs.c_str(), rhs.c_str()) < 0 ? 1 : 0;
		}
	};
	*/

	struct ciLessLibC
	{
		bool operator()(const std::string &lhs, const std::string &rhs) const {
			return strcasecmp(lhs.c_str(), rhs.c_str()) < 0 ? 1 : 0;
		}
	};

	/**
	*@brief Reorders the first vector in ascending order; the second vector is reordered on the basis of the first reordering
	*/	
	template<class T1, class T2>
	bool ReorderPairsOfVectors(std::vector<T1>& v1, std::vector<T2>& v2);

	/**
	*@brief Compacts the first vector and the second one accordingly
	*/	
	template<class T1, class T2>
	bool CleanPairsOfVectors(std::vector<T1>& v1, std::vector<T2>& v2);

	/**
	*@brief Compares element-by-element two different vectors and returns true only if all the elements are equal
	*/	
	template<typename T>
	bool compare_vectors(const std::vector<T>& v1, const std::vector<T>& v2);

	/**
	*@brief Recognizes the three kinetic parameters in the first line of a reaction in CHEMKIN format
	*/	
	bool AnalyzeReactionLine(std::string const&  firstline, double &A, double &beta, double &E, std::string& reconstructedline);
	
	/**
	*@brief Separates the reactant and product sides in a reaction written in CHEMKIN format
	*/
	bool SeparateReactantSideFromProductSide(const std::string& reconstructedline, std::string &line_reactants, std::string &line_products, bool &iReversible);

	/**
	*@brief Separates the names of the species from the stoichiometric coefficients
	*/
	bool SeparateSpeciesAndStoichiometricCoefficients(std::string& line, std::vector<std::string>& species_names, std::vector<double> &coefficients);
	
	/**
	*@brief Separates the name of the species from the stoichiometric coefficient
	*/
	bool SeparateStoichiometricCoefficientsFromSpecies(const std::string &line, std::string &species, double &number);
	
	/**
	*@brief Recognize a key-word (according to the CHEMKIN syntax rules)
	*/	
	bool LookForKeyWord(const std::string& line, std::string& keyword);
	
	/**
	*@brief Reads coefficients
	*/		
	bool ReadCoefficients(const std::string tag, std::string& line, std::vector<double>& coefficients);
	
	/**
	*@brief Analyzes an option line with n numerical coefficients (in CHEMKIN format) 
	*/		
	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n, std::vector<double>& coefficients);

	/**
	*@brief Analyzes an option line with n1 or n2 numerical coefficients (in CHEMKIN format) 
	*/	
	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n1, const int n2, std::vector<double>& coefficients);

	/**
	*@brief Analyzes an option line with n strings (in CHEMKIN format) 
	*/
	bool ReadReactionKeyWordPlusWords(const std::string tag, std::string& line, const int n, std::vector<std::string>& coefficients);

	/**
	*@brief Analyzes an option line with n1 or n2 strings (in CHEMKIN format) 
	*/
	bool ReadReactionKeyWordPlusWords(const std::string tag, std::string& line, const int n1, const int n2, std::vector<std::string>& coefficients);
	
	/**
	*@brief Analyzes an option line with any number of strings (in CHEMKIN format) 
	*/
	bool ReadReactionKeyWordPlusWords(const std::string tag, std::string& line, std::vector<std::string>& words);

	/**
	*@brief TODO
	*/	
	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n, std::string& word, std::vector<double>& coefficients);

	/**
	*@brief TODO
	*/
	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n1, const int n2, std::string& word, std::vector<double>& coefficients);
	
	/**
	*@brief Returns the string corresponding to the reaction in a readable format
	*/	
	void GetKineticConstantString(const double A, const double beta, const double E, const double lambda, std::string& line_kForwardStringSI, std::string& line_kForwardStringCGS);
	
	/**
	*@brief Returns the units of the kinetic constants in a readable format (SI units)
	*/	
	std::string GetUnitsOfKineticConstantsSI(const double lambda);
	
	/**
	*@brief Returns the units of the kinetic constants in a readable format (CGS units)
	*/
	std::string GetUnitsOfKineticConstantsCGS(const double lambda);

}

#include "KineticsUtilityFunctions.hpp"


#endif	/* OpenSMOKE_KineticsUtilityFunctions_H */


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

#include <boost/regex.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>

namespace OpenSMOKE_Utilities
{
	template<typename T>
	bool compare_vectors(const std::vector<T>& v1, const std::vector<T>& v2)
	{
		if (v1.size() != v2.size())
			return false;
	
		for(unsigned int i=0;i<v1.size();i++)
			if (v1[i] != v2[i])
				return false;

		return true;
	}

	template <typename T>
	std::vector<T> sort_from_ref( std::vector<T> const& in, std::vector< std::pair< size_t, myiterdouble> > const& reference) 
	{
		std::vector<T> ret(in.size());

		size_t const size = in.size();
		for (size_t i = 0; i < size; ++i)
			ret[i] = in[reference[i].first];

		return ret;
	}

	// Reorder indices
	template<class T1, class T2>
	bool ReorderPairsOfVectors(std::vector<T1>& v1, std::vector<T2>& v2)
	{
		if (v1.size() != v2.size())
			return false;
	
		std::vector<std::pair<T1, T2> > vec(v1.size());
		for(unsigned int i=0;i<v1.size();i++)
		{
			vec[i].first  = v1[i];
			vec[i].second = v2[i];
		}
	
		std::sort(vec.begin(), vec.end(), sort_pair_second<T1, T2>());

		for(unsigned int i=0;i<v1.size();i++)
		{
			v1[i] = vec[i].first;
			v2[i] = vec[i].second;
		}

		return true;
	}

	template<class T1, class T2>
	bool CleanPairsOfVectors(std::vector<T1>& v1, std::vector<T2>& v2)
	{
		if (v1.size() != v2.size())
			return false;
	
		std::vector<std::pair<T1, T2> > vec(v1.size());
		for (std::size_t i = v1.size() - 1; i>0; i--)
		{
			if (v1[i] == v1[i-1])
			{
				v2[i-1] += v2[i];
				v1[i] = -1;
				v2[i] = -1;
			}
		}

		v1.resize(std::remove(v1.begin(), v1.end(), -1) - v1.begin());
		v2.resize(std::remove(v2.begin(), v2.end(), -1) - v2.begin());

		if (v1.size() != v2.size())
			return false;

		return true;
	}

	bool AnalyzeReactionLine(std::string const&  firstline, double &A, double &beta, double &E, std::string& reconstructedline)
	{
		// Separating reaction string from parameters
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_only_blanks;
		boost::char_separator<char> sep(" ");
		tokenizer_only_blanks tokens(firstline, sep);
		const std::size_t count = std::distance(tokens.begin(), tokens.end());
		if (count<4)
		{
			std::cout << "Error in defining the kinetic parameters: wrong number of arguments";
			return false;
		}

		try
		{
			tokenizer_only_blanks::iterator tok_iter = tokens.begin();
			std::advance(tok_iter, count-3);
			A = boost::lexical_cast<double>(*tok_iter);
			tok_iter++;
			beta = boost::lexical_cast<double>(*tok_iter);
			tok_iter++;
			E = boost::lexical_cast<double>(*tok_iter);
		}
		catch(boost::bad_lexical_cast &)
		{
			std::cout << "\nNumerical conversion failure at the first line.\nThe kinetic parameters are not written properly (they are not numbers)" << std::endl;
			return false;
		}  

		tokenizer_only_blanks::iterator tok_iter = tokens.begin();
		for(unsigned int i=1;i<=count-3;i++)
		{
			reconstructedline += *tok_iter;
			tok_iter++;
		}

		return true;
	}

	bool SeparateSpeciesAndStoichiometricCoefficients(std::string&  line, std::vector<std::string>& species_names, std::vector<double> &coefficients)
	{
		// Pressure dependent reactions
		boost::replace_all(line, "(+M)", "+$PDR$");
		boost::replace_all(line, "(+m)", "+$PDR$");
	
		// Third body species in pressure dependent reactions
		{
			std::string string_rule = "[\\(][\\+].*[\\)]";

			//boost::regex rule("[\(][\+].*[\)]", boost::regex_constants::icase);
			boost::regex rule(string_rule.c_str(), boost::regex_constants::icase);
			boost::smatch matches;

			if (boost::regex_search(line, matches, rule)) 
			{
				for (unsigned int i = 0; i < matches.size(); ++i)
				{
					std::string name_original = matches[i];
					std::string name = matches[i];
					boost::replace_all(name, "(+", "+$PDS$");
					name = name.substr(0,name.length()-1);

					boost::replace_all(line, name_original, name);
				}
			}
		}
	
		// Adding end
		line += "\n";
	
		// Ionic species
		boost::replace_all(line, "++",   "$PLUS$+");
		boost::replace_all(line, "+\n",  "$PLUS$");
	
		// Third body reactions
		boost::replace_all(line, "+M\n", "+$TBR$");
		boost::replace_all(line, "+m\n", "+$TBR$");

		// Removing end
		boost::replace_all(line, "\n", "");

		// Tokenize
		std::vector<std::string> names;
		{
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_plus;
			boost::char_separator<char> sep("+");
			tokenizer_plus tokens(line, sep);
			const std::size_t count = std::distance(tokens.begin(), tokens.end());
			names.resize(count);
			std::vector<std::string>::iterator names_iter=names.begin();
			for (tokenizer_plus::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
			{
				*names_iter = *tok_iter;
				boost::replace_all(*names_iter, "$PLUS$", "+");
				names_iter++;
			}  
		}


		coefficients.resize(names.size());
		species_names.resize(names.size());
		unsigned int i=0;
		for(std::vector<std::string>::iterator names_iter=names.begin();names_iter != names.end(); ++names_iter)
		{
			SeparateStoichiometricCoefficientsFromSpecies(*names_iter, species_names[i], coefficients[i]);
			i++;
		}

		return true;
	}

	bool SeparateStoichiometricCoefficientsFromSpecies(const std::string &line, std::string &species, double &number)
	{
		if (isdigit(line[0]) == false && line[0]!='.')
		{
			species = line;
			number = 1.;
			return true;
		}
		else
		{
			for(unsigned int j=0;j<line.size();j++)
			{
				if (isdigit(line[j]) == false && line[j]!='.')
				{
					species = line.substr(j, line.size()-j);
					std::string number_string = line.substr(0,j);
			
					try
					{
						number = boost::lexical_cast<double>(number_string);
					}
					catch(boost::bad_lexical_cast &)
					{
						std::cout << "Numerical conversion failure at the first line. The stoichiometric coefficients are not written properly (they are not numbers)" << std::endl;
						return false;
					}

					return true;
				}
			}

			return false;
		}
	}

	std::string GetUnitsOfKineticConstantsSI(const double lambda)
	{
		std::string units;

			 if (lambda == 1.)	units += "[1/s]";
		else if (lambda == 2.)	units += "[m3/kmol/s]";
		else if (lambda == 3.)	units += "[m6/kmol2/s]";
		else if (lambda == 4.)	units += "[m9/kmol3/s]";
		else
		{
			std::stringstream exponent;
			std::stringstream exponent3;
			exponent << lambda - 1.;
			exponent3 << 3.*(lambda - 1.);
			units += "[m" + exponent3.str() + "/kmol" + exponent.str() + "/s]";
		}

		return units;
	}

	std::string GetUnitsOfKineticConstantsCGS(const double lambda)
	{
		std::string units;

			 if (lambda == 1.)	units += "[1/s]";
		else if (lambda == 2.)	units += "[cm3/mol/s]";
		else if (lambda == 3.)	units += "[cm6/mol2/s]";
		else if (lambda == 4.)	units += "[cm9/mol3/s]";
		else
		{
			std::stringstream exponent;
			std::stringstream exponent3;
			exponent << lambda - 1.;
			exponent3 << 3.*(lambda - 1.);
			units += "[cm" + exponent3.str() + "/mol" + exponent.str() + "/s]";
		}

		return units;
	}

	void GetKineticConstantString(const double A, const double beta, const double E, const double lambda, std::string& line_kForwardStringSI, std::string& line_kForwardStringCGS)
	{
		line_kForwardStringSI  = "k = ";
		line_kForwardStringCGS = "k = ";

		// Frequency factor
		{
			std::stringstream ANumber_;
			ANumber_ << std::scientific << std::setprecision(6) << A;
			line_kForwardStringSI += ANumber_.str();
		}
		{
			std::stringstream ANumber_;
			ANumber_ << std::scientific << std::setprecision(6) << A/std::pow(1.e3, 1.-lambda);
			line_kForwardStringCGS += ANumber_.str();
		}
		
		// Temperature exponent
		if (beta != 0.) 
		{
			std::stringstream betaNumber_;
			betaNumber_ << std::fixed << std::setprecision(3) << beta;
			line_kForwardStringSI  += " T^" + betaNumber_.str();
			line_kForwardStringCGS += " T^" + betaNumber_.str();
		}

		// Activation energy
		if (E != 0.)    
		{
			std::stringstream ENumber_;
			ENumber_ << std::scientific << std::setprecision(6) << -E;
			line_kForwardStringSI += " exp(" + ENumber_.str() + "/RT)";
		}
		if (E != 0.)    
		{
			std::stringstream ENumber_;
			ENumber_ << std::fixed << std::setprecision(2) << -E/(1.e3*Conversions::J_from_cal);
			line_kForwardStringCGS += " exp(" + ENumber_.str() + "/RT)";
		}
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n1, const int n2, std::vector<double>& coefficients)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);
						
		if (std::distance (tokens_slash.begin(), tokens_slash.end()) != 2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (1)!" << std::endl;
			return false;
		}
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;
	
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());
		if ( n!=n1 && n!=n2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (2)!" << std::endl;
			return false;
		}
		try
		{
			coefficients.resize(n);
			tokenizer_blank::iterator tok_blank = tokens_blank.begin();
			for(std::size_t i=0;i<n;i++)
			{
				coefficients[i] = boost::lexical_cast<double>(*tok_blank);
				++tok_blank;
			}
		}
		catch(boost::bad_lexical_cast &)
		{
			std::cout << tag << "The coefficients are not written properly (they are not numbers)" << std::endl;
			return false;
		}  

		return true;
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n1, const int n2, std::vector<std::string>& coefficients)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);

		if (std::distance(tokens_slash.begin(), tokens_slash.end()) != 2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (1)!" << std::endl;
			return false;
		}
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());
		if (n != n1 && n != n2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (2)!" << std::endl;
			return false;
		}

		coefficients.resize(n);
		tokenizer_blank::iterator tok_blank = tokens_blank.begin();
		for (std::size_t i = 0; i<n; i++)
		{
			coefficients[i] = *tok_blank;
			++tok_blank;
		}

		return true;
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n1, const int n2, std::string& word, std::vector<double>& coefficients)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);
						
		if (std::distance (tokens_slash.begin(), tokens_slash.end()) != 2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (1)!" << std::endl;
			return false;
		}
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;
	
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());
		if ( n!=n1 && n!=n2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (2)!" << std::endl;
			return false;
		}
		try
		{
			coefficients.resize(n-1);
			tokenizer_blank::iterator tok_blank = tokens_blank.begin();
			word = *tok_blank;
			++tok_blank;
			for(std::size_t i=0;i<n-1;i++)
			{
				coefficients[i] = boost::lexical_cast<double>(*tok_blank);
				++tok_blank;
			}
		}
		catch(boost::bad_lexical_cast &)
		{
			std::cout << tag << "The coefficients are not written properly (they are not numbers)" << std::endl;
			return false;
		}  

		return true;
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n1, const int n2, std::string& word, std::vector<std::string>& coefficients)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);

		if (std::distance(tokens_slash.begin(), tokens_slash.end()) != 2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (1)!" << std::endl;
			return false;
		}
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());
		if (n != n1 && n != n2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (2)!" << std::endl;
			return false;
		}

		coefficients.resize(n - 1);
		tokenizer_blank::iterator tok_blank = tokens_blank.begin();
		word = *tok_blank;
		++tok_blank;
		for (std::size_t i = 0; i<n - 1; i++)
		{
			coefficients[i] = boost::lexical_cast<std::string>(*tok_blank);
			++tok_blank;
		}

		return true;
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n, std::vector<double>& coefficients)
	{
		return ReadReactionKeyWordPlusCoefficients(tag, line, n, n, coefficients);
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n, std::vector<std::string>& coefficients)
	{
		return ReadReactionKeyWordPlusCoefficients(tag, line, n, n, coefficients);
	}

	bool ReadReactionKeyWordPlusWords(const std::string tag, std::string& line, const int n, std::vector<std::string>& words)
	{
		return ReadReactionKeyWordPlusWords(tag, line, n, n, words);
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n, std::string& word, std::vector<double>& coefficients)
	{
		return ReadReactionKeyWordPlusCoefficients(tag, line, n, n, word, coefficients);
	}

	bool ReadReactionKeyWordPlusCoefficients(const std::string tag, std::string& line, const int n, std::string& word, std::vector<std::string>& coefficients)
	{
		return ReadReactionKeyWordPlusCoefficients(tag, line, n, n, word, coefficients);
	}

	bool ReadCoefficients(const std::string tag, std::string& line, std::vector<double>& coefficients)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);
						
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;
	
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());
		try
		{
			coefficients.resize(n);
			tokenizer_blank::iterator tok_blank = tokens_blank.begin();
			for (std::size_t i = 0; i<n; i++)
			{
				coefficients[i] = boost::lexical_cast<double>(*tok_blank);
				++tok_blank;
			}
		}
		catch(boost::bad_lexical_cast &)
		{
			std::cout << tag << "The coefficients are not written properly (they are not numbers)" << std::endl;
			return false;
		}  

		return true;
	}

	bool ReadReactionKeyWordPlusWords(const std::string tag, std::string& line, const int n1, const int n2, std::vector<std::string>& words)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);
						
		if (std::distance (tokens_slash.begin(), tokens_slash.end()) != 2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (1)!" << std::endl;
			return false;
		}
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;
	
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());
		if ( n!=n1 && n!=n2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (2)!" << std::endl;
			return false;
		}

		words.resize(n);
		tokenizer_blank::iterator tok_blank = tokens_blank.begin();
		for (std::size_t i = 0; i<n; i++)
		{
			words[i] = *tok_blank;
			++tok_blank;
		}
	
		return true;
	}

	bool ReadReactionKeyWordPlusWords(const std::string tag, std::string& line, std::vector<std::string>& words)
	{
		boost::algorithm::trim(line);
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_slash;
		boost::char_separator<char> sep_slash("/");
		tokenizer_slash tokens_slash(line, sep_slash);
						
		if (std::distance (tokens_slash.begin(), tokens_slash.end()) != 2)
		{
			std::cout << tag << "Syntax error: Wrong number of arguments (1)!" << std::endl;
			return false;
		}
		tokenizer_slash::iterator tok_slash = tokens_slash.begin();
		++tok_slash;
		std::string numbers = *tok_slash;
	
		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_blank;
		boost::char_separator<char> sep_blank(" ");
		tokenizer_blank tokens_blank(numbers, sep_blank);
		const std::size_t n = std::distance(tokens_blank.begin(), tokens_blank.end());

		words.resize(n);
		tokenizer_blank::iterator tok_blank = tokens_blank.begin();
		for (std::size_t i = 0; i<n; i++)
		{
			words[i] = *tok_blank;
			++tok_blank;
		}
	
		return true;
	}


	bool LookForKeyWord(const std::string& line, std::string& keyword)
	{
		size_t found=line.find("/");
		if (found!=std::string::npos)
		{
			keyword = line.substr(0,found);
			boost::algorithm::trim(keyword);
		}
		else
		{
			typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_plus;
			boost::char_separator<char> sep(" ");
			tokenizer_plus tokens(line, sep);
			tokenizer_plus::iterator tok_iter = tokens.begin();
			keyword=*tok_iter; 
		}
		return true;
	}

	bool SeparateReactantSideFromProductSide(const std::string& reconstructedline, std::string &line_reactants, std::string &line_products, bool &iReversible)
	{
		size_t found=reconstructedline.find("<=>");
		if (found!=reconstructedline.npos)
		{
			iReversible = true;
			line_reactants = reconstructedline.substr(0,found);
			line_products = reconstructedline.substr(found+3,reconstructedline.npos);
			return true;
		}
		else
		{
			size_t found=reconstructedline.find("=>");
			if (found!=reconstructedline.npos)
			{
				iReversible = false;
				line_reactants = reconstructedline.substr(0,found);
				line_products = reconstructedline.substr(found+2,reconstructedline.npos);
				return true;
			}
			else
			{
				size_t found=reconstructedline.find("=");
				if (found!=reconstructedline.npos)
				{
					iReversible = true;
					line_reactants = reconstructedline.substr(0,found);
					line_products = reconstructedline.substr(found+1,reconstructedline.npos);
					return true;
				}
				else
				{
					std::cout << "Missing one of these symbols in the reaction: <=> || => || =" << std::endl;
					return false;
				}
			}
		}
	}

}

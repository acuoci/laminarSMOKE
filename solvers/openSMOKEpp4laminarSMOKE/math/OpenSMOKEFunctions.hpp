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

#include <stdlib.h>
#include <time.h>
#include <boost/math/constants/constants.hpp> 
#include <boost/date_time.hpp>
#include "boost/date_time/gregorian/gregorian.hpp"

#if OPENSMOKE_USE_MKL == 1
	#include <mkl.h>
#elif OPENSMOKE_USE_OPENBLAS == 1
	#include "lapacke.h"
#endif

namespace OpenSMOKE
{
	static const float  OPENSMOKE_MACH_EPS_FLOAT = MachEpsFloat();
	static const double OPENSMOKE_MACH_EPS_DOUBLE = MachEps();

	void ErrorMessage(const std::string functionName, const std::string errorMessage)
	{
		std::cout << "Function:     " << functionName << std::endl;
		std::cout << "Fatal error:  " << errorMessage << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
	}

	int FatalErrorMessage(const std::string errorMessage)
	{
		std::cout << "Fatal error:  " << errorMessage << std::endl;
		std::cout << "Press enter to exit..." << std::endl;
		getchar();
		exit(OPENSMOKE_FATAL_ERROR_EXIT);
		return OPENSMOKE_FATAL_ERROR_EXIT;
	}

	float MachEpsFloat()
	{
		float macheps = 1.F;
		float eps = 2.F;
	
		while(eps != 1.F)
		{
			macheps /= 2.F;
			eps = float(1.F + macheps); // in simple rounded form
		}
	
		if(macheps < 1.e-7)
			return 1.192092896e-07F;
		else
			return macheps*2.F;
	}

	double MachEps()
	{
		double macheps = 1.;
		double eps = 2.;
		while(eps != 1.)
		{
			macheps /= 2.;
			eps = 1. + macheps; // in double rounded form
		}
		return macheps*2.;
	}

	float SqrtSumSqr(const int n, float *x)
	{	
		double norm = 0.;
		for(int i=0;i<n;i++)
		{
			double aux = x[i];
			norm += aux*aux;
		}
		norm = std::sqrt(norm);
		
		if(norm > OPENSMOKE_BIG_FLOAT)
			norm = OPENSMOKE_BIG_FLOAT;
		
		return float(norm);
	}

	double SqrtSumSqr(const int n, double *x)
	{
		double aux;
		double xmax = 0.;
		double xmin = OPENSMOKE_BIG_DOUBLE;
		for(int j=0;j<n;j++)
		{
			aux = std::fabs(x[j]);
			if(xmax<aux)	xmax = aux;
			if(xmin>aux)	xmin = aux;
		}
		
		if(xmax == 0.)
			return xmax;
		
		if (xmin == 0.)
			xmin = OPENSMOKE_TINY_DOUBLE;
		
		aux = std::sqrt(OPENSMOKE_BIG_DOUBLE/((double)n));
		
		// to avoid the problems of
		#if OPENSMOKE_LONG_DOUBLE == OPENSMOKE_DOUBLE
			if	( xmax < aux &&												// overflow
				  xmax > OPENSMOKE_TINY_DOUBLE/OPENSMOKE_MACH_EPS_FLOAT)	 // small numbers
		#else
			long double longaux = (long double)xmax/(long double)xmin;
			if(	xmax < aux &&									// overflow
				xmax > OPENSMOKE_TINY_DOUBLE/OPENSMOKE_MACH_EPS_FLOAT &&	// small numbers
				longaux < 1./OPENSMOKE_MACH_EPS_FLOAT)				// sort
		#endif
			{
				double norm = 0.;	// without problems: double
				for(int i = 0;i < n;i++)
				{
					aux = x[i];
					norm += aux*aux;
				}
				return std::sqrt(norm);
			}
			else	// if there are problems it works in long double
			{
				#if OPENSMOKE_LONG_DOUBLE == OPENSMOKE_DOUBLE
		
					double *y = new double[n];
					
					if(!y)
					{
						std::cout << "Not enough memory in double SqrtSumSqr(int n,double *x) function" << std::endl;
						std::cout << "Press enter to exit..." << std::endl;
						getchar();
						exit(OPENSMOKE_FATAL_ERROR_EXIT);
					}
					
					for(int i=0;i<n;i++)
						y[i] = x[i]/xmax;
					double norm = 0.;
					for(int i=0;i<n;i++)
						norm += y[i]*y[i];
					delete y;

					return xmax*std::sqrt(norm);
		
				#else
					long double norm = 0.;
					for(int i = 0;i < n;i++)
					{
						longaux = x[i];
						norm += longaux*longaux;
					}
					if(norm < OPENSMOKE_BIG_DOUBLE && norm > OPENSMOKE_TINY_DOUBLE)
						return std::sqrt(double(norm));
					longaux = (long double)xmax*(long double)n;
					norm /= longaux; // avoids overflow
					norm /= longaux;
					norm = longaux*std::sqrt(double(norm)); // renormalises
					
					// avoids overflow
					if(norm > OPENSMOKE_BIG_DOUBLE) norm = OPENSMOKE_BIG_DOUBLE;
						return double(norm);
				#endif
		}
	}
        
    double OpenSMOKEClock(void)
	{
		return (double)(std::clock())/CLOCKS_PER_SEC;
	}

    double OpenSMOKEGetCpuTime(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return dsecnd();
		#else
            return OpenSMOKEClock();
		#endif
    }

	#if defined(_WIN32) || defined(_WIN64) 
		unsigned __int64 OpenSMOKEGetCpuClocks(void)
	#else
		unsigned long int OpenSMOKEGetCpuClocks(void)
	#endif
    {
		#if OPENSMOKE_USE_MKL == 1
			unsigned MKL_INT64 clocks;
			mkl_get_cpu_clocks( &clocks );
			return clocks;
		#else
            return std::clock();
		#endif	
	}

    double OpenSMOKEGetCpuFrequency(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return mkl_get_cpu_frequency();
		#else
            return 0;
		#endif
	}

    double OpenSMOKEGetMaxCpuFrequency(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return 0;
			//return mkl_get_max_cpu_frequency();
		#else
            return 0;
		#endif
	}

    double OpenSMOKEGetCpuClocksFrequency(void)
    {
		#if OPENSMOKE_USE_MKL == 1
			return 0;
			//return mkl_get_clocks_frequency();
		#else
            return 0;
		#endif
	}

	void OpenOutputFileXML(std::ofstream &fXML, const boost::filesystem::path output_file_xml)
	{
		if (!fXML.is_open())
			fXML.open(output_file_xml.c_str(), std::ios::out);

		if (!fXML.is_open())
		{
			std::string msg = std::string("The file ") + std::string(output_file_xml.string()).c_str() + std::string(" could not be open!");
			FatalErrorMessage(msg.c_str());
		}
	}

	void OpenOutputFileASCII(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii)
	{
		if (!fASCII.is_open())
			fASCII.open(output_file_ascii.c_str(), std::ios::out);

		if (!fASCII.is_open())
		{
			std::string msg = std::string("The file ") + std::string(output_file_ascii.string()).c_str() + std::string(" could not be open!");
			FatalErrorMessage(msg.c_str());
		}

		fASCII.setf(std::ios::scientific);
	}

	void OpenOutputFileASCII_Append(std::ofstream &fASCII, const boost::filesystem::path output_file_ascii)
	{
		if (!fASCII.is_open())
			fASCII.open(output_file_ascii.c_str(), std::ios::out | std::ios::app);

		if (!fASCII.is_open())
		{
			std::string msg = std::string("The file ") + std::string(output_file_ascii.string()).c_str() + std::string(" could not be open!");
			FatalErrorMessage(msg.c_str());
		}

		fASCII.setf(std::ios::scientific);
	}

	void CheckKineticsFolder(const boost::filesystem::path& path_output)
	{
		if (boost::filesystem::is_directory(path_output) == true)
		{
			// In case the folder does not exist
			if (boost::filesystem::exists(path_output) == false)
			{
				std::string message = "The folder " + path_output.string() + " you specified through the @KineticsFolder option does not exist!";
				FatalErrorMessage(message);
			}

			// In case the kinetics.xml file does not exist
			{
				const boost::filesystem::path kinetics_file = path_output / "kinetics.xml";
				if (boost::filesystem::exists(kinetics_file) == false)
				{
					std::string message = "The " + kinetics_file.string() + " file does not exist!\n";
					message += "Please check the content of the folder you specified through the @KineticsFolder option";
					FatalErrorMessage(message);
				}
			}

			// In case the reaction_names.xml file does not exist
			{
				const boost::filesystem::path reaction_names_file = path_output / "reaction_names.xml";
				if (boost::filesystem::exists(reaction_names_file) == false)
				{
					std::string message = "The " + reaction_names_file.string() + " file does not exist!\n";
					message += "Please check the content of the folder you specified through the @KineticsFolder option";
					FatalErrorMessage(message);
				}
			}
		}
		else
		{
			std::string message = "The folder " + path_output.string() + " you specified through the @KineticsFolder option does not exist!";
			FatalErrorMessage(message);
		}
	}

	void OpenInputFileXML(rapidxml::xml_document<>& doc, std::vector<char>& xml_copy, const boost::filesystem::path& file_name)
	{
		if (!boost::filesystem::exists(file_name))
		{
			std::string message = "The " + file_name.string() + " file does not exist";
			FatalErrorMessage(message);
		}

		std::ifstream fInput;
		fInput.open(std::string(file_name.string()).c_str(), std::ios::in);

		if (!fInput.is_open())
		{
			std::string msg = std::string("The ") + std::string(file_name.string()).c_str() + std::string(" file cannot be open!");
			FatalErrorMessage(msg.c_str());
		}

		const std::string string_file = std::string(std::istreambuf_iterator<char>(fInput), std::istreambuf_iterator<char>());
		fInput.close();

		xml_copy.assign(string_file.begin(), string_file.end());
		xml_copy.push_back('\0');

		doc.parse<rapidxml::parse_declaration_node | rapidxml::parse_no_data_nodes>(&xml_copy[0]);
	}

	void OpenInputFileASCII(std::ifstream &fASCII, const boost::filesystem::path input_file_ascii)
	{
		if (!boost::filesystem::exists(input_file_ascii))
		{
			std::string message = "The " + input_file_ascii.string() + " file does not exist";
			FatalErrorMessage(message);
		}

		fASCII.open(input_file_ascii.c_str(), std::ios::in);

		if (!fASCII.is_open())
		{
			std::string msg = std::string("The ") + std::string(input_file_ascii.string()).c_str() + std::string(" file cannot be open!");
			FatalErrorMessage(msg.c_str());
		}
	}

	void SetXMLFile(std::ostream& xml_string)
	{
		xml_string << std::setprecision(8);
		xml_string.setf(std::ios::scientific);

		xml_string << "<?xml version=\"1.0\" encoding=\"utf-8\"?>" << std::endl;
		xml_string << "<opensmoke version=\"0.1a\">" << std::endl;
	}

	void PrintTagOnASCIILabel(const unsigned int width, std::ostream &fOut, const std::string tag, unsigned int &counter)
	{
		std::stringstream number;
		number << counter++;
		std::string label = tag + "(" + number.str() + ")";
		fOut << std::setw(width) << std::left << label;
	}

	void CreateDirectory(const boost::filesystem::path& path_output)
	{
		boost::filesystem::create_directories(path_output);

		if (!boost::filesystem::exists(path_output))
		{
			std::cout << "Error while attempting to create the " << path_output.string() << " directory." << std::endl;
			std::cout << "Please check the permissions you have on the folder." << std::endl;
			std::cout << "Press enter to exit..." << std::endl;
			getchar();
			exit(OPENSMOKE_FATAL_ERROR_EXIT);
		}
	}

	std::string GetCurrentTime()
	{
		std::string current_time = boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());
		current_time.erase(0, 9);
		current_time.insert(2, ":");
		current_time.insert(5, ":");

		return current_time;
	}

	std::string GetRawCurrentTime()
	{
		std::string current_time = boost::posix_time::to_iso_string(boost::posix_time::second_clock::local_time());
		current_time.erase(0, 9);
		return current_time;
	}

	std::string GetCurrentDate()
	{
		std::string current_date = boost::gregorian::to_iso_string(boost::posix_time::second_clock::local_time().date());
		current_date.insert(4, "-");
		current_date.insert(7, "-");

		return current_date;
	}

	std::string SplitStringIntoSeveralLines(std::string source, const std::size_t width, const std::string whitespace)
	{
		std::size_t  currIndex = width - 1;
		std::size_t  sizeToElim;
		while (currIndex < source.length())
		{
			currIndex = source.find_last_of(whitespace, currIndex + 1);
			if (currIndex == std::string::npos)
				break;
			currIndex = source.find_last_not_of(whitespace, currIndex);
			if (currIndex == std::string::npos)
				break;
			sizeToElim = source.find_first_not_of(whitespace, currIndex + 1) - currIndex - 1;
			source.replace(currIndex + 1, sizeToElim, "\n");
			currIndex += (width + 1);
		}
		return source;
	}

	boost::filesystem::path GetExecutableFileName(char** argv)
	{
		#ifdef __linux
		{
			char buf[1024];
			memset(buf, 0, sizeof(buf));
			if (!readlink("/proc/self/exe", buf, sizeof(buf)-1))
			{
				perror("readlink");
				exit(OPENSMOKE_FATAL_ERROR_EXIT);
			}
			boost::filesystem::path executable_file = buf;
			return executable_file;
		}
		#elif __APPLE__ 
		{
			char path[1024];
			uint32_t size = sizeof(path);
			if (_NSGetExecutablePath(path, &size) == 0)
			{
				boost::filesystem::path executable_file = path;
				return executable_file;
			}
			else
			{
				perror("readlink");
				exit(OPENSMOKE_FATAL_ERROR_EXIT);
			}
		}
		#else
		{
			return boost::filesystem::system_complete(argv[0]);
		}
		#endif
	}

	void OpenSMOKE_logo(const std::string application_name, const std::string author_name)
	{
		std::string current_time = __TIME__;
		std::string current_date = __DATE__;
		std::string author_complete = "Author: " + author_name;
		std::string compilation_time = "Compilation date: " + current_date + " at " + current_time;
		std::string version = "Version: "; version += __OPENSMOKE_VERSION__;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << std::endl;
		std::cout << "          ___                   ____  __  __  ___  _  _______                 " << std::endl;
		std::cout << "         / _ \\ _ __   ___ _ __ / ___||  \\/  |/ _ \\| |/ / ____| _     _        " << std::endl;
		std::cout << "        | | | | '_ \\ / _ \\ '_ \\\\___ \\| |\\/| | | | | ' /|  _| _| |_ _| |_      " << std::endl;
		std::cout << "        | |_| | |_) |  __/ | | |___) | |  | | |_| | . \\| |__|_   _|_   _|     " << std::endl;
		std::cout << "         \\___/| .__/ \\___|_| |_|____/|_|  |_|\\___/|_|\\_\\_____||_|   |_|       " << std::endl;
		std::cout << "              |_|                                                             " << std::endl;
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "           Department of Chemistry, Materials and Chemical Engineering        " << std::endl;
		std::cout << "                              Politecnico di Milano                           " << std::endl;
		std::cout << "                         http://www.opensmoke.polimi.it/                      " << std::endl;
		std::cout << "                      http://creckmodeling.chem.polimi.it/                    " << std::endl;
		std::cout << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(application_name.size()) / 2; i++)	std::cout << " ";
		std::cout << application_name << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(version.size()) / 2; i++)	std::cout << " ";
		std::cout << version << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(author_complete.size()) / 2; i++)	std::cout << " ";
		std::cout << author_complete << std::endl;
		for (unsigned int i = 1; i <= 39 - (unsigned int)(compilation_time.size()) / 2; i++)	std::cout << " ";
		std::cout << compilation_time << std::endl;
		std::cout << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                                  WARNING                                    " << std::endl;
		std::cout << "   This version of OpenSMOKE++ Suite can be used for educational purposes    " << std::endl;
		std::cout << "              only and cannot be distributed to third parties.               " << std::endl;
		std::cout << "       The software is and remains the sole property of Alberto Cuoci.       " << std::endl;
		std::cout << "      Whenever the OpenSMOKE++ Suite is used to produce any publication,     " << std::endl;
		std::cout << "       a detailed reference to the OpenSMOKE++ code should be reported       " << std::endl;
		std::cout << "                            (see User's Guide).                              " << std::endl;
		std::cout << "    Use for commercial purposes is not permitted. For any commercial issue   " << std::endl;
		std::cout << "         please contact Alberto Cuoci (email: alberto.cuoci@polimi.it)       " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;

		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
		std::cout << "                            LIMITED WARRANTY                                 " << std::endl;
		std::cout << "     This software is provided \"as is\" and without warranties as to        " << std::endl;
		std::cout << "  performance of merchantability or any other warranties whether expressed   " << std::endl;
		std::cout << "    or implied. Because of the various hardware and software environments    " << std::endl;
		std::cout << "   into which this library may be installed, no warranty of fitness for a    " << std::endl;
		std::cout << "   particular purpose is offered. The user must assume the entire risk of    " << std::endl;
		std::cout << "                          using  the library.                                " << std::endl;
		std::cout << "-----------------------------------------------------------------------------" << std::endl;
		std::cout << "" << std::endl;
	}

	void OpenSMOKE_logo(std::ofstream& fOut, const std::string application_name)
	{
		std::string current_time = __TIME__;
		std::string current_date = __DATE__;
		std::string compilation_time = "Compilation date: " + current_date + " at " + current_time;
		std::string version = "Version: "; version += __OPENSMOKE_VERSION__;

		fOut << "***********************************************************************************************************" << std::endl;
		fOut << "***********************************************************************************************************" << std::endl;
		fOut << std::endl;
		fOut << "                     ___                   ____  __  __  ___  _  _______                 " << std::endl;
		fOut << "                    / _ \\ _ __   ___ _ __ / ___||  \\/  |/ _ \\| |/ / ____| _     _        " << std::endl;
		fOut << "                   | | | | '_ \\ / _ \\ '_ \\\\___ \\| |\\/| | | | | ' /|  _| _| |_ _| |_      " << std::endl;
		fOut << "                   | |_| | |_) |  __/ | | |___) | |  | | |_| | . \\| |__|_   _|_   _|     " << std::endl;
		fOut << "                    \\___/| .__/ \\___|_| |_|____/|_|  |_|\\___/|_|\\_\\_____||_|   |_|       " << std::endl;
		fOut << "                         |_|                                                             " << std::endl;
		fOut << std::endl;
		fOut << std::endl;
		fOut << "                      Department of Chemistry, Materials and Chemical Engineering        " << std::endl;
		fOut << "                                         Politecnico di Milano                           " << std::endl;
		fOut << "                                 http://creckmodeling.chem.polimi.it/                    " << std::endl;
		fOut << std::endl;
		for (unsigned int i = 1; i <= 50 - (unsigned int)(application_name.size()) / 2; i++)	fOut << " ";
		fOut << application_name << std::endl;
		for (unsigned int i = 1; i <= 50 - (unsigned int)(version.size()) / 2; i++)	fOut << " ";
		fOut << version << std::endl;
		for (unsigned int i = 1; i <= 50 - (unsigned int)(compilation_time.size()) / 2; i++)	fOut << " ";
		fOut << compilation_time << std::endl;
		fOut << std::endl;
		fOut << "***********************************************************************************************************" << std::endl;
		fOut << "***********************************************************************************************************" << std::endl;
		fOut << std::endl;
	}

	unsigned int CalculateSpeciesFieldWidth(const std::string species_name, const int max_number_species)
	{
		const double coefficient = 2.3*max_number_species;
		const std::size_t max_number_of_digits = boost::lexical_cast<std::size_t>(std::ceil(std::log10(coefficient)));
		const std::size_t safety_digits = 2;
		const std::size_t minimum_number = 16;

		// name + _ + W + ( max_digits ) + safety 
		std::size_t number_to_compare = species_name.size() + 3 + max_number_of_digits + 1 + safety_digits;
		return boost::lexical_cast<unsigned int>(std::max(minimum_number, number_to_compare));
	}

	void ImportReactionNames(const boost::filesystem::path file_name, const unsigned int n, std::vector<std::string>& reaction_strings)
	{
		if (!boost::filesystem::exists(file_name))
		{
			std::string message = "The " + file_name.string() + " file does not exist";
			FatalErrorMessage(message);
		}

		rapidxml::xml_document<> local_xml;
		std::vector<char> local_xml_input_string;
		OpenSMOKE::OpenInputFileXML(local_xml, local_xml_input_string, file_name);

		// Names of reactions
		{
			rapidxml::xml_node<>* reaction_names_node = local_xml.first_node("opensmoke")->first_node("reaction-names");
			std::stringstream values(reaction_names_node->value());

			reaction_strings.reserve(n);
			for (unsigned int j = 0; j<n; j++)
			{
				std::string reaction_string;
				values >> reaction_string;
				reaction_strings.push_back(reaction_string);
			}
		}
	}

	std::vector<double> CubicRootsReal(const double a, const double b, const double c, const double d)
	{
		std::vector<double> solution;

		const double alfa = b / a;
		const double beta = c / a;
		const double gamma = d / a;

		const double p = beta - alfa * alfa / 3.;
		const double q = 2. * alfa * alfa * alfa / 27. - alfa * beta / 3. + gamma;

		const double D = q * q / 4. + p * p * p / 27.;

		if (D > 0.)
		{
			solution.resize(1);
			solution[0] = cbrt(-0.5 * q + std::sqrt(D)) +
				cbrt(-0.5 * q - std::sqrt(D)) - alfa / 3.;
		}
		else if (D == 0.)
		{
			solution.resize(2);
			solution[0] = -2. * cbrt(q / 2.) - alfa / 3.;
			solution[1] = cbrt(q / 2.) - alfa / 3.;
		}
		else if (D < 0.)
		{
			solution.resize(3);
			double r = std::sqrt(-p * p * p / 27.);
			double costheta = -q / (2. * r);
			double theta = std::acos(costheta);
			solution[0] = 2. * cbrt(r) * std::cos(theta / 3.) - alfa / 3.;
			solution[1] = 2. * cbrt(r) * std::cos((2. * boost::math::constants::pi<double>() + theta) / 3.) - alfa / 3.;
			solution[2] = 2. * cbrt(r) * std::cos((4. * boost::math::constants::pi<double>() + theta) / 3.) - alfa / 3.;
		}

		return solution;
	}

	std::vector<double> Normalize(const std::vector<double>& x)
	{
		double x_sum = 0.;

		for (unsigned int i = 0; i < x.size(); i++)
			x_sum += x[i];

		std::vector<double> x_tmp(x.size());
		for (unsigned int i = 0; i < x.size(); i++)
			x_tmp[i] = x[i] / x_sum;

		return x_tmp;
	}

	double Median(const std::vector<double>& vec)
	{
		std::vector<double> vec_sorted = vec;
		double median;

		std::sort(vec_sorted.begin(), vec_sorted.begin() + vec_sorted.size());


		if (vec_sorted.size() % 2 == 0)
		{
			median = vec_sorted[vec_sorted.size() / 2 - 1];
			median += vec_sorted[vec_sorted.size() / 2];
			median /= 2;
		}

		else if (vec_sorted.size() % 2 == 1)
			median = vec_sorted[vec_sorted.size() / 2];

		vec_sorted.clear();

		return median;
	}

	double MedianAbsoluteDeviation(const std::vector<double>& vec)
	{
		if (vec.size() <= 1)
			OpenSMOKE::FatalErrorMessage("Variance cannot be evaluated if a vector has less than 2 elements!");

		const double med = Median(vec);
		double mad_value = 0.;

		std::vector<double> deviations(vec.size());
		for (unsigned int i = 0; i < vec.size(); i++)
			deviations[i] = std::fabs(vec[i] - med);

		mad_value = Median(deviations);

		return mad_value;
	}

	void CheckAndCorrectSumOfFractions(std::vector<double>& x, const double boundary_eps, const double eps)
	{
		for (unsigned int i = 0; i < x.size(); i++)
		{
			if (x[i] < -boundary_eps)
			{
				std::cerr << "Fatal Error: Negative mole fractions! " << std::endl;
				exit(-1);
			}
			else if (x[i] > 1. + boundary_eps)
			{
				std::cerr << "Fatal Error: Mole fractions higher than 1! " << std::endl;
				exit(-1);
			}
		}

		const double sum = std::accumulate(x.begin(), x.end(), 0.);
		if (sum > 1. + eps || sum < 1. - eps)
		{
			std::cerr << "Fatal Error: Sum of mass fractions is " << sum << std::endl;
			exit(-1);
		}
		else
		{
			for (unsigned int i = 0; i < x.size(); i++)
				x[i] /= sum;
		}
	}

	unsigned int NumberOfLinesInFile(const boost::filesystem::path& path)
	{
		std::ifstream file(path.c_str());

		return static_cast<unsigned int>(std::count(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(), '\n'));
	}

	unsigned int CountWordsInString(std::stringstream& str)
	{
		unsigned int n_elements = 0;

		std::string element;
		while (str >> element)
			++n_elements;

		return n_elements;
	}

	void SparsityPatternTridiagonal(const unsigned int number_equations, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		const unsigned int nonzeros = 3 * number_equations - 2;
		rows.resize(nonzeros);
		cols.resize(nonzeros);

		int count = 0;

		// Main diagonal
		for (unsigned int i = 1; i <= number_equations; i++)
		{
			rows[count] = i;
			cols[count] = i;
			count++;
		}

		// Adjacent diagonals
		for (unsigned int i = 1; i <= number_equations - 1; i++)
		{
			// upper diagonal
			rows[count] = i;
			cols[count] = i + 1;
			count++;

			// lower diagonal
			rows[count] = i + 1;
			cols[count] = i;
			count++;
		}
	}

	void SparsityPatternPentadiagonal(const unsigned int number_equations, const unsigned int width, std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		const unsigned int nonzeros = 5 * number_equations - 2 * width - 2;
		rows.resize(nonzeros);
		cols.resize(nonzeros);

		int count = 0;

		// Main diagonal
		for (unsigned int i = 1; i <= number_equations; i++)
		{
			rows[count] = i;
			cols[count] = i;
			count++;
		}

		// Adjacent diagonals
		for (unsigned int i = 1; i <= number_equations - 1; i++)
		{
			// upper diagonal
			rows[count] = i;
			cols[count] = i + 1;
			count++;

			// lower diagonal
			rows[count] = i + 1;
			cols[count] = i;
			count++;
		}

		// External diagonals
		for (unsigned int i = 1; i <= number_equations - width; i++)
		{
			// upper external diagonal
			rows[count] = i;
			cols[count] = i + width;
			count++;

			// lower external diagonal
			rows[count] = i + width;
			cols[count] = i;
			count++;
		}
	}

	void SparsityPatternBlock(	const unsigned int number_blocks, const unsigned int block_size,
								const std::vector<unsigned int>& rows_single, const std::vector<unsigned int>& cols_single,
								std::vector<unsigned int>& rows, std::vector<unsigned int>& cols)
	{
		const unsigned int block_squared = block_size*block_size;
		const unsigned int nonzeros = rows_single.size()*block_squared;
		rows.resize(nonzeros);
		cols.resize(nonzeros);

		int count = 0;
		for (unsigned int k = 0; k < rows_single.size(); k++)
		{
			for (unsigned int i = 0; i < block_size; i++)
			{
				for (unsigned int j = 0; j < block_size; j++)
				{
					rows[count] = (rows_single[k] - 1) * block_size + i + 1;
					cols[count] = (cols_single[k] - 1) * block_size + j + 1;
					count++;
				}
			}
		}
	}
}

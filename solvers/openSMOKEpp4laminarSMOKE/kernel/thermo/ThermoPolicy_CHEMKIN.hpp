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

#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <string>
#include <Eigen/Dense>

int QUADROOTS(double p[5], double r[3][5]);
int CUBICROOTS(double p[], double r[3][5]);
int BIQUADROOTS(double p[5], double r[3][5]);
void FindCubicRoots(const double A, const double B, const double C, const double D, unsigned int& n, double x[3]);
void FindBiquadricRoots(const double A, const double B, const double C, const double D, const double E, unsigned int& n, double x[4]);
std::string ScientificNotationWithFixedExponentDigits(const double number, const int expSize);

bool IsInteger(const double number)
{
  return std::floor(number) == number;
}

namespace OpenSMOKE
{

	ThermoPolicy_CHEMKIN::ThermoPolicy_CHEMKIN() {
	}

	ThermoPolicy_CHEMKIN::ThermoPolicy_CHEMKIN(const ThermoPolicy_CHEMKIN& orig)
	{
		name_thermo_  = orig.name_thermo_;
		phase_ = orig.phase_;
		thermodynamics_parameters_ = orig.thermodynamics_parameters_;
		atomic_composition_ = orig.atomic_composition_;
	}

	ThermoPolicy_CHEMKIN::~ThermoPolicy_CHEMKIN() {
	}

	void ThermoPolicy_CHEMKIN::CopyThermodynamics(ThermoPolicy_CHEMKIN& rhs) const 
	{
		rhs.name_thermo_ = name_thermo_;
		rhs.phase_ = phase_;
		rhs.atomic_composition_ = atomic_composition_;
		rhs.thermodynamics_parameters_ = thermodynamics_parameters_;
	}

	void ThermoPolicy_CHEMKIN::AtomicComposition(OpenSMOKE::AtomicComposition *atomic) const
	{ 
		*atomic = atomic_composition_; 
	}

	bool ThermoPolicy_CHEMKIN::ReadThermodynamicsFromStrings(const std::vector<std::string> &line, const OpenSMOKE::OpenSMOKEVectorDouble& temperature_limits_default_)
	{
			std::vector<std::string>        element_names_;
			std::vector<double>             element_coefficients_;
			OpenSMOKE::OpenSMOKEVectorDouble	thermodynamics_coefficients_;
			OpenSMOKE::OpenSMOKEVectorDouble    temperature_limits_;

			try
			{
				if (line[0].size()<80)
					throw "Expected 1 at column 80";
				if (line[0].at(79) != '1')
					throw "Expected 1 at column 80";
			}
			catch(char const *msg)
			{
				std::cout << msg << std::endl;
				return false;
			}
        
			// Analyze the first line
			int iFirstLineThermodynamics=1;
			{
				int offsets_firstline[] = {16,2,6,20,1,10,10,8,5,1,1,20};
				boost::offset_separator f_firstline(offsets_firstline, offsets_firstline+12);
				boost::tokenizer<boost::offset_separator> tok(line[0],f_firstline);

				std::vector<std::string> words(12);
				int j=0;
				for(boost::tokenizer<boost::offset_separator>::iterator beg=tok.begin(); beg!=tok.end();++beg)
					words[j++] = *beg;

				name_thermo_  = words[0];
				boost::to_upper(name_thermo_);
				std::string date_  = words[2];
				phase_ = words[4].at(0);
				std::string atomic_ = "";
				std::string additional_atomic_ = "";

				atomic_ = words[3]+words[8];

				if (words[11].size() > 0)
				{
					if (words[11].at(0) == '&')
					{
						additional_atomic_ = line[1];
						iFirstLineThermodynamics += 1;
					}
					else
					{
                                                // If additional characters are reported after the 80th position
                                                // they are ignored
                                            
						//atomic_ += words[11];
					}
				}

				// Recovering the elemental composition
				{
						std::vector<std::string>        element_numbers;
						int offsets_atomic[] = {2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3,2,3};
						boost::offset_separator f_atomic(offsets_atomic, offsets_atomic+18);
						boost::tokenizer<boost::offset_separator> tok(atomic_,f_atomic);

						for(boost::tokenizer<boost::offset_separator>::iterator beg=tok.begin(); beg!=tok.end();++beg)
						{  
							std::string name=*beg;
							beg++;
							std::string number=*beg;
							boost::algorithm::trim(name);
							boost::algorithm::trim(number);

							if (name.size()!=0 && std::isdigit(name[0])==false)
							{
								boost::to_upper(name);
								element_names_.push_back(name);
								element_numbers.push_back(number);
							}
						}

						if (additional_atomic_.size()!=0)
						{
								typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_only_blanks;
								boost::char_separator<char> sep(" ");
								tokenizer_only_blanks tokens(additional_atomic_, sep);
								for (tokenizer_only_blanks::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
								{
									std::string name_aux = boost::to_upper_copy(*tok_iter);
									element_names_.push_back(name_aux);
									tok_iter++;
									element_numbers.push_back(*tok_iter);
								}   
						}
               
						try
						{
							element_coefficients_.resize(element_numbers.size());
							for(unsigned int j=0;j<element_numbers.size();j++)
							{
									boost::algorithm::trim(element_numbers[j]);
									element_coefficients_[j] = boost::lexical_cast<double>(element_numbers[j]);
							}
						}
						catch(boost::bad_lexical_cast &)
						{
							std::cout << "Numerical conversion failure at the first line. The atomic coefficients are not written properly (they are not numbers)" << std::endl;
							return false;
						}  
				}
            
				// Checking temperatures
				try
				{
					boost::algorithm::trim(words[5]);
					boost::algorithm::trim(words[6]);
					boost::algorithm::trim(words[7]);
					ChangeDimensions(3, &temperature_limits_, true);
                                        temperature_limits_[1] = temperature_limits_default_[1];
                                        temperature_limits_[2] = temperature_limits_default_[2];
                                        temperature_limits_[3] = temperature_limits_default_[3];
                                        
					if (!words[5].empty())
						temperature_limits_[1] = boost::lexical_cast<double>(words[5]);  // low
					if (!words[6].empty())
						temperature_limits_[3] = boost::lexical_cast<double>(words[6]);  // high
					if (!words[7].empty())
						temperature_limits_[2] = boost::lexical_cast<double>(words[7]);  // medium
				}
				catch(boost::bad_lexical_cast &)
				{
					std::cout << "Numerical conversion failure at the first line. The temperatures are not written properly." << std::endl;
					return false;
				}

				try
				{ 
					if ( (temperature_limits_[1] > temperature_limits_[2]) || (temperature_limits_[2] > temperature_limits_[3]) )
						throw "Wrong sequence of temperature limits. Please provide temperatures in the following order: Tmin Tmax Tmedium.";
				}
				catch(char const *msg)
				{
					std::cout << msg << std::endl;
					return false;
				}
        
			}
        

			std::vector<std::string>        thermodynamics_numbers;
			{
				typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_only_blanks;
				boost::char_separator<char> sep(" ");
				tokenizer_only_blanks tokens(line[iFirstLineThermodynamics], sep);
            
				// Additional thermodynamics is provided
				if (*tokens.begin()=="TEMP")
				{
					std::vector<std::string> temperature_numbers;
					ChangeDimensions(0,&temperature_limits_, true);
					for (tokenizer_only_blanks::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
					{
						if (tok_iter != tokens.begin())
							temperature_numbers.push_back(*tok_iter);
					}
					try
					{
						for (unsigned int j=0;j<temperature_numbers.size();j++)
						{
							boost::algorithm::trim(temperature_numbers[j]);
							temperature_limits_.Append(boost::lexical_cast<double>(temperature_numbers[j]));
						}
					}
					catch(boost::bad_lexical_cast &)
					{
						std::cout << "Numerical conversion failure of temperature specified through the TEMP keyword." << std::endl;
						return false;
					}
                
					for (unsigned int j=0;j<temperature_numbers.size()-1;j++)
					{
							int offsets_first_line[] = {15,15,15,15,15};
							int offsets_second_line[] = {15,15};
							boost::offset_separator f_first_line(offsets_first_line, offsets_first_line+5, false);
							boost::offset_separator f_second_line(offsets_second_line, offsets_second_line+2, false);
							boost::tokenizer<boost::offset_separator> tok1(line[iFirstLineThermodynamics+1+j*2],f_first_line);
							boost::tokenizer<boost::offset_separator> tok2(line[iFirstLineThermodynamics+2+j*2],f_second_line);
							for(boost::tokenizer<boost::offset_separator>::iterator beg=tok1.begin(); beg!=tok1.end();++beg)
									thermodynamics_numbers.push_back(*beg);
							for(boost::tokenizer<boost::offset_separator>::iterator beg=tok2.begin(); beg!=tok2.end();++beg)
									thermodynamics_numbers.push_back(*beg);                        
					}
				}
				else // Standard thermodynamics is provided
				{
					thermodynamics_numbers.resize(14);

					try
					{
						if (line[iFirstLineThermodynamics].size()<80)
							throw "Expected 2 at column 80";
						if (line[iFirstLineThermodynamics+1].size()<80)
							throw "Expected 3 at column 80";
						if (line[iFirstLineThermodynamics+2].size()<80)
							throw "Expected 4 at column 80";
					}
					catch(char const *msg)
					{
						std::cout << msg << std::endl;
						return false;
					}
                
					try
					{
						if (line[iFirstLineThermodynamics].at(79) != '2')
							throw "Expected 2 at column 80";
						if (line[iFirstLineThermodynamics+1].at(79) != '3')
							throw "Expected 3 at column 80";
						if (line[iFirstLineThermodynamics+2].at(79) != '4')
							throw "Expected 4 at column 80";
					}
					catch(char const *msg)
					{
						std::cout << msg << std::endl;
						return false;
					}
                   
					int offsets_second_third_line[] = {15,15,15,15,15};
					int offsets_fourth_line[] = {15,15,15,15};
					boost::offset_separator f_second_third_line(offsets_second_third_line, offsets_second_third_line+5, false);
					boost::offset_separator f_fourth_line(offsets_fourth_line, offsets_fourth_line+4, false);

					boost::tokenizer<boost::offset_separator> tok2(line[iFirstLineThermodynamics],f_second_third_line);
					boost::tokenizer<boost::offset_separator> tok3(line[iFirstLineThermodynamics+1],f_second_third_line);
					boost::tokenizer<boost::offset_separator> tok4(line[iFirstLineThermodynamics+2],f_fourth_line);

					int j=0; 
					for(boost::tokenizer<boost::offset_separator>::iterator beg=tok2.begin(); beg!=tok2.end();++beg)
						thermodynamics_numbers[j++] = *beg;
					for(boost::tokenizer<boost::offset_separator>::iterator beg=tok3.begin(); beg!=tok3.end();++beg)
						thermodynamics_numbers[j++] = *beg;
					for(boost::tokenizer<boost::offset_separator>::iterator beg=tok4.begin(); beg!=tok4.end();++beg)
						thermodynamics_numbers[j++] = *beg;    
				}
            
				ChangeDimensions(boost::lexical_cast<int>(thermodynamics_numbers.size()), &thermodynamics_coefficients_, true);
            
				// Conversion from strings to numbers
				try
				{
					for(unsigned int j=0;j<thermodynamics_numbers.size();j++)
					{
						boost::algorithm::replace_first(thermodynamics_numbers[j], "e 0", "E+0");
						boost::algorithm::replace_first(thermodynamics_numbers[j], "E 0", "E+0");
						boost::algorithm::trim(thermodynamics_numbers[j]);
						thermodynamics_coefficients_[j+1]  = boost::lexical_cast<double>(thermodynamics_numbers[j]);
					}
				}
				catch(boost::bad_lexical_cast &)
				{
					std::cout << "Numerical conversion failure for the thermodynamic coefficients. Probably they are not written properly (they are not numbers)" << std::endl;
					return false;
				}
			}
        
			// This is not CHEMKIN standard, but many kinetic schemes are written in this way
			{
					typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_only_blanks;
					boost::char_separator<char> sep(" ");
					tokenizer_only_blanks tok(name_thermo_, sep);
					name_thermo_  = *tok.begin();
					boost::algorithm::trim(name_thermo_);
			}
		
			ChangeDimensions(temperature_limits_.Size() + thermodynamics_coefficients_.Size(), &thermodynamics_parameters_, true);
			for(int i=1;i<=temperature_limits_.Size();i++)
				thermodynamics_parameters_[i] = temperature_limits_[i];
			for(int i=temperature_limits_.Size()+1;i<=thermodynamics_parameters_.Size();i++)
				thermodynamics_parameters_[i] = thermodynamics_coefficients_[i-temperature_limits_.Size()];

			atomic_composition_(element_names_, element_coefficients_);
			if (atomic_composition_.CheckAtomicComposition()==false)
				return false;

			return true;
	}

	double ThermoPolicy_CHEMKIN::cp(const double T) const
	{
			double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

			if(T > thermodynamics_parameters_[2])
				return ( thermodynamics_parameters_[4] + T*(thermodynamics_parameters_[5] + T*(thermodynamics_parameters_[6] + T*(thermodynamics_parameters_[7] + T*thermodynamics_parameters_[8]))) )*RGAS_over_MW_;
			else
				return ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13] + T*(thermodynamics_parameters_[14] + T*thermodynamics_parameters_[15]))) )*RGAS_over_MW_;
	}
   
	double ThermoPolicy_CHEMKIN::cv(const double T) const
	{
		return cp(T)-PhysicalConstants::R_J_kmol/MolecularWeight();
	}

	double ThermoPolicy_CHEMKIN::enthalpy(const double T) const
	{
			double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

			if(T > thermodynamics_parameters_[2])
				return ( ( thermodynamics_parameters_[4]  + T*(thermodynamics_parameters_[5]/2. +  T*(thermodynamics_parameters_[6]/3.  + T*(thermodynamics_parameters_[7]/4.  + T*thermodynamics_parameters_[8]/5.))) ) + thermodynamics_parameters_[9]/T)  *T * RGAS_over_MW_;
			else
				return ( ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12]/2. + T*(thermodynamics_parameters_[13]/3. + T*(thermodynamics_parameters_[14]/4. + T*thermodynamics_parameters_[15]/5.)))) + thermodynamics_parameters_[16]/T) *T * RGAS_over_MW_;
	}

	double ThermoPolicy_CHEMKIN::entropy(const double T) const
	{
			double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

			if(T > thermodynamics_parameters_[2])
				return ( ( thermodynamics_parameters_[4]*std::log(T)  + T*(thermodynamics_parameters_[5] +  T*(thermodynamics_parameters_[6]/2.  + T*(thermodynamics_parameters_[7]/3.  + T*thermodynamics_parameters_[8]/4.))) ) + thermodynamics_parameters_[10]) * RGAS_over_MW_;
			else
				return ( ( thermodynamics_parameters_[11]*std::log(T) + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13]/2. + T*(thermodynamics_parameters_[14]/3. + T*thermodynamics_parameters_[15]/4.)))) + thermodynamics_parameters_[17]) * RGAS_over_MW_;
	}

	double ThermoPolicy_CHEMKIN::gibbs_energy(const double T) const
	{
		return enthalpy(T) - T*entropy(T);
	}

	bool ThermoPolicy_CHEMKIN::WriteOnBinaryFile(std::ofstream& fOutput )
	{
		if(!fOutput.write( reinterpret_cast < char * > (&name_thermo_), sizeof(name_thermo_)))
					return false;
	
		if(!fOutput.write( reinterpret_cast < char * > (&phase_), sizeof(phase_)))
				return false;
	
		atomic_composition_.Save(fOutput, OPENSMOKE_BINARY_FILE);

		thermodynamics_parameters_.Save(fOutput, OPENSMOKE_BINARY_FILE);

		return true;
	}

	bool ThermoPolicy_CHEMKIN::ReadFromBinaryFile(std::ifstream& fInput)
	{
		if(!fInput.read(reinterpret_cast < char * > (&name_thermo_), sizeof(name_thermo_)))
					return false;

		if(!fInput.read(reinterpret_cast < char * > (&phase_), sizeof(phase_)))
					return false;
	
		atomic_composition_.Load(fInput, OPENSMOKE_BINARY_FILE);

		thermodynamics_parameters_.Load(fInput, OPENSMOKE_BINARY_FILE);

		return true;
	}

	bool ThermoPolicy_CHEMKIN::WriteOnASCIIFileOldStyle(std::ofstream& fOutput) const
	{
		double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

		// High temperature Cp
		fOutput << 5 << std::endl;
		for (int i=4;i<=8;i++)
			fOutput << std::setprecision(16) << thermodynamics_parameters_[i]*RGAS_over_MW_ << " ";
		fOutput << std::endl;

		// Low temperature Cp
		fOutput << 5 << std::endl;
		for (int i=11;i<=15;i++)
			fOutput << std::setprecision(16) << thermodynamics_parameters_[i]*RGAS_over_MW_ << " ";
		fOutput << std::endl;

		// High temperature H
		fOutput << 6 << std::endl;
		fOutput << std::setprecision(16) << thermodynamics_parameters_[4] << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[5]/2. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[6]/3. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[7]/4. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[8]/5. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[9] << " ";
		fOutput << std::endl;

		// Low temperature H
		fOutput << 6 << std::endl;
		fOutput << std::setprecision(16) << thermodynamics_parameters_[11] << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[12]/2. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[13]/3. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[14]/4. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[15]/5. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[16] << " ";
		fOutput << std::endl;

		// High temperature S
		fOutput << 6 << std::endl;
		fOutput << std::setprecision(16) << thermodynamics_parameters_[4] << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[5] << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[6]/2. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[7]/3. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[8]/4. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[10] << " ";
		fOutput << std::endl;

		// Low temperature S
		fOutput << 6 << std::endl;
		fOutput << std::setprecision(16) << thermodynamics_parameters_[11] << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[12] << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[13]/2. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[14]/3. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[15]/4. << " ";
		fOutput << std::setprecision(16) << thermodynamics_parameters_[17] << " ";
		fOutput << std::endl;

		return true;
	}

	bool ThermoPolicy_CHEMKIN::WriteOnASCIIFile(std::ostream& fOutput) const
	{
		// High temperature
		for (int i=4;i<=10;i++)
			fOutput << std::setprecision(18) << thermodynamics_parameters_[i] << " ";
		
		// Low temperature
		for (int i=11;i<=17;i++)
			fOutput << std::setprecision(18) << thermodynamics_parameters_[i] << " ";

		// Temperatures (L/H/M)
		fOutput << std::setprecision(18) << thermodynamics_parameters_[1] << " ";
		fOutput << std::setprecision(18) << thermodynamics_parameters_[3] << " ";
		fOutput << std::setprecision(18) << thermodynamics_parameters_[2] << " ";

		// Molecular weight
		fOutput << std::setprecision(18) << MolecularWeight() << std::endl;

		return true;
	}

/*	bool ThermoPolicy_CHEMKIN::WriteOnASCIIFile(std::stringstream &fOutput) const
	{
		// High temperature
		for (int i=4;i<=10;i++)
			fOutput << std::setprecision(8) << thermodynamics_parameters_[i] << " ";
		
		// Low temperature
		for (int i=11;i<=17;i++)
			fOutput << std::setprecision(8) << thermodynamics_parameters_[i] << " ";

		// Temperatures (L/H/M)
		fOutput << std::setprecision(8) << thermodynamics_parameters_[1] << " ";
		fOutput << std::setprecision(8) << thermodynamics_parameters_[3] << " ";
		fOutput << std::setprecision(8) << thermodynamics_parameters_[2] << " ";

		// Molecular weight
		fOutput << std::setprecision(12) << MolecularWeight() << std::endl;

		return true;
	}*/

	void ThermoPolicy_CHEMKIN::ThermodynamicsStatus(std::ostream& fOut, const std::vector<double>& T) const
	{
		const double mw = MolecularWeight();
		double conv_h	= mw/Conversions::J_from_kcal/1000.;
		double conv_s   = mw/Conversions::J_from_kcal;
		double conv_cp  = mw/Conversions::J_from_kcal;

		fOut << " ====================================================================================================" << std::endl;
		fOut << "   THERMO TABLE FOR MOLECULE " << name_thermo_                                                         << std::endl;
		fOut << " ====================================================================================================" << std::endl;
		fOut << "   Elemental Composition:" << std::endl;

		for(unsigned int i=0;i<atomic_composition_.element_names().size();i++)
			if (atomic_composition_.element_coefficients()[i] != 0.)	
				fOut << "     " << atomic_composition_.element_names()[i] << " " << atomic_composition_.element_coefficients()[i] << std::endl;

		fOut << "   Formation Enthalpy at 298K =          " << enthalpy(298.)*conv_h		<< " kcal/mol"	<< std::endl;
		fOut << "   Formation Free Gibbs Energy at 298K = " << gibbs_energy(298.)*conv_h	<< " kcal/mol"	<< std::endl;
		fOut << "   Molecular Weight =                    " << mw							<< " kg/kmol"	<< std::endl;
		fOut << " -----------------------------------------------------------------------------------------------------" << std::endl;
		fOut << "   Temperature       Cp             H             G            S             DH            DG         " << std::endl;          
		fOut << "       [K]       [cal/mol/K]    [kcal/mol]    [kcal/mol]  [cal/(mol-K)]  [kcal/mol]     [kcal/mol]    " << std::endl;        
		fOut << " -----------------------------------------------------------------------------------------------------" << std::endl;

		for(unsigned int i=0;i<T.size();i++)
		{
			fOut << "     " << std::setw(14) << std::left << std::setprecision(4) << std::fixed << T[i];
			fOut            << std::setw(14) << std::left << std::setprecision(4) << std::fixed << cp(T[i])*conv_cp;
			fOut            << std::setw(14) << std::left << std::setprecision(4) << std::fixed << enthalpy(T[i])*conv_h;
			fOut            << std::setw(14) << std::left << std::setprecision(4) << std::fixed << gibbs_energy(T[i])*conv_h;
			fOut            << std::setw(14) << std::left << std::setprecision(4) << std::fixed << entropy(T[i])*conv_s;
			fOut            << std::setw(14) << std::left << std::setprecision(4) << std::fixed << (enthalpy(T[i])-enthalpy(298.))*conv_h;
			fOut            << std::setw(14) << std::left << std::setprecision(4) << std::fixed << (gibbs_energy(T[i])-gibbs_energy(298.))*conv_h;
			fOut			<< std::endl;
		}

		fOut << " -----------------------------------------------------------------------------------------------------" << std::endl;
		fOut << std::endl << std::endl;
	}

	void ThermoPolicy_CHEMKIN::ThermodynamicsStatus(std::ostream& fOut) const
	{
		const double mw = MolecularWeight();
		double conv_h = mw / Conversions::J_from_kcal / 1000.;
		double conv_s = mw / Conversions::J_from_kcal;
		double conv_cp = mw / Conversions::J_from_kcal;

		// Species name
		fOut << std::setw(22) << std::left << std::setprecision(4) << std::fixed << name_thermo_;

		// Molecular weight
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << mw;

		// Thermo data @ 298K
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << cp(298.)*conv_cp;
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << enthalpy(298.)*conv_h;
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << gibbs_energy(298.)*conv_h;
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << entropy(298.)*conv_s;

		// Thermo data @ 1000K
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << cp(1000.)*conv_cp;
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << enthalpy(1000.)*conv_h;
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << gibbs_energy(1000.)*conv_h;
		fOut << std::setw(14) << std::right << std::setprecision(4) << std::fixed << entropy(1000.)*conv_s;

		fOut << std::endl;
	}

	int ThermoPolicy_CHEMKIN::CheckThermodynamicConsistency(std::ostream& fout)
	{
		double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

		double max_error = 1.e-3;
		double error_to_correct = 1.e-6;

		double T = thermodynamics_parameters_[2];

		double cp_plus   = ( thermodynamics_parameters_[4] + T*(thermodynamics_parameters_[5] + T*(thermodynamics_parameters_[6] + T*(thermodynamics_parameters_[7] + T*thermodynamics_parameters_[8]))) )*RGAS_over_MW_;
		double cp_minus  = ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13] + T*(thermodynamics_parameters_[14] + T*thermodynamics_parameters_[15]))) )*RGAS_over_MW_;

		double h_plus    = ( ( thermodynamics_parameters_[4]  + T*(thermodynamics_parameters_[5]/2. +  T*(thermodynamics_parameters_[6]/3.  + T*(thermodynamics_parameters_[7]/4.  + T*thermodynamics_parameters_[8]/5.))) ) + thermodynamics_parameters_[9]/T)  *T * RGAS_over_MW_;
		double h_minus   = ( ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12]/2. + T*(thermodynamics_parameters_[13]/3. + T*(thermodynamics_parameters_[14]/4. + T*thermodynamics_parameters_[15]/5.)))) + thermodynamics_parameters_[16]/T) *T * RGAS_over_MW_;

		double s_plus    = ( ( thermodynamics_parameters_[4]*std::log(T)  + T*(thermodynamics_parameters_[5] +  T*(thermodynamics_parameters_[6]/2.  + T*(thermodynamics_parameters_[7]/3.  + T*thermodynamics_parameters_[8]/4.))) ) + thermodynamics_parameters_[10]) * RGAS_over_MW_;
		double s_minus   = ( ( thermodynamics_parameters_[11]*std::log(T) + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13]/2. + T*(thermodynamics_parameters_[14]/3. + T*thermodynamics_parameters_[15]/4.)))) + thermodynamics_parameters_[17]) * RGAS_over_MW_;

		double dcp_plus   = (  thermodynamics_parameters_[5] +  T*(2.*thermodynamics_parameters_[6] +  T*(3.*thermodynamics_parameters_[7] +  4.*T*thermodynamics_parameters_[8])) )*RGAS_over_MW_;
		double dcp_minus  = (  thermodynamics_parameters_[12] + T*(2.*thermodynamics_parameters_[13] + T*(3.*thermodynamics_parameters_[14] + 4.*T*thermodynamics_parameters_[15])) )*RGAS_over_MW_;

		double error_cp  = std::fabs(cp_plus - cp_minus)/(0.50*(cp_plus+cp_minus));
		double error_h   = std::fabs(h_plus - h_minus)/(0.50*(h_plus+h_minus));
		double error_s   = std::fabs(s_plus - s_minus)/(0.50*(s_plus+s_minus));
		double error_dcp = std::fabs(dcp_plus - dcp_minus)/(0.50*(dcp_plus+dcp_minus));

		// Error: the thermodynamic properties are not consistent
		if (error_cp > max_error || error_h > max_error || error_s > max_error)
		{
			fout << "Species: " << name_thermo_ << std::endl;
			fout << "The thermodynamic properties are not consistent. Please check the thermodynamic coefficients in order to remove any inconsistencies." << std::endl;
			fout << "Thermodynamic function        (-)           (+)          err(%)" << std::endl;
			fout << "Cp [J/kg/K]             " << cp_minus  << "  " << cp_plus  << "  " << error_cp*100. << std::endl;
			fout << "H  [J/kg]               " << h_minus   << "  " << h_plus   << "  " << error_h*100. << std::endl;
			fout << "S  [J/kg/K]             " << s_minus   << "  " << s_plus   << "  " << error_s*100. << std::endl;
			fout << "dCp/dT [J/kg/K^2]       " << dcp_minus << "  " << dcp_plus << "  " << error_dcp*100. << std::endl;
			fout << std::endl;
			return -1;
		}

		// The coefficients are corrected to make the thermodynamic properties at the intermediate temperature perfectly equal!
		// These corrections are usually very small and are not reported to the user
		// The corrections are performed only on the low temperature coefficients
		{
			thermodynamics_parameters_[11] = cp_plus/RGAS_over_MW_ - (cp_minus/RGAS_over_MW_ - thermodynamics_parameters_[11]);

			double h_minus   = ( ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12]/2. + T*(thermodynamics_parameters_[13]/3. + T*(thermodynamics_parameters_[14]/4. + T*thermodynamics_parameters_[15]/5.)))) + thermodynamics_parameters_[16]/T) *T * RGAS_over_MW_;
			thermodynamics_parameters_[16] = h_plus/RGAS_over_MW_ - (h_minus/RGAS_over_MW_ - thermodynamics_parameters_[16]);

			double s_minus   = ( ( thermodynamics_parameters_[11]*std::log(T) + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13]/2. + T*(thermodynamics_parameters_[14]/3. + T*thermodynamics_parameters_[15]/4.)))) + thermodynamics_parameters_[17]) * RGAS_over_MW_;
			thermodynamics_parameters_[17] = s_plus/RGAS_over_MW_ - (s_minus/RGAS_over_MW_ - thermodynamics_parameters_[17]);
		}

		// Only if the error is larger than error_to_correct the information is returned to the user
		if (error_cp > error_to_correct || error_h > error_to_correct || error_s > error_to_correct)
			return 0;
		else
			return 1;
	}

	int ThermoPolicy_CHEMKIN::CheckForAnomalies(std::ostream& fout)
	{
		// Low temperature Cp
		{
			unsigned int n_real;
			double x[3];
			FindCubicRoots(4.*thermodynamics_parameters_[15], 3.*thermodynamics_parameters_[14], 2.*thermodynamics_parameters_[13], thermodynamics_parameters_[12], n_real, x);
			
			bool local_maxima = false;
			for(unsigned int i=0;i<n_real;i++)
				if( (x[i] > thermodynamics_parameters_[1]) && (x[i] < thermodynamics_parameters_[2]) )
					local_maxima = true;

			if (local_maxima == true)
			{
				fout << "LT Cp for " << name_thermo_ << " defined in [";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[1];
				fout << "-";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[2];
				fout << "] has local maxima at: ";
				for(unsigned int i=0;i<n_real;i++)
					if ( (x[i] > thermodynamics_parameters_[1]) && (x[i] < thermodynamics_parameters_[2]) )fout << x[i] << " ";
				fout << std::endl;
			}
		}

		// High temperature Cp
		{
			unsigned int n_real;
			double x[3];
			FindCubicRoots(4.*thermodynamics_parameters_[8], 3.*thermodynamics_parameters_[7], 2.*thermodynamics_parameters_[6], thermodynamics_parameters_[5], n_real, x);
			
			bool local_maxima = false;
			for(unsigned int i=0;i<n_real;i++)
				if( (x[i] > thermodynamics_parameters_[2]) && (x[i] < thermodynamics_parameters_[3]) )
					local_maxima = true;

			if (local_maxima == true)
			{
				fout << "HT Cp for " << name_thermo_ << " defined in [";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[2];
				fout << "-";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[3];
				fout << "] has local maxima at: ";
				for(unsigned int i=0;i<n_real;i++)
					if( (x[i] > thermodynamics_parameters_[2]) && (x[i] < thermodynamics_parameters_[3]) ) fout << x[i] << " ";
				fout << std::endl;
			}
		}

		// Low temperature H
		{
			unsigned int n_real;
			double x[4];
			FindBiquadricRoots(thermodynamics_parameters_[15], thermodynamics_parameters_[14], thermodynamics_parameters_[13], thermodynamics_parameters_[12], thermodynamics_parameters_[11], n_real, x);

			bool local_maxima = false;
			for(unsigned int i=0;i<n_real;i++)
				if( (x[i] > thermodynamics_parameters_[1]) && (x[i] < thermodynamics_parameters_[2]) )
					local_maxima = true;

			if (local_maxima == true)
			{
				fout << "LT enthalpy for " << name_thermo_ << " defined in [";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[1];
				fout << "-";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[2];
				fout << "] has local maxima at: ";
				for(unsigned int i=0;i<n_real;i++)
					if( (x[i] > thermodynamics_parameters_[1]) && (x[i] < thermodynamics_parameters_[2]) ) fout << x[i] << " ";
				fout << std::endl;
			}
		}

		// High temperature H
		{
			unsigned int n_real;
			double x[4];
			FindBiquadricRoots(thermodynamics_parameters_[8], thermodynamics_parameters_[7], thermodynamics_parameters_[6], thermodynamics_parameters_[5], thermodynamics_parameters_[4], n_real, x);

			bool local_maxima = false;
			for(unsigned int i=0;i<n_real;i++)
				if( (x[i] > thermodynamics_parameters_[2]) && (x[i] < thermodynamics_parameters_[3]) )
					local_maxima = true;

			if (local_maxima == true)
			{
				fout << "HT enthalpy for " << name_thermo_ << " defined in [";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[2];
				fout << "-";
				fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[3];
				fout << "] has local maxima at: ";
				for(unsigned int i=0;i<n_real;i++)
					if( (x[i] > thermodynamics_parameters_[2]) && (x[i] < thermodynamics_parameters_[3]) ) fout << x[i] << " ";
				fout << std::endl;
			}
		}

		// High temperature S
		{
			double T = thermodynamics_parameters_[1];
			double s_previous = ( ( thermodynamics_parameters_[11]*std::log(T) + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13]/2. + T*(thermodynamics_parameters_[14]/3. + T*thermodynamics_parameters_[15]/4.)))) + thermodynamics_parameters_[17]);

			for (;;)
			{
				T += 1.;
				if (T>thermodynamics_parameters_[2])
					break;

				const double s = ( ( thermodynamics_parameters_[11]*std::log(T) + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13]/2. + T*(thermodynamics_parameters_[14]/3. + T*thermodynamics_parameters_[15]/4.)))) + thermodynamics_parameters_[17]);
				if (s < s_previous*(1.+1e-32))
				{
					fout << "LT entropy for " << name_thermo_ << " defined in [";
					fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[1];
					fout << "-";
					fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[2];
					fout << "] has local maxima at: ";
					fout << T << " ";
					fout << std::endl;
					break;
				}
				s_previous = s;
			}
		}

		// High temperature S
		{
			double T = thermodynamics_parameters_[2];
			double s_previous = ( ( thermodynamics_parameters_[4]*std::log(T)  + T*(thermodynamics_parameters_[5] +  T*(thermodynamics_parameters_[6]/2.  + T*(thermodynamics_parameters_[7]/3.  + T*thermodynamics_parameters_[8]/4.))) ) + thermodynamics_parameters_[10]);

			for (;;)
			{
				T += 1.;
				if (T>thermodynamics_parameters_[3])
					break;

				const double s = ( ( thermodynamics_parameters_[4]*std::log(T)  + T*(thermodynamics_parameters_[5] +  T*(thermodynamics_parameters_[6]/2.  + T*(thermodynamics_parameters_[7]/3.  + T*thermodynamics_parameters_[8]/4.))) ) + thermodynamics_parameters_[10]);
				if (s < s_previous*(1.+1e-32))
				{
					fout << "HT entropy for " << name_thermo_ << " defined in [";
					fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[2];
					fout << "-";
					fout << std::fixed << std::setprecision(2) << thermodynamics_parameters_[3];
					fout << "] has local maxima at: ";
					fout << T << " ";
					fout << std::endl;
					break;
				}
				s_previous = s;
			}
		}

		return 0;
	}

	int ThermoPolicy_CHEMKIN::ReportStatusSpecificHeatCoefficient(std::ostream& fout)
	{
		const double T = thermodynamics_parameters_[2];
		const double cp_plus   = ( thermodynamics_parameters_[4] + T*(thermodynamics_parameters_[5] + T*(thermodynamics_parameters_[6] + T*(thermodynamics_parameters_[7] + T*thermodynamics_parameters_[8]))) );
		const double cp_minus  = ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13] + T*(thermodynamics_parameters_[14] + T*thermodynamics_parameters_[15]))) );

		const double dcp_plus   = (  thermodynamics_parameters_[5] +  T*(2.*thermodynamics_parameters_[6] +  T*(3.*thermodynamics_parameters_[7] +  4.*T*thermodynamics_parameters_[8])) );
		const double dcp_minus  = (  thermodynamics_parameters_[12] + T*(2.*thermodynamics_parameters_[13] + T*(3.*thermodynamics_parameters_[14] + 4.*T*thermodynamics_parameters_[15])) );

		const double dcp2_plus   = (  2.*thermodynamics_parameters_[6] +  T*(6.*thermodynamics_parameters_[7] +  12.*T*thermodynamics_parameters_[8]) );
		const double dcp2_minus  = (  2.*thermodynamics_parameters_[13] + T*(6.*thermodynamics_parameters_[14] + 12.*T*thermodynamics_parameters_[15]) );

		const double dcp3_plus   = (  6.*thermodynamics_parameters_[7] +  24.*T*thermodynamics_parameters_[8] );
		const double dcp3_minus  = (  6.*thermodynamics_parameters_[14] + 24.*T*thermodynamics_parameters_[15] );

		const double dcp4_plus   = 24.*thermodynamics_parameters_[8];
		const double dcp4_minus  = 24.*thermodynamics_parameters_[15];

		const double error_cp  = std::fabs(cp_plus - cp_minus)/(0.50*(cp_plus+cp_minus));
		const double error_dcp = std::fabs(dcp_plus - dcp_minus)/std::fabs(0.50*(dcp_plus+dcp_minus)+1e-32);
		const double error_dcp2 = std::fabs(dcp2_plus - dcp2_minus)/std::fabs(0.50*(dcp2_plus+dcp2_minus)+1e-32);
		const double error_dcp3 = std::fabs(dcp3_plus - dcp3_minus)/std::fabs(0.50*(dcp3_plus+dcp3_minus)+1e-32);
		const double error_dcp4 = std::fabs(dcp4_plus - dcp4_minus)/std::fabs(0.50*(dcp4_plus+dcp4_minus)+1e-32);

		fout << std::setw(18) << std::left << name_thermo_;
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[1];
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[2];
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[3];

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << cp_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << cp_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_cp*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dcp*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp2_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp2_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dcp2*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp3_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp3_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dcp3*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp4_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dcp4_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dcp4*100.;

		fout << std::endl;

		return 0;
	}

	int ThermoPolicy_CHEMKIN::ReportStatusEnthalpy(std::ostream& fout)
	{
		const double T = thermodynamics_parameters_[2];
		const double h_plus    = ( ( thermodynamics_parameters_[4]  + T*(thermodynamics_parameters_[5]/2. +  T*(thermodynamics_parameters_[6]/3.  + T*(thermodynamics_parameters_[7]/4.  + T*thermodynamics_parameters_[8]/5.))) ) + thermodynamics_parameters_[9]/T) ;
		const double h_minus   = ( ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12]/2. + T*(thermodynamics_parameters_[13]/3. + T*(thermodynamics_parameters_[14]/4. + T*thermodynamics_parameters_[15]/5.)))) + thermodynamics_parameters_[16]/T) ;

		const double dh_plus   = ( thermodynamics_parameters_[4] + T*(thermodynamics_parameters_[5] + T*(thermodynamics_parameters_[6] + T*(thermodynamics_parameters_[7] + T*thermodynamics_parameters_[8]))) ) / T;
		const double dh_minus  = ( thermodynamics_parameters_[11] + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13] + T*(thermodynamics_parameters_[14] + T*thermodynamics_parameters_[15]))) ) /T;

		const double dh2_plus   = (  thermodynamics_parameters_[5] +  T*(2.*thermodynamics_parameters_[6] +  T*(3.*thermodynamics_parameters_[7] +  4.*T*thermodynamics_parameters_[8])) ) / T;
		const double dh2_minus  = (  thermodynamics_parameters_[12] + T*(2.*thermodynamics_parameters_[13] + T*(3.*thermodynamics_parameters_[14] + 4.*T*thermodynamics_parameters_[15])) ) / T;

		const double dh3_plus   = (  2.*thermodynamics_parameters_[6] +  T*(6.*thermodynamics_parameters_[7] +  12.*T*thermodynamics_parameters_[8]) ) / T;
		const double dh3_minus  = (  2.*thermodynamics_parameters_[13] + T*(6.*thermodynamics_parameters_[14] + 12.*T*thermodynamics_parameters_[15]) ) / T;

		const double dh4_plus   = (  6.*thermodynamics_parameters_[7] +  24.*T*thermodynamics_parameters_[8] ) / T;
		const double dh4_minus  = (  6.*thermodynamics_parameters_[14] + 24.*T*thermodynamics_parameters_[15] ) / T;

		const double error_h  = std::fabs(h_plus - h_minus)/(0.50*(h_plus+h_minus));
		const double error_dh = std::fabs(dh_plus - dh_minus)/std::fabs(0.50*(dh_plus+dh_minus)+1e-32);
		const double error_dh2 = std::fabs(dh2_plus - dh2_minus)/std::fabs(0.50*(dh2_plus+dh2_minus)+1e-32);
		const double error_dh3 = std::fabs(dh3_plus - dh3_minus)/std::fabs(0.50*(dh3_plus+dh3_minus)+1e-32);
		const double error_dh4 = std::fabs(dh4_plus - dh4_minus)/std::fabs(0.50*(dh4_plus+dh4_minus)+1e-32);

		fout << std::setw(18) << std::left << name_thermo_;
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[1];
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[2];
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[3];

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << h_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << h_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_h*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dh*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh2_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh2_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dh2*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh3_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh3_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dh3*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh4_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << dh4_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_dh4*100.;

		fout << std::endl;

		return 0;
	}

	int ThermoPolicy_CHEMKIN::ReportStatusEntropy(std::ostream& fout)
	{
		const double T = thermodynamics_parameters_[2];
		const double s_plus    = ( ( thermodynamics_parameters_[4]*std::log(T)  + T*(thermodynamics_parameters_[5] +  T*(thermodynamics_parameters_[6]/2.  + T*(thermodynamics_parameters_[7]/3.  + T*thermodynamics_parameters_[8]/4.))) ) + thermodynamics_parameters_[10]);
		const double s_minus   = ( ( thermodynamics_parameters_[11]*std::log(T) + T*(thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13]/2. + T*(thermodynamics_parameters_[14]/3. + T*thermodynamics_parameters_[15]/4.)))) + thermodynamics_parameters_[17]);

		const double ds_plus    =  thermodynamics_parameters_[4]/T  + thermodynamics_parameters_[5] +  T*(thermodynamics_parameters_[6]  + T*(thermodynamics_parameters_[7]  + T*thermodynamics_parameters_[8] ) );
		const double ds_minus   =  thermodynamics_parameters_[11]/T + thermodynamics_parameters_[12] + T*(thermodynamics_parameters_[13] + T*(thermodynamics_parameters_[14] + T*thermodynamics_parameters_[15] ) ) ;

		const double ds2_plus    =  -thermodynamics_parameters_[4]/T/T  + thermodynamics_parameters_[6]  + T*(2.*thermodynamics_parameters_[7]  + 3.*T*thermodynamics_parameters_[8] ) ;
		const double ds2_minus   =  -thermodynamics_parameters_[11]/T/T + thermodynamics_parameters_[13] + T*(2.*thermodynamics_parameters_[14] + 3.*T*thermodynamics_parameters_[15] );

		const double ds3_plus    =  2.*thermodynamics_parameters_[4]/T/T/T  + 2.*thermodynamics_parameters_[7]  + 6.*T*thermodynamics_parameters_[8];
		const double ds3_minus   =  2.*thermodynamics_parameters_[11]/T/T/T + 2.*thermodynamics_parameters_[14] + 6.*T*thermodynamics_parameters_[15];

		const double ds4_plus    =  -6.*thermodynamics_parameters_[4]/T/T/T/T  + 6.*thermodynamics_parameters_[8];
		const double ds4_minus   =  -6.*thermodynamics_parameters_[11]/T/T/T/T + 6.*thermodynamics_parameters_[15];

		const double error_s  = std::fabs(s_plus - s_minus)/(0.50*(s_plus+s_minus));
		const double error_ds = std::fabs(ds_plus - ds_minus)/std::fabs(0.50*(ds_plus+ds_minus)+1e-32);
		const double error_ds2 = std::fabs(ds2_plus - ds2_minus)/std::fabs(0.50*(ds2_plus+ds2_minus)+1e-32);
		const double error_ds3 = std::fabs(ds3_plus - ds3_minus)/std::fabs(0.50*(ds3_plus+ds3_minus)+1e-32);
		const double error_ds4 = std::fabs(ds4_plus - ds4_minus)/std::fabs(0.50*(ds4_plus+ds4_minus)+1e-32);

		fout << std::setw(18) << std::left << name_thermo_;
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[1];
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[2];
		fout << std::setw(9) << std::fixed << std::setprecision(2) << std::left << thermodynamics_parameters_[3];

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << s_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << s_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_s*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_ds*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds2_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds2_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_ds2*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds3_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds3_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_ds3*100.;

		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds4_minus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << ds4_plus;
		fout << std::setw(18) << std::scientific << std::setprecision(8) << std::left << error_ds4*100.;

		fout << std::endl;

		return 0;
	}

	
	void ThermoPolicy_CHEMKIN::ReformulationOfThermodynamics(std::ostream& fout, const unsigned int policy_intermediate_temperature, const double intermediate_temperature, const double max_temperature)
	{
		unsigned int n = 15;
		OpenSMOKEVectorDouble temperatures(2*n);
		temperatures[1] = thermodynamics_parameters_[1];
		for(unsigned int i=2;i<=n;i++)
			temperatures[i] = temperatures[i-1] + (thermodynamics_parameters_[2]-thermodynamics_parameters_[1])/double(n-1);
		for(unsigned int i=n+1;i<=2*n;i++)
			temperatures[i] = temperatures[i-1] + (max_temperature-thermodynamics_parameters_[2])/double(n);

		// Evaluating the profile to fit
		Eigen::MatrixXd y(temperatures.Size(), 1);
		for(int i=1;i<=temperatures.Size();i++)
		{
			const double T = temperatures[i];
			const double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

			y(i-1,0) = cp(temperatures[i]) / RGAS_over_MW_;
		}

		int N = 1;
		Eigen::MatrixXd XTX(5+N,5+N);
		Eigen::MatrixXd XT(5+N, temperatures.Size());
		Eigen::MatrixXd X(temperatures.Size(), 5+N);
		Eigen::MatrixXd Y(5+N, 1);

		// Assembling X and XT Matrices (Fixed part)
		for(int i=0;i<temperatures.Size();i++)
		{
			const double T = temperatures[i+1];
			X(i,0) = 1.;
			X(i,1) = T;
			X(i,2) = T*T;
			X(i,3) = T*T*T;
			X(i,4) = T*T*T*T;
		}

		double Tknot = 700.;
		double T_min_relative_error = Tknot;
		double min_relative_error = 1.e16;

		double ALT, BLT, CLT, DLT, ELT;
		double AHT, BHT, CHT, DHT, EHT;

		for(;;)
		{
			if (policy_intermediate_temperature == 0)			// user-defined intermediate temperature
				Tknot = intermediate_temperature;

			for(int i=0;i<temperatures.Size();i++)
			{
				const double T = temperatures[i+1];

				if (T<=Tknot)
				{
					X(i,5) = 0.;
					if (N>1)	X(i,6) = 0.;
					if (N>2)	X(i,7) = 0.;
				}
				else
				{
					double difference = T-Tknot;
					X(i,5) = difference*difference*difference*difference;
					if (N>1)	X(i,6) = difference*difference*difference;
					if (N>2)	X(i,7) = difference*difference;
				}
			}

			// Calculating the LS matrix
			XT = X.transpose();
			XTX = XT*X;
			Y = XT*y;

			Eigen::MatrixXd fittedParameters;
			fittedParameters = XTX.partialPivLu().solve(Y);

			const double Tknot2 = Tknot*Tknot;
			const double Tknot3 = Tknot2*Tknot;
			const double Tknot4 = Tknot2*Tknot2;

			ALT = fittedParameters(0,0);
			BLT = fittedParameters(1,0);
			CLT = fittedParameters(2,0);
			DLT = fittedParameters(3,0);
			ELT = fittedParameters(4,0);

			AHT = ALT + fittedParameters(5,0)*Tknot4;
			BHT = BLT - 4.*fittedParameters(5,0)*Tknot3;
			CHT = CLT + 6.*fittedParameters(5,0)*Tknot2;
			DHT = DLT - 4.*fittedParameters(5,0)*Tknot;
			EHT = ELT + fittedParameters(5,0);

			if (N>1)
			{
				AHT +=  - fittedParameters(6,0)*Tknot3;
				BHT +=  3.*fittedParameters(6,0)*Tknot2;
				CHT += -3.*fittedParameters(6,0)*Tknot;
				DHT += fittedParameters(6,0);
			}
			if (N>2)
			{
				AHT += fittedParameters(7,0)*Tknot2;
				BHT += -2.*fittedParameters(7,0)*Tknot;
				CHT += fittedParameters(7,0);
			}

			double relative_error = 0.;
			for(int i=1;i<=temperatures.Size();i++)
			{
				const double T = temperatures[i];
				double cp_new;

				if (T<=Tknot)
					cp_new  = ( ALT + T*(BLT + T*(CLT + T*(DLT + T*ELT))) );
				else
					cp_new  = ( AHT + T*(BHT + T*(CHT + T*(DHT + T*EHT))) );

				relative_error += std::fabs(cp_new-y(i-1,0))/y(i-1,0);

			//	double h_new, s_new;
			//	if (T<=Tknot)
			//	{
			//		h_new  = ( ALT + T*(BLT/2. + T*(CLT/3. + T*(DLT/4. + T*ELT/5.))) + FLT/T);
			//		s_new  = ( ALT*std::log(T) + T*(BLT + T*(CLT/2. + T*(DLT/3. + T*ELT/4.))) + GLT);
			//	}
			//	else
			//	{
			//		h_new  = ( AHT + T*(BHT/2. + T*(CHT/3. + T*(DHT/4. + T*EHT/5.))) + FHT/T);
			//		s_new  = ( AHT*std::log(T) + T*(BHT + T*(CHT/2. + T*(DHT/3. + T*EHT/4.))) + GHT);
			//	}

			}

			if (policy_intermediate_temperature == 0)
				break;
				
			if (relative_error < min_relative_error)
			{
				min_relative_error = relative_error;
				T_min_relative_error = Tknot;
			}
			
			Tknot += 10.;
			if (Tknot > 1800.)
			{
				ReformulationOfThermodynamics(fout, 0, T_min_relative_error, max_temperature);
				return;
			}
		}

		{
			const double RGAS_over_MW_=PhysicalConstants::R_J_mol/(MolecularWeight()*1.e-3);

			const double Tknot2 = Tknot*Tknot;
			const double Tknot3 = Tknot*Tknot2;
			const double Tknot4 = Tknot2*Tknot2;
			const double FLT = ( enthalpy(Tknot)/RGAS_over_MW_/Tknot - (ALT+BLT/2.*Tknot+CLT/3.*Tknot2+DLT/4.*Tknot3+ELT/5.*Tknot4) ) *Tknot;
			const double FHT = ( enthalpy(Tknot)/RGAS_over_MW_/Tknot - (AHT+BHT/2.*Tknot+CHT/3.*Tknot2+DHT/4.*Tknot3+EHT/5.*Tknot4) ) *Tknot;

			double GLT = ( entropy(Tknot)/RGAS_over_MW_ - (ALT*std::log(Tknot)+BLT*Tknot+CLT/2.*Tknot2+DLT/3.*Tknot3+ELT/4.*Tknot4) );
			double GHT = ( entropy(Tknot)/RGAS_over_MW_ - (AHT*std::log(Tknot)+BHT*Tknot+CHT/2.*Tknot2+DHT/3.*Tknot3+EHT/4.*Tknot4) );

			// Write on a file
			fout << std::setw(16) << std::left << name_thermo_;
			fout << std::setw(8)  << std::left << " ";

			std::vector<unsigned int> list_of_nonstandard_atomic_composition;
			std::vector<unsigned int> list_of_standard_atomic_composition;

			for(unsigned int i=0;i<atomic_composition_.element_coefficients().size();i++)
			{
				if (atomic_composition_.element_coefficients()[i] != 0.)	
				{
					if ( IsInteger(atomic_composition_.element_coefficients()[i]) == false || int(atomic_composition_.element_coefficients()[i]) > 99)
						list_of_nonstandard_atomic_composition.push_back(i);
					else
						list_of_standard_atomic_composition.push_back(i);
				}
			}

			if (list_of_standard_atomic_composition.size()<=4 && list_of_nonstandard_atomic_composition.size() == 0)
			{
				for(unsigned int j=0;j<list_of_standard_atomic_composition.size();j++)
				{
					unsigned int i = list_of_standard_atomic_composition[j];
					fout << std::setw(2)  << std::left  << atomic_composition_.element_names()[i];
					fout << std::setw(3)  << std::right << int(atomic_composition_.element_coefficients()[i]);
				}
				for(std::size_t j=list_of_standard_atomic_composition.size();j<4;j++)
				{
						fout << std::setw(2)  << " ";
						fout << std::setw(3)  << " ";
				}
			}
			else
			{
				list_of_nonstandard_atomic_composition.insert(list_of_nonstandard_atomic_composition.end(), list_of_standard_atomic_composition.begin(), list_of_standard_atomic_composition.end());
				for(unsigned int j=0;j<4;j++)
				{
						fout << std::setw(2)  << " ";
						fout << std::setw(3)  << " ";
				}
			}
		
			fout << phase_;
			fout << std::setw(10)  << std::right << std::fixed << std::setprecision(2) << thermodynamics_parameters_[1];
			fout << std::setw(10)  << std::right << std::fixed << std::setprecision(2) << max_temperature;
			fout << std::setw(8)   << std::right << std::fixed << std::setprecision(2) << Tknot;
			fout << "      1";
			if (list_of_nonstandard_atomic_composition.size() != 0)
				fout << "&";
			fout << std::endl;

			if (list_of_nonstandard_atomic_composition.size() != 0)
			{
				for (unsigned int j = 0; j<list_of_nonstandard_atomic_composition.size(); j++)
				{
					unsigned int i = list_of_nonstandard_atomic_composition[j];
					fout << atomic_composition_.element_names()[i] << " " << atomic_composition_.element_coefficients()[i] << " ";
				}
				fout << std::endl;
			}

			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(AHT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(BHT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(CHT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(DHT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(EHT,2);
			fout << "    2" << std::endl;
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(FHT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(GHT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(ALT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(BLT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(CLT,2);
			fout << "    3" << std::endl;
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(DLT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(ELT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(FLT,2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(GLT,2);
			fout << std::setw(15) << " ";
			fout << "    4" << std::endl;
		}
	}

	void ThermoPolicy_CHEMKIN::ReformulationOfThermodynamicsFixedIntermediateTemperature(std::ostream& fout, const unsigned int policy_intermediate_temperature, const double intermediate_temperature, const double max_temperature)
	{
		const double Tknot = 1000.;

		unsigned int n = 15;
		OpenSMOKEVectorDouble temperatures(2 * n);
		temperatures[1] = thermodynamics_parameters_[1];
		for (unsigned int i = 2; i <= n; i++)
			temperatures[i] = temperatures[i - 1] + (thermodynamics_parameters_[2] - thermodynamics_parameters_[1]) / double(n - 1);
		for (unsigned int i = n + 1; i <= 2 * n; i++)
			temperatures[i] = temperatures[i - 1] + (max_temperature - thermodynamics_parameters_[2]) / double(n);

		// Evaluating the profile to fit
		Eigen::MatrixXd y(temperatures.Size(), 1);
		for (int i = 1; i <= temperatures.Size(); i++)
		{
			const double T = temperatures[i];
			const double RGAS_over_MW_ = PhysicalConstants::R_J_mol / (MolecularWeight()*1.e-3);

			y(i - 1, 0) = cp(temperatures[i]) / RGAS_over_MW_;
		}

		int N = 1;
		Eigen::MatrixXd XTX(5 + N, 5 + N);
		Eigen::MatrixXd XT(5 + N, temperatures.Size());
		Eigen::MatrixXd X(temperatures.Size(), 5 + N);
		Eigen::MatrixXd Y(5 + N, 1);

		// Assembling X and XT Matrices (Fixed part)
		for (int i = 0; i<temperatures.Size(); i++)
		{
			const double T = temperatures[i + 1];
			X(i, 0) = 1.;
			X(i, 1) = T;
			X(i, 2) = T*T;
			X(i, 3) = T*T*T;
			X(i, 4) = T*T*T*T;
		}

		double ALT, BLT, CLT, DLT, ELT;
		double AHT, BHT, CHT, DHT, EHT;

		{
			for (int i = 0; i<temperatures.Size(); i++)
			{
				const double T = temperatures[i + 1];

				if (T <= Tknot)
				{
					X(i, 5) = 0.;
					if (N>1)	X(i, 6) = 0.;
					if (N>2)	X(i, 7) = 0.;
				}
				else
				{
					double difference = T - Tknot;
					X(i, 5) = difference*difference*difference*difference;
					if (N>1)	X(i, 6) = difference*difference*difference;
					if (N>2)	X(i, 7) = difference*difference;
				}
			}

			// Calculating the LS matrix
			XT = X.transpose();
			XTX = XT*X;
			Y = XT*y;

			Eigen::MatrixXd fittedParameters;
			fittedParameters = XTX.partialPivLu().solve(Y);

			const double Tknot2 = Tknot*Tknot;
			const double Tknot3 = Tknot2*Tknot;
			const double Tknot4 = Tknot2*Tknot2;

			ALT = fittedParameters(0, 0);
			BLT = fittedParameters(1, 0);
			CLT = fittedParameters(2, 0);
			DLT = fittedParameters(3, 0);
			ELT = fittedParameters(4, 0);

			AHT = ALT + fittedParameters(5, 0)*Tknot4;
			BHT = BLT - 4.*fittedParameters(5, 0)*Tknot3;
			CHT = CLT + 6.*fittedParameters(5, 0)*Tknot2;
			DHT = DLT - 4.*fittedParameters(5, 0)*Tknot;
			EHT = ELT + fittedParameters(5, 0);

			if (N>1)
			{
				AHT += -fittedParameters(6, 0)*Tknot3;
				BHT += 3.*fittedParameters(6, 0)*Tknot2;
				CHT += -3.*fittedParameters(6, 0)*Tknot;
				DHT += fittedParameters(6, 0);
			}
			if (N>2)
			{
				AHT += fittedParameters(7, 0)*Tknot2;
				BHT += -2.*fittedParameters(7, 0)*Tknot;
				CHT += fittedParameters(7, 0);
			}

			double relative_error = 0.;
			for (int i = 1; i <= temperatures.Size(); i++)
			{
				const double T = temperatures[i];
				double cp_new;

				if (T <= Tknot)
					cp_new = (ALT + T*(BLT + T*(CLT + T*(DLT + T*ELT))));
				else
					cp_new = (AHT + T*(BHT + T*(CHT + T*(DHT + T*EHT))));
			}
		}

		{
			const double RGAS_over_MW_ = PhysicalConstants::R_J_mol / (MolecularWeight()*1.e-3);

			const double Tknot2 = Tknot*Tknot;
			const double Tknot3 = Tknot*Tknot2;
			const double Tknot4 = Tknot2*Tknot2;
			const double FLT = (enthalpy(Tknot) / RGAS_over_MW_ / Tknot - (ALT + BLT / 2.*Tknot + CLT / 3.*Tknot2 + DLT / 4.*Tknot3 + ELT / 5.*Tknot4)) *Tknot;
			const double FHT = (enthalpy(Tknot) / RGAS_over_MW_ / Tknot - (AHT + BHT / 2.*Tknot + CHT / 3.*Tknot2 + DHT / 4.*Tknot3 + EHT / 5.*Tknot4)) *Tknot;

			double GLT = (entropy(Tknot) / RGAS_over_MW_ - (ALT*std::log(Tknot) + BLT*Tknot + CLT / 2.*Tknot2 + DLT / 3.*Tknot3 + ELT / 4.*Tknot4));
			double GHT = (entropy(Tknot) / RGAS_over_MW_ - (AHT*std::log(Tknot) + BHT*Tknot + CHT / 2.*Tknot2 + DHT / 3.*Tknot3 + EHT / 4.*Tknot4));

			// Write on a file
			fout << std::setw(16) << std::left << name_thermo_;
			fout << std::setw(8) << std::left << " ";

			std::vector<unsigned int> list_of_nonstandard_atomic_composition;
			std::vector<unsigned int> list_of_standard_atomic_composition;

			for (unsigned int i = 0; i<atomic_composition_.element_coefficients().size(); i++)
			{
				if (atomic_composition_.element_coefficients()[i] != 0.)
				{
					if (IsInteger(atomic_composition_.element_coefficients()[i]) == false || int(atomic_composition_.element_coefficients()[i]) > 99)
						list_of_nonstandard_atomic_composition.push_back(i);
					else
						list_of_standard_atomic_composition.push_back(i);
				}
			}

			if (list_of_standard_atomic_composition.size() <= 4 && list_of_nonstandard_atomic_composition.size() == 0)
			{
				for (unsigned int j = 0; j<list_of_standard_atomic_composition.size(); j++)
				{
					unsigned int i = list_of_standard_atomic_composition[j];
					fout << std::setw(2) << std::left << atomic_composition_.element_names()[i];
					fout << std::setw(3) << std::right << int(atomic_composition_.element_coefficients()[i]);
				}
				for (std::size_t j = list_of_standard_atomic_composition.size(); j<4; j++)
				{
					fout << std::setw(2) << " ";
					fout << std::setw(3) << " ";
				}
			}
			else
			{
				list_of_nonstandard_atomic_composition.insert(list_of_nonstandard_atomic_composition.end(), list_of_standard_atomic_composition.begin(), list_of_standard_atomic_composition.end());
				for (unsigned int j = 0; j<4; j++)
				{
					fout << std::setw(2) << " ";
					fout << std::setw(3) << " ";
				}
			}

			fout << phase_;
			fout << std::setw(10) << std::right << std::fixed << std::setprecision(2) << thermodynamics_parameters_[1];
			fout << std::setw(10) << std::right << std::fixed << std::setprecision(2) << max_temperature;
			fout << std::setw(8) << std::right << std::fixed << std::setprecision(2) << Tknot;
			fout << "      1";
			if (list_of_nonstandard_atomic_composition.size() != 0)
				fout << "&";
			fout << std::endl;

			if (list_of_nonstandard_atomic_composition.size() != 0)
			{
				for (unsigned int j = 0; j<list_of_nonstandard_atomic_composition.size(); j++)
				{
					unsigned int i = list_of_nonstandard_atomic_composition[j];
					fout << atomic_composition_.element_names()[i] << " " << atomic_composition_.element_coefficients()[i] << " ";
				}
				fout << std::endl;
			}

			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(AHT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(BHT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(CHT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(DHT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(EHT, 2);
			fout << "    2" << std::endl;
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(FHT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(GHT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(ALT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(BLT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(CLT, 2);
			fout << "    3" << std::endl;
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(DLT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(ELT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(FLT, 2);
			fout << std::setw(15) << std::right << ScientificNotationWithFixedExponentDigits(GLT, 2);
			fout << std::setw(15) << " ";
			fout << "    4" << std::endl;
		}
	}
}



void FindBiquadricRoots(const double A, const double B, const double C, const double D, const double E, unsigned int& n, double x[4])
{
	x[0] = x[1] = x[2] = 0.;
	n = 0;

	if (A != 0)		// Biquadric equation
	{
		double p[5], r[3][5];
		p[0] = A; 
		p[1] = B; 
		p[2] = C; 
		p[3] = D;
		p[4] = E;

		BIQUADROOTS(p, r);

		if (r[2][1] == 0.)	x[n++] = r[1][1];
		if (r[2][2] == 0.)	x[n++] = r[1][2];
		if (r[2][3] == 0.)	x[n++] = r[1][3];
		if (r[2][4] == 0.)	x[n++] = r[1][4];
	}
	else
	{
		double x_cubic[3];
		FindCubicRoots(B, C, D, E, n, x_cubic);
		x[0] = x_cubic[0]; x[1] = x_cubic[1]; x[2] = x_cubic[2];
	}
}

void FindCubicRoots(const double A, const double B, const double C, const double D, unsigned int& n, double x[3])
{
	x[0] = x[1] = x[2] = 0.;
	n = 0;

	if (A != 0)		// Cubic equation
	{
		double p[5], r[3][5];
		p[0] = A; 
		p[1] = B; 
		p[2] = C; 
		p[3] = D;

		CUBICROOTS(p, r);

		if (r[2][1] == 0.)	x[n++] = r[1][1];
		if (r[2][2] == 0.)	x[n++] = r[1][2];
		if (r[2][3] == 0.)	x[n++] = r[1][3];
	}
	else if (B != 0.)	// Quadratic equation
	{
		const double delta = C*C-4.*D*B;
		if (delta >= 0)
		{
			n=2;
			x[0] = (std::sqrt(delta)-C)/(2.*B);
			x[1] = (-std::sqrt(delta)-C)/(2.*B);
		}
	}
	else if (C !=0)		// Linear equation
	{
		n=1;
		x[0] = -D/C;
	}
}


/* CACM Algorithm 326
   Roots of low order polynomials
   Author: Terence R.F.Nonweiler
   CACM  (Apr 1968) p269
   Translated into c and programmed by M.Dow
   ANUSF, Australian National University, Canberra, Australia
   m.dow@anu.edu.au
*/

int QUADROOTS(double p[5], double r[3][5])
{
/*
Array r[3][5]  p[5]
Roots of poly p[0] x^2 + p[1] x+p[2]=0
x=r[1][k] + i r[2][k]  k=1,2
*/
  double b,c,d;
  b=-p[1]/p[0]/2;
  c=p[2]/p[0];
  d=b*b-c;
  if(d>0)
  {
    if(b>0) b=r[1][2]=std::sqrt(d)+b;
    else    b=r[1][2]=-std::sqrt(d)+b;
    r[1][1]=c/b; r[2][1]=r[2][2]=0;
  }
  else
  {
    d=r[2][1]=std::sqrt(-d); r[2][2]=-d;
    r[1][1]=r[1][2]=b;
  }
  return(0);
}

int CUBICROOTS(double p[], double r[3][5])  
{
/*
Array r[3][5]  p[5]
Roots of poly p[0] x^3 + p[1] x^2...+p[3]=0
x=r[1][k] + i r[2][k]  k=1,...,3
Assumes 0<arctan(x)<pi/2 for x>0
*/

  double s,t,b,c,d;
  int k;
  if(p[0]!=1)
   for(k=1;k<4;k++) p[k]=p[k]/p[0]; p[0]=1;
   s=p[1]/3.0; t=s*p[1];
   b=0.5*(s*(t/1.5-p[2])+p[3]); t=(t-p[2])/3.0;
   c=t*t*t; d=b*b-c;
   if(d>=0)
   {
     d=std::pow((std::sqrt(d)+std::fabs(b)),1.0/3.0);
     if(d!=0)
     {
       if(b>0) b=-d;
       else b=d;
       c=t/b;
     }
     d=r[2][2]=std::sqrt(0.75)*(b-c); b=b+c;
     c=r[1][2]=-0.5*b-s;
     if((b>0 && s<=0) || (b<0 && s>0))
     {
       r[1][1]=c; r[2][1]=-d; r[1][3]=b-s;
       r[2][3]=0;
     }
     else
     {
       r[1][1]=b-s; r[2][1]=0; r[1][3]=c;
       r[2][3]=-d;
     }
   }  /* end 2 equal or complex roots */
   else
   {
     if(b==0)
       d=std::atan(1.0)/1.5;
     else
       d=std::atan(std::sqrt(-d)/std::fabs(b))/3.0;
     if(b<0)
       b=std::sqrt(t)*2.0;
     else
       b=-2.0*std::sqrt(t);
     c=std::cos(d)*b; t=-std::sqrt(0.75)*std::sin(d)*b-0.5*c;
     d=-t-c-s; c=c-s; t=t-s;
     if(std::fabs(c)>std::fabs(t))
       r[1][3]=c;
     else
       {
       r[1][3]=t; t=c;
       }
     if(std::fabs(d)>std::fabs(t))
       r[1][2]=d;
     else
     {
       r[1][2]=t; t=d;
     }
     r[1][1]=t;
     for(k=1;k<4;k++) r[2][k]=0;
   }
     return(0);
}

int BIQUADROOTS(double p[5], double r[3][5])
/* add _ if calling from fortran */
/*
Array r[3][5]  p[5]
Roots of poly p[0] x^4 + p[1] x^3...+p[4]=0
x=r[1][k] + i r[2][k]  k=1,...,4
*/
{
  double a,b,c,d,e;
  int k,j;
  if(p[0] != 1.0)
  {
    for(k=1;k<5;k++) p[k]=p[k]/p[0];
    p[0]=1;
  }
  e=0.25*p[1];
  b=2*e;
  c=b*b;
  d=0.75*c;
  b=p[3]+b*(c-p[2]);
  a=p[2]-d;
  c=p[4]+e*(e*a-p[3]);
  a=a-d;
  p[1]=0.5*a;
  p[2]=(p[1]*p[1]-c)*0.25;
  p[3]=b*b/(-64.0);
  if(p[3]<-1e-6)
  {
    CUBICROOTS(p,r);
    for(k=1;k<4;k++)
    {
      if(r[2][k]==0 && r[1][k]>0)
      {
	d=r[1][k]*4; a=a+d;
	if(a>=0 && b>=0)
	  p[1]=std::sqrt(d);
	else if(a<=0 && b<=0)
	  p[1]=std::sqrt(d);
	else p[1]=-std::sqrt(d);
	b=0.5*(a+b/p[1]);
	goto QUAD;
      }
    }
  }
  if(p[2]<0)
  {
    b=std::sqrt(c); d=b+b-a;
    p[1]=0; if(d>0) p[1]=std::sqrt(d);
  }
  else
  {
    if(p[1]>0)
      b=std::sqrt(p[2])*2.0+p[1];
    else
      b=-std::sqrt(p[2])*2.0+p[1];
    if(b!=0)
      p[1]=0;
    else
      {
	for(k=1;k<5;k++)
	{
	  r[1][k]=-e;
	  r[2][k]=0;
	}
	goto END;
      }
    }
QUAD:p[2]=c/b; QUADROOTS(p,r);
  for(k=1;k<3;k++)
  for(j=1;j<3;j++) r[j][k+2]=r[j][k];
  p[1]=-p[1]; p[2]=b; QUADROOTS(p,r);
  for(k=1;k<5;k++) r[1][k]=r[1][k]-e;
  END:;
  return(0);
}

std::string ScientificNotationWithFixedExponentDigits(const double number, const int expSize)
{
	std::ostringstream oss;
	std::string output;
	oss << std::scientific << std::setprecision(8) << number;
	std::size_t ePos = oss.str().find("e");
	std::size_t dPos = oss.str().find(".");
	if(ePos == 0)
	{
		//no exponent
		return "";
	}
	else if(dPos == 0)
	{
		//not decimal
		return "";
	}
	else
	{
		output = oss.str().substr(0, ePos) + oss.str().substr(ePos, 2);
		if(oss.str().size()-expSize > ePos+1)
			output += oss.str().substr(oss.str().size()-expSize, oss.str().size());
		else
		{
		 //expSize too big (or bug -> e used but no exponent?)
		}
		return output;
	}
}

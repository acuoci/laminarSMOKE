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

#include "kernel/thermo/AtomicComposition.h"
#include "kernel/thermo/AtomicElementMap.h"

namespace OpenSMOKE
{
	AtomicComposition::AtomicComposition()
	{
	}

	AtomicComposition::AtomicComposition(const AtomicComposition& orig) :
	names_(orig.names_), coefficients_(orig.coefficients_)
	{
	}

	AtomicComposition::AtomicComposition(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		// TODO
		//	names_.Load(fInput, fileFormat);
		//	coefficients_.Load(fInput, fileFormat);
		CheckForFatalError(false);
	}
		
	AtomicComposition::AtomicComposition( std::vector<std::string> const& names, std::vector<double> const& coefficients ) :
	names_(names), coefficients_(coefficients)
	{
	}    

	void AtomicComposition::operator() ( std::vector<std::string>& names, std::vector<double>& coefficients )
	{
		names_ = names;
		coefficients_ = coefficients;
	}

	AtomicComposition& AtomicComposition::operator=(const AtomicComposition& rhs)
	{
		names_ = rhs.names_;
		coefficients_ = rhs.coefficients_;
		return *this;
	}

	double AtomicComposition::mw() const
	{
		double mw_ = 0.;
                
		for(unsigned int i=0;i<coefficients_.size();i++)
				mw_ += AtomicWeights[names_[i]]*coefficients_[i];
		return mw_;
	}

	double AtomicComposition::GetCoefficient(const std::string name)
	{
		double coefficient = AtomicWeights[name];
		if (coefficient<=0)
			ErrorMessage(typeid(AtomicComposition).name(), "Not available atomic element: " + name);

		return coefficient;
	}

	void AtomicComposition::Save(std::ofstream& fOutput, const OpenSMOKE_File_Format fileFormat)
	{
		if (fileFormat == OPENSMOKE_BINARY_FILE)
		{
			unsigned int n = (unsigned int)(names_.size());
		
			fOutput.write(reinterpret_cast<char *>(&n),sizeof(n));
            for(unsigned int i=0;i<n;i++)
                fOutput.write(reinterpret_cast<char*>(&names_[i]),sizeof(std::string));
		}
		else
		{
			// TODO
			CheckForFatalError(false);
		}
	}

	void AtomicComposition::Load(std::ifstream& fInput, const OpenSMOKE_File_Format fileFormat)
	{
		int n;

		if(!fInput.read(reinterpret_cast<char*> (&n), sizeof(int)))
		{
			// TODO
			std::cout << "I was unable to read from binary file" << std::endl;
			CheckForFatalError(false);
		}

		// Reading vector elements
		names_.resize(n);
		for(int i = 0; i < n; i++)
                fInput.read(reinterpret_cast<char *>(&(names_[i])), sizeof(std::string));
	}

	void AtomicComposition::Status(std::ostream& fOut)
	{
		fOut << coefficients_.size() << std::endl;
		for(unsigned int i=0;i<names_.size();i++)
			fOut << names_[i] << "\t" << coefficients_[i] << std::endl;
	}

	bool AtomicComposition::CheckAtomicComposition() const
	{	
		for(unsigned int i=0;i<names_.size();i++)
		{
			AtomicElementMap::iterator it=AtomicWeights.find(names_[i]);
			if( it == AtomicWeights.end())
			{
				std::cout << "The following element is not available: " << names_[i] << std::endl;
				return false;
			}
		}
		return true;
	}

	const std::vector<std::string>&  AtomicComposition::element_names() const
	{
		return names_;
	}

	const std::vector<double>& AtomicComposition::element_coefficients() const 
	{
		return coefficients_;
	}
}
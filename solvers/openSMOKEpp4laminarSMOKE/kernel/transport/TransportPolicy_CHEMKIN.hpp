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

namespace OpenSMOKE
{

	TransportPolicy_CHEMKIN::TransportPolicy_CHEMKIN() {
	}

	TransportPolicy_CHEMKIN::TransportPolicy_CHEMKIN(const TransportPolicy_CHEMKIN& orig) 
	{
		name_transport_ = orig.name_transport_;
		transport_double_parameters_ = orig.transport_double_parameters_;
		transport_int_parameters_ = orig.transport_int_parameters_;
	}

	TransportPolicy_CHEMKIN::~TransportPolicy_CHEMKIN() {
	}

	void TransportPolicy_CHEMKIN::CopyTransportProperties(TransportPolicy_CHEMKIN& rhs) const 
	{
		rhs.name_transport_ = name_transport_;
		rhs.transport_double_parameters_ = transport_double_parameters_;
		rhs.transport_int_parameters_ = transport_int_parameters_;
	}

	bool TransportPolicy_CHEMKIN::ReadTransportFromStrings(const std::vector<std::string> &lines)
	{
		std::string line = lines[0];
		boost::trim(line);

		transport_int_parameters_.resize(1);
		transport_double_parameters_.resize(5);

		typedef boost::tokenizer<boost::char_separator<char> >  tokenizer_only_blanks;
		boost::char_separator<char> sep(" ");
		tokenizer_only_blanks tokens(line, sep);
        
		const std::size_t count = std::distance(tokens.begin(), tokens.end());
        
		try
		{
			if (count != 7)
				throw count;
		}
		catch(unsigned int k)
		{
			std::cout << "Wrong number of coefficients! Expected 7, Found " << k << std::endl;
			return false;
		}
        
		std::vector<std::string> words(7);
        
		int i=0;
		for (tokenizer_only_blanks::iterator tok_iter = tokens.begin();tok_iter != tokens.end(); ++tok_iter)
		{
			words[i] = *tok_iter;
			i++;
		}
             
		
		name_transport_ = words[0];
		boost::to_upper(name_transport_);
        
		try
		{
			boost::algorithm::trim(words[1]);
			boost::algorithm::trim(words[2]);
			boost::algorithm::trim(words[3]);
			boost::algorithm::trim(words[4]);
			boost::algorithm::trim(words[5]);
			boost::algorithm::trim(words[6]);

			transport_int_parameters_[0]=int(boost::lexical_cast<double>(words[1]));		// shape factor
			transport_double_parameters_[0]=boost::lexical_cast<double>(words[2]);		    // epsylon_over_kb_
			transport_double_parameters_[1]=boost::lexical_cast<double>(words[3]);			// sigma_
			transport_double_parameters_[2]=boost::lexical_cast<double>(words[4]);			// mu_
			transport_double_parameters_[3]=boost::lexical_cast<double>(words[5]);			// alfa_
			transport_double_parameters_[4]=boost::lexical_cast<double>(words[6]);			// zRot298_
		}
		catch(boost::bad_lexical_cast &)
		{
			std::cout << "Numerical conversion failure. Probably the transport properties coefficients are not written properly (they are not numbers)" << std::endl;
			return false;
		} 
        
		return true;
	}

	int TransportPolicy_CHEMKIN::shape_factor()	const 
	{
		return transport_int_parameters_[0];
	}
	
	double TransportPolicy_CHEMKIN::epsylon_over_kb() const 
	{
		return transport_double_parameters_[0];
	}
	
	double TransportPolicy_CHEMKIN::sigma()	const 
	{
		return transport_double_parameters_[1];
	}
	
	double TransportPolicy_CHEMKIN::mu() const 
	{
		return transport_double_parameters_[2];
	}    
	
	double TransportPolicy_CHEMKIN::alfa() const 
	{
		return transport_double_parameters_[3];
	}
	
	double TransportPolicy_CHEMKIN::zRot298() const 
	{
		return transport_double_parameters_[4];
	}

	void TransportPolicy_CHEMKIN::TransportStatus(std::ostream& fOut) const
	{
		fOut << "Shape factor:  " << transport_int_parameters_[0] << std::endl;
		fOut << "Eps/kb:        " << transport_double_parameters_[0] << std::endl;
		fOut << "Sigma:         " << transport_double_parameters_[1] << std::endl;
		fOut << "Mu:            " << transport_double_parameters_[2] << std::endl;
		fOut << "Alfa:          " << transport_double_parameters_[3] << std::endl;
		fOut << "zRot298:       " << transport_double_parameters_[4] << std::endl;   
	}
}
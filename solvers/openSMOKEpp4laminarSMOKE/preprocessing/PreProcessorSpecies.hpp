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

namespace OpenSMOKE
{
	template<typename Species>
	PreProcessorSpecies<Species>::PreProcessorSpecies() {
	}

	template<typename Species>
	PreProcessorSpecies<Species>::PreProcessorSpecies(const PreProcessorSpecies& orig) {
	}

	template<typename Species>
	PreProcessorSpecies<Species>::~PreProcessorSpecies() {
	}

	template<typename Species>
	template<typename Thermo, typename Transport, typename Kinetics >
	PreProcessorSpecies<Species>::PreProcessorSpecies(const Thermo&  thermo, const Transport& transport, const Kinetics& kinetics, std::ostream& flog) :
		flog_(flog)
	{

		this->species_.resize(kinetics.names_species().size());
		this->names_.push_back("species");

		if (kinetics.names_species().size() > 20)
			std::cout << " * Extracting " << kinetics.names_species().size() << " species from databases: ";

                typename Transport::map_species transport_species=transport.species();
		int count=0;
		for(unsigned int i=0;i<kinetics.names_species().size();i++)
		{
			if (kinetics.names_species().size() > 20)
			{
				if (i%(kinetics.names_species().size()/20) == 0) 
					std::cout << "%";
				if (i==kinetics.names_species().size()-1)
					std::cout << std::endl;
			}
			
			bool iThermoFound = true;
			typename Thermo::map_species thermo_species=thermo.species();
			typename Thermo::map_species::iterator thermo_it=thermo_species.find(kinetics.names_species()[i]);
			if( thermo_it == thermo_species.end())
			{
					iThermoFound = false;
					std::cout << "Error: This species is not available in thermodynamic database: " << kinetics.names_species()[i] << std::endl;
			}

			bool iTransportFound = true;
			
			typename Transport::map_species::iterator transport_it=transport_species.find(kinetics.names_species()[i]);
			if( transport_it == transport_species.end())
			{
				iTransportFound = false;
				std::cout << "Error: This species is not available in transport database: " << kinetics.names_species()[i] << std::endl;
			}
			
			if (iThermoFound == true && iTransportFound == true)
			{
				count++;
				this->species_[i](thermo_species[kinetics.names_species()[i]], transport_species[kinetics.names_species()[i]]);
				this->names_.push_back(kinetics.names_species()[i]);
			}
		}

		CheckForFatalError( (count == kinetics.names_species().size()) );
		atomicTable_.CreateAtomicElementLists(this->species_);
	}

	template<typename Species>
	template<typename Thermo, typename Transport >
	PreProcessorSpecies<Species>::PreProcessorSpecies(const Thermo&  thermo, const Transport& transport,std::vector<std::string> kinetics, std::ostream& flog) :
		flog_(flog)
	{

		typename Transport::map_species transport_species = transport.species();

		this->species_.resize(transport_species.size());
		this->names_.push_back("species");

		if( transport_species.size() > 20)
			std::cout << " * Extracting " << transport_species.size() << " species from databases: ";

		
		int count = 0;
		int i = 0;
		for (typename Transport::map_species::const_iterator iterator = transport_species.begin(); iterator != transport_species.end(); iterator++)
		{
			bool iThermoFound = true;
			typename Thermo::map_species thermo_species = thermo.species();
			typename Thermo::map_species::iterator thermo_it = thermo_species.find(iterator->second.name_transport());
			if (thermo_it == thermo_species.end())
			{
				iThermoFound = false;
				std::cout << "Error: This species is not available in thermodynamic database: " << iterator->second.name_transport() << std::endl;
			}

			if (iThermoFound == true)
			{
				count++;
				this->species_[i](thermo_species[iterator->second.name_transport()], transport_species[iterator->second.name_transport()]);
				this->names_.push_back(iterator->second.name_transport());
			}

			i++;
		}

		CheckForFatalError((count == transport_species.size()));
		atomicTable_.CreateAtomicElementLists(this->species_);
	}

	template<typename Species>
	template<typename Thermo, typename Kinetics >
	PreProcessorSpecies<Species>::PreProcessorSpecies(const Thermo&  thermo, const Kinetics& kinetics, std::ostream& flog) :
		flog_(flog)
	{
		this->species_.resize(kinetics.names_species().size());
                this->names_.resize(kinetics.names_species().size()+1);
		this->names_[0]="species";

		if (kinetics.names_species().size() > 20)
			std::cout << " * Extracting " << kinetics.names_species().size() << "/" << thermo.species().size() << " species from thermodynamic database: ";
        
                typename Thermo::map_species thermo_species=thermo.species();
		int count=0;
		for(unsigned int i=0;i<kinetics.names_species().size();i++)
		{
			if (kinetics.names_species().size() > 20)
			{
				if (i%(kinetics.names_species().size()/20) == 0) 
					std::cout << "%";
			}
                        
                        if (i==kinetics.names_species().size()-1)
				std::cout << std::endl;
			
			bool iThermoFound = true;
			
			typename Thermo::map_species::iterator thermo_it=thermo_species.find(kinetics.names_species()[i]);
			if( thermo_it == thermo_species.end())
			{
				if (kinetics.names_species()[i] == "R" || kinetics.names_species()[i] == "RH")
				{
					iThermoFound == true;
				}
				else
				{
					iThermoFound = false;
					std::cout << "Error: This species is not available in thermodynamic database: " << kinetics.names_species()[i] << std::endl;
				}
			}
			
			if (iThermoFound == true)
			{
				count++;
				this->species_[i](thermo_species[kinetics.names_species()[i]]);
                                this->names_[i+1] = kinetics.names_species()[i];
			}
		}

		CheckForFatalError( (count == kinetics.names_species().size()) );
		atomicTable_.CreateAtomicElementLists(this->species_);
	}

	template<typename Species>
	template<typename Thermo, typename Kinetics >
	PreProcessorSpecies<Species>::PreProcessorSpecies(const Thermo&  thermo, const std::vector<std::string> additional_species, const Kinetics& kinetics, std::ostream& flog) :
		flog_(flog)
	{
		std::vector<std::string> vector_of_species(additional_species.size() + kinetics.names_additional_species().size());
		for(unsigned int i=0;i<additional_species.size();++i)
			vector_of_species[i] = additional_species[i];
		
		for(unsigned int i=0;i<kinetics.names_additional_species().size();++i)
			vector_of_species[i+additional_species.size()] = kinetics.names_additional_species()[i];
		
		const std::size_t total_number_of_species = vector_of_species.size();
		this->species_.resize( total_number_of_species );
		this->names_.push_back("species");
		
		if (total_number_of_species > 20)
			std::cout << " * Extracting " << total_number_of_species << " species from thermodynamic database: ";
	
		int count=0;
		for(unsigned int i=0;i<total_number_of_species;i++)
		{
			if (total_number_of_species > 20)
			{
				if (i%(total_number_of_species/20) == 0) 
					std::cout << "%";
				if (i==total_number_of_species-1)
					std::cout << std::endl;
			}
			
			bool iThermoFound = true;
			typename Thermo::map_species thermo_species=thermo.species();
			typename Thermo::map_species::iterator thermo_it=thermo_species.find(vector_of_species[i]);
			if( thermo_it == thermo_species.end())
			{
					iThermoFound = false;
					std::cout << "Error: This species is not available in thermodynamic database: " << vector_of_species[i] << std::endl;
			}
			
			if (iThermoFound == true)
			{
				count++;
				this->species_[i](thermo_species[vector_of_species[i]]);
				this->names_.push_back(vector_of_species[i]);
			}
		}
		
		CheckForFatalError( (count == vector_of_species.size()) );
		atomicTable_.CreateAtomicElementLists(this->species_);
	}

	template<typename Species>
	template<typename Thermo >
	PreProcessorSpecies<Species>::PreProcessorSpecies(const Thermo&  thermo, std::ostream& flog) :
		flog_(flog)
	{
		typename Thermo::map_species thermo_species=thermo.species();

		this->species_.resize(thermo_species.size());
		this->names_.push_back("species");

		unsigned int i=0;
		for (typename Thermo::map_species::const_iterator it = thermo_species.begin(); it != thermo_species.end(); ++it)
		{
			std::string name = it->first;
			this->species_[i](thermo_species[name]);
			this->names_.push_back(name);
			i++;
		}

		atomicTable_.CreateAtomicElementLists(this->species_);
	}

	template<typename Species>
	template<typename Thermo, typename Transport >
	PreProcessorSpecies<Species>::PreProcessorSpecies(std::ostream& flog, const Thermo&  thermo, const Transport& transport) :
		flog_(flog)
	{
		typename Thermo::map_species thermo_species=thermo.species();

		this->species_.resize(thermo_species.size());
		this->names_.push_back("species");

		unsigned int i=0;
		for (typename Thermo::map_species::const_iterator it = thermo_species.begin(); it != thermo_species.end(); ++it)
		{
			std::string name = it->first;
			
			bool iTransportFound = true;
			typename Transport::map_species transport_species=transport.species();
			typename Transport::map_species::iterator transport_it=transport_species.find(name);
			if( transport_it == transport_species.end())
			{
				iTransportFound = false;
				std::cout << "Error: This species is not available in transport database: " << name << std::endl;
			}
			
			if (iTransportFound == true)
			{			
				this->species_[i](thermo_species[name], transport_species[name]);
				this->names_.push_back(name);
				i++;
			}
		}

		CheckForFatalError( (i == thermo_species.size()) );

		atomicTable_.CreateAtomicElementLists(this->species_);
	}

	template<typename Species>
	void PreProcessorSpecies<Species>:: WriteASCIIFileOldStyle(const std::string file_name) const
	{
		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		this->WriteThermodynamicDataOnASCIIFileOldStyle(fOutput);

		// Element matrix coefficients
		fOutput << atomicTable_.element_weights_list().size() << " ";
		fOutput << this->species().size() << std::endl;
		for (int i=0;i<atomicTable_.element_weights_list().size();i++)
		{
			for (unsigned int j=0;j<this->species().size();j++)
				fOutput << atomicTable_.element_coefficients_list()(j,i) << " ";
			fOutput << std::endl;
		}
	
		// Element molecular weights
		fOutput << atomicTable_.element_weights_list().size() << std::endl;
		for (int i=0;i<atomicTable_.element_weights_list().size();i++)
			fOutput << atomicTable_.element_weights_list()(i) << " ";
		fOutput << std::endl;

		// Element names
		for (unsigned int i=0;i<atomicTable_.element_names_list().size();i++)
			fOutput << atomicTable_.element_names_list()[i] << std::endl;

		// Transport
		this->WriteTransportDataOnASCIIFileOldStyle(fOutput);

		fOutput.close();
	}

	template<typename Species>
	void PreProcessorSpecies<Species>::WriteASCIIFile(const std::string file_name) const
	{
		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		// Thermodynamics
		this->WriteThermodynamicsDataOnASCIIFile(fOutput);

		// Transport properties
		this->WriteTransportDataOnASCIIFile(fOutput);

		fOutput.close();
	}

	template<typename Species>
	void PreProcessorSpecies<Species>::WriteXMLFile(std::stringstream& xml_string) const
	{
		xml_string << "<NumberOfElements>" << std::endl;
		xml_string << atomicTable_.element_coefficients_list().cols() << std::endl;
		xml_string << "</NumberOfElements>" << std::endl;

		xml_string << "<NamesOfElements>" << std::endl;
		for (unsigned int i = 0; i < atomicTable_.element_names_list().size(); i++)
			xml_string << atomicTable_.element_names_list()[i] << " ";
		xml_string << std::endl;
		xml_string << "</NamesOfElements>" << std::endl;

		xml_string << "<AtomicComposition>" << std::endl;
		for (unsigned int i = 0; i < this->species_.size(); i++)
		{
			for (int j = 0; j < atomicTable_.element_weights_list().size(); j++)
			{
				if (atomicTable_.element_coefficients_list()(i, j) == 0.)
					xml_string << "0. ";
				else
					xml_string << std::setprecision(10) << atomicTable_.element_coefficients_list()(i, j) << " ";
			}
			xml_string << std::endl;
		}
		xml_string << "</AtomicComposition>" << std::endl;

		// Names of species
		this->WriteSpeciesOnXMLFile(xml_string);

		// Thermodynamics
		this->WriteThermodynamicsDataOnXMLFile(xml_string);

		// Transport properties
		this->WriteTransportDataOnXMLFile(xml_string);
	}

	template<typename Species>
	bool PreProcessorSpecies<Species>::SpeciesBundling(std::stringstream& xml_string, const double epsilon)
	{
		return this->WriteSpeciesBundlingOnXMLFile(xml_string, epsilon);
	}
	
	template<typename Species>
	bool PreProcessorSpecies<Species>::ThermodynamicConsistency()
	{
		return this->CheckThermodynamicConsistency(flog_);
	}

	template<typename Species>
	bool PreProcessorSpecies<Species>::WriteElementTableOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput << "------------------------------"	<< std::endl;
		fOutput << "  Elements     Atomic Weight  "	<< std::endl;
		fOutput << "------------------------------"	<< std::endl;
		for (int i=0;i<atomicTable_.element_weights_list().size();i++)
		{
			fOutput << "  " << i+1 << ".\t";
 			fOutput << std::setw(9) << std::left << atomicTable_.element_names_list()[i];
			fOutput << std::fixed   << std::setw(12) << std::setprecision(8) << std::right << atomicTable_.element_weights_list()[i];
			fOutput << std::endl;
		}
		fOutput << "------------------------------"	<< std::endl;
		fOutput << std::endl << std::endl;

		fOutput << "------------------------------------------------------------------------------------";
		fOutput << std::setfill('-');
		fOutput << std::setw(8*(atomicTable_.element_weights_list().size()-1)+2) << "-" << std::endl;
		fOutput << std::setfill(' ');
		fOutput << "                                                Molecular        Temperature        Elements" << std::endl;
		fOutput << "  Species                  Phase  Charge         weight         Low      High";
		for (int i=0;i<atomicTable_.element_weights_list().size();i++)
			fOutput << std::setw(8) << std::right << atomicTable_.element_names_list()[i];
		fOutput << std::endl;
		fOutput << "------------------------------------------------------------------------------------";
		fOutput << std::setfill('-');
		fOutput << std::setw(8*(atomicTable_.element_weights_list().size()-1)+2) << "-" << std::endl;
		fOutput << std::setfill(' ');
		
		for (unsigned int i=0;i<this->species_.size();i++)
		{
			fOutput << std::right << std::setw(5) << i+1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << this->species_[i].name();
			fOutput << std::setw(3)  << std::right << this->species_[i].phase();
			fOutput << std::setw(7)  << std::right << 0;
			fOutput << std::setw(20) << std::fixed << std::right << std::setprecision(6)  << this->species_[i].MolecularWeight();
			fOutput << std::setw(10) << std::right << std::setprecision(2)  << this->species_[i].Tlow();
			fOutput << std::setw(10) << std::right << std::setprecision(2)  << this->species_[i].Thigh();
			for (int j=0;j<atomicTable_.element_weights_list().size();j++)
				fOutput << std::setw(8) << std::setprecision(0) << std::fixed << std::right << atomicTable_.element_coefficients_list()(i,j);
			fOutput << std::endl;
		}
		fOutput << "------------------------------------------------------------------------------------";
		fOutput << std::setfill('-');
		fOutput << std::setw(8*(atomicTable_.element_weights_list().size()-1)+2) << "-" << std::endl;
		fOutput << std::setfill(' ');
		fOutput << std::endl;

		this->WriteReorderedElementTableOnASCIIFile(fOutput);

		this->WriteTransportTableOnASCIIFile(fOutput);

		return true;
	}

	template<typename Species>
	bool PreProcessorSpecies<Species>::WriteReorderedElementTableOnASCIIFile(std::ostream& fOutput) const
	{
		unsigned int nrows = this->species_.size();
		unsigned int ncols = atomicTable_.element_weights_list().size();

		std::vector<size_t> indices_reordered(nrows);
		for (unsigned int i = 0; i < this->species_.size(); i++)
			indices_reordered[i] = i;

		Eigen::MatrixXd elements_reordered(nrows, ncols);
		for (unsigned int i = 0; i < this->species_.size(); i++)
			for (int j = 0; j < atomicTable_.element_weights_list().size(); j++)
				elements_reordered(i, j) = atomicTable_.element_coefficients_list()(i, j);

		for(int jcol=ncols-1;jcol>=0;jcol--)
		{
			std::vector<size_t> indices_current = indices_reordered;
			Eigen::MatrixXd elements_current = elements_reordered;

			std::vector<double> element(nrows);
			for (unsigned int i = 0; i < this->species_.size(); i++)
				element[i] = elements_current(i, jcol);

			std::vector<size_t> indices = SortAndTrackIndicesIncreasing(element);
			for (unsigned int i = 0; i < this->species_.size(); i++)
			{
				elements_reordered.row(i) = elements_current.row(indices[i]);
				indices_reordered[i] = indices_current[indices[i]];
			}
		}

		fOutput << "------------------------------------------------------------------------------------";
		fOutput << std::setfill('-');
		fOutput << std::setw(8 * (atomicTable_.element_weights_list().size() - 1) + 2) << "-" << std::endl;
		fOutput << std::setfill(' ');
		fOutput << "                                                Molecular        Temperature        Elements" << std::endl;
		fOutput << "  Species                  Phase  Charge         weight         Low      High";
		for (int i = 0; i<atomicTable_.element_weights_list().size(); i++)
			fOutput << std::setw(8) << std::right << atomicTable_.element_names_list()[i];
		fOutput << std::endl;
		fOutput << "------------------------------------------------------------------------------------";
		fOutput << std::setfill('-');
		fOutput << std::setw(8 * (atomicTable_.element_weights_list().size() - 1) + 2) << "-" << std::endl;
		fOutput << std::setfill(' ');

		for (unsigned int ii = 0; ii<this->species_.size(); ii++)
		{
			unsigned int i = indices_reordered[ii];

			fOutput << std::right << std::setw(5) << i + 1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left << this->species_[i].name();
			fOutput << std::setw(3) << std::right << this->species_[i].phase();
			fOutput << std::setw(7) << std::right << 0;
			fOutput << std::setw(20) << std::fixed << std::right << std::setprecision(6) << this->species_[i].MolecularWeight();
			fOutput << std::setw(10) << std::right << std::setprecision(2) << this->species_[i].Tlow();
			fOutput << std::setw(10) << std::right << std::setprecision(2) << this->species_[i].Thigh();
			for (int j = 0; j<atomicTable_.element_weights_list().size(); j++)
				fOutput << std::setw(8) << std::setprecision(0) << std::fixed << std::right << atomicTable_.element_coefficients_list()(i, j);
			fOutput << std::endl;
		}
		fOutput << "------------------------------------------------------------------------------------";
		fOutput << std::setfill('-');
		fOutput << std::setw(8 * (atomicTable_.element_weights_list().size() - 1) + 2) << "-" << std::endl;
		fOutput << std::setfill(' ');
		fOutput << std::endl;

		this->WriteTransportTableOnASCIIFile(fOutput);

		return true;
	}
}



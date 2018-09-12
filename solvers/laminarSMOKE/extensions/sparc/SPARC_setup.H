/*-----------------------------------------------------------------------*\
|                                                                         |
|                    ╔═══╦═╗╔═╦═══╦╗╔═╦═══╗                               |
|                    ║╔═╗║║╚╝║║╔═╗║║║╔╣╔══╝                               | 
|   ╔╗╔══╦╗╔╦╦═╗╔══╦═╣╚══╣╔╗╔╗║║ ║║╚╝╝║╚══╗                               |
|   ║║║╔╗║╚╝╠╣╔╗╣╔╗║╔╩══╗║║║║║║║ ║║╔╗║║╔══╝                               |
|   ║╚╣╔╗║║║║║║║║╔╗║║║╚═╝║║║║║║╚═╝║║║╚╣╚══╗                               |
|   ╚═╩╝╚╩╩╩╩╩╝╚╩╝╚╩╝╚═══╩╝╚╝╚╩═══╩╝╚═╩═══╝                               |
|                                                                         |
|                                                                         |
|   Authors: A. Cuoci, G. D'Alessio, A. Parente                           |
|                                                                         |
|   Contacts: Alberto Cuoci                                               |
|   email: alberto.cuoci@polimi.it                                        |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano (Italy)                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of laminarSMOKE solver.                             |
|                                                                         |
|   License                                                               |
|                                                                         |
|   Copyright(C) 2018 A. Cuoci                                            |
|   laminarSMOKE is free software: you can redistribute it and/or modify  |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   laminarSMOKE is distributed in the hope that it will be useful,       |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with laminarSMOKE. If not, see <http://www.gnu.org/licenses/>.  |
|                                                                         |
\*-----------------------------------------------------------------------*/

#if SPARC==1

// Clustering
int nclusters = 0;
std::vector<Eigen::VectorXi> classifier_indices;
std::vector<OpenSMOKE::ThermodynamicsMap_CHEMKIN*> vector_thermodynamicsMapXML;
std::vector<OpenSMOKE::KineticsMap_CHEMKIN*> vector_kineticsMapXML;
std::vector<OpenSMOKE::TransportPropertiesMap_CHEMKIN*> vector_transportMapXML; 

// VQ2
const dictionary& vq2Dictionary = solverOptions.subDict("VQ2");
const Switch iVQ2 = Switch(vq2Dictionary.lookup(word("vq2")));
const Switch iVQ2dubug = Switch(vq2Dictionary.lookup(word("debug")));
autoPtr<std::ofstream> fVQ2;
label vq2IntegralUpdate  = 10;
label vq2IntegralCounter = vq2IntegralUpdate;
ClassifierVQ2 VQ2;

if (iVQ2 == true)
{
	// This folder contains a 'kinetics' and a 'data' subfolders
	Foam::string main_folder_name = vq2Dictionary.lookup("folder");	
	boost::filesystem::path path_main_folder = main_folder_name;
	boost::filesystem::path path_kinetics_folder = path_main_folder / "kinetics";
	boost::filesystem::path path_data_folder = path_main_folder / "data";	

	int scaling = 0;
	word scaling_type(vq2Dictionary.lookup("scaling"));
	if (scaling_type == "AUTO")		scaling = 0;
	else if (scaling_type == "VAST")	scaling = 1;
	else if (scaling_type == "PARETO")	scaling = 2;
	else
	{
		Info << "Wrong scaling option: AUTO || VAST || PARETO" << endl;
		abort();
	}

	label eigens  = vq2Dictionary.lookupOrDefault("eigens", 1);

	std::cout << "Opening files containing VQ2 data..." << std::endl;
	VQ2.Setup(path_data_folder.c_str(), eigens, scaling);
	nclusters = VQ2.nc();
	
	if (VQ2.nv() != (thermodynamicsMapXML->NumberOfSpecies()+1))
		OpenSMOKE::ErrorMessage("Importing data for VQ2", "Mismatch between the data and the kinetic mechanism");

	// Read thermodynamics and kinetics maps
	vector_thermodynamicsMapXML.resize(nclusters);
	vector_transportMapXML.resize(nclusters);
	vector_kineticsMapXML.resize(nclusters);
	classifier_indices.resize(nclusters);
	
	for (int j = 0; j < nclusters; j++)
	{
		std::stringstream label; label << j;
		const std::string folder_name = "kinetics." + label.str();

		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

		// Mapping
		{
			unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
			classifier_indices[j].resize(NClocal);
			classifier_indices[j].setConstant(-1);
			for (unsigned int i=0;i<NClocal;i++)
				for (unsigned int k=0;k<NC;k++)
					if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
						classifier_indices[j](i) = k;

			for (unsigned int i=0;i<NClocal;i++)
				if (classifier_indices[j](i) == -1)
				{
					Info << "Kinetic mechanism #" << i << endl;
					Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
					abort();
				}
		}
	}

	if ( !boost::filesystem::exists("vq2.out") )
        {
                fVQ2.reset(new std::ofstream("vq2.out", std::ios::out));
                fVQ2().setf(std::ios::scientific);

                fVQ2() << std::setw(20) << "time";
                fVQ2() << std::setw(20) << "mean_species";
                fVQ2() << std::setw(20) << "max_species";
                fVQ2() << std::setw(20) << "mean_reactions";
                fVQ2() << std::setw(20) << "max_reactions";
                fVQ2() << std::setw(20) << "meanCPUchem";
                fVQ2() << std::setw(20) << "maxCPUchem";
                fVQ2() << std::endl;
        }
        else
        {
                fVQ2.reset(new std::ofstream("vq2.out", std::ios::app));
                fVQ2().setf(std::ios::scientific);
        }
}

// SOM
int npca = 0;
Eigen::MatrixXd pca_weights;
Eigen::VectorXd pca_mu;
Eigen::VectorXd pca_sigma;
int nf;
Eigen::VectorXi listretspecies;
const dictionary& somDictionary = solverOptions.subDict("SelfOrganizingMap");
const Switch iSOM = Switch(somDictionary.lookup(word("som")));
autoPtr<std::ofstream> fSom;
label somIntegralUpdate  = 10;
label somIntegralCounter = somIntegralUpdate;
if (iSOM == true)
{
	Foam::string file_name = somDictionary.lookup("xml");	
	boost::filesystem::path path_xml_input_file = file_name;
	Foam::string folder_name = somDictionary.lookup("folder");	
	boost::filesystem::path path_kinetics_folder = folder_name;
	
	std::cout << "Opening XML file containing SOM data..." << std::endl;

	// Open the XML file
	rapidxml::xml_document<> doc;
	std::vector<char> xml_string;
	OpenSMOKE::OpenInputFileXML(doc, xml_string, path_xml_input_file);
	rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");

	// Read number of clusters
	{
		rapidxml::xml_node<>* nclusters_node = opensmoke_node->first_node("classes");
		try
		{
			nclusters = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nclusters_node->value())));
			std::cout << " * Number of clusters: " << nclusters << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of clusters.");
		}
	}

	// Read number of principal components
	{
		rapidxml::xml_node<>* npca_node = opensmoke_node->first_node("principal-components");
		try
		{
			npca = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(npca_node->value())));
			std::cout << " * Number of principal components: " << npca << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number principal components.");
		}
	}

	// Read number retained species 
	{
		rapidxml::xml_node<>* nretspecies_node = opensmoke_node->first_node("number-retained-species");
		try
		{
			nf = boost::lexical_cast<unsigned int>(boost::trim_copy(std::string(nretspecies_node->value())));
			std::cout << " * Number of retained species: " << nf << std::endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the number of retained species.");
		}
	}

	// Read list retained species 
	if (nf != thermodynamicsMapXML->NumberOfSpecies())
	{
		rapidxml::xml_node<>* listretspecies_node = opensmoke_node->first_node("list-retained-species");
		try
		{
			listretspecies.resize(nf);

			std::stringstream fInput;
			fInput << listretspecies_node->value();

			for (unsigned int j = 0; j < nf; j++)
				fInput >> listretspecies(j);	// (1-index based)

			Info << "Retained species: " << nf << endl;
			for (unsigned int j = 0; j < nf; j++)
				Info << thermodynamicsMapXML->NamesOfSpecies()[listretspecies(j) - 1] << endl;
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading the list of retained species.");
		}
	}

	pca_weights.resize(nf+1,npca);
	pca_mu.resize(nf+1);
	pca_sigma.resize(nf+1);

	// Read mean values of filtered variables
	{
		rapidxml::xml_node<>* mu_node = opensmoke_node->first_node("mu");
		try
		{
			std::stringstream fInput;
			fInput << mu_node->value();

			for (unsigned int j = 0; j < nf+1; j++)
				fInput >> pca_mu(j);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading means of filtered variables.");
		}
	}

	// Read std deviations of filtered variables
	{
		rapidxml::xml_node<>* sigma_node = opensmoke_node->first_node("sigma");
		try
		{
			std::stringstream fInput;
			fInput << sigma_node->value();

			for (unsigned int j = 0; j < nf + 1; j++)
				fInput >> pca_sigma(j);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading std deviations of filtered variables.");
		}
	}

	// Read PCA weights
	{
		rapidxml::xml_node<>* w_node = opensmoke_node->first_node("weights");
		try
		{
			std::stringstream fInput;
			fInput << w_node->value();

			for (unsigned int j = 0; j < nf + 1; j++)
				for (unsigned int i = 0; i < npca; i++)
					fInput >> pca_weights(j,i);
		}
		catch (...)
		{
			OpenSMOKE::ErrorMessage("Importing data from MATLAB(R) preprocessing", "Error in reading weights of PCA.");
		}
	}

	// Read thermodynamics and kinetics maps
	vector_thermodynamicsMapXML.resize(nclusters);
	vector_transportMapXML.resize(nclusters);
	vector_kineticsMapXML.resize(nclusters);
	classifier_indices.resize(nclusters);
	for (int j = 0; j < nclusters; j++)
	{
		std::stringstream label; label << j;
		const std::string folder_name = "kinetics." + label.str();

		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string, path_kinetics_folder / folder_name / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		vector_thermodynamicsMapXML[j] = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc);
		vector_transportMapXML[j] = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc);
		vector_kineticsMapXML[j] = new OpenSMOKE::KineticsMap_CHEMKIN(*vector_thermodynamicsMapXML[j], doc);
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << "Time to read XML file: " << tEnd - tStart << std::endl;

		// Mapping
		{
			unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			unsigned int NClocal = vector_thermodynamicsMapXML[j]->NumberOfSpecies();
			classifier_indices[j].resize(NClocal);
			classifier_indices[j].setConstant(-1);
			for (unsigned int i=0;i<NClocal;i++)
				for (unsigned int k=0;k<NC;k++)
					if (thermodynamicsMapXML->NamesOfSpecies()[k] == vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i])
						classifier_indices[j](i) = k;

			for (unsigned int i=0;i<NClocal;i++)
				if (classifier_indices[j](i) == -1)
				{
					Info << "Kinetic mechanism #" << i << endl;
					Info << "Species " << vector_thermodynamicsMapXML[j]->NamesOfSpecies()[i] << " is not included in the main kinetic mechanism" << endl;
					abort();
				}
		}
	}

	if ( !boost::filesystem::exists("som.out") )
                        {
                                fSom.reset(new std::ofstream("som.out", std::ios::out));
                                fSom().setf(std::ios::scientific);

                                fSom() << std::setw(20) << "time";
                                fSom() << std::setw(20) << "mean_species";
                                fSom() << std::setw(20) << "max_species";
                                fSom() << std::setw(20) << "mean_reactions";
                                fSom() << std::setw(20) << "max_reactions";
                                fSom() << std::setw(20) << "meanCPUchem";
                                fSom() << std::setw(20) << "maxCPUchem";
                                fSom() << std::endl;
                        }
                        else
                        {
                                fSom.reset(new std::ofstream("som.out", std::ios::app));
                                fSom().setf(std::ios::scientific);
                        }	
}

#endif

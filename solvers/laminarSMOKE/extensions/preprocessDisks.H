Switch diskPreprocess = true;
PtrList<volScalarField> sourceFromDisk;

if (diskSourceTerms == true)
{
	List<int> list_disks;
	std::vector< Eigen::VectorXi > disk_topology_indices;
	std::vector< Eigen::VectorXd > disk_source_terms;

	const dictionary& diskSourceTermsDictionary = solverOptions.subDict("DiskSourceTerms");
	
	Foam::string folder_disks = diskSourceTermsDictionary.lookup("folder");
	boost::filesystem::path path_folder_disks = folder_disks;

	diskPreprocess = Switch(diskSourceTermsDictionary.lookupOrDefault(word("preProcess"), word("on")));

	list_disks = readList<int>(diskSourceTermsDictionary.lookup("disks"));
	Switch exclude_negative_disk_source_terms = Switch(diskSourceTermsDictionary.lookupOrDefault(word("excludeNegativeSourceTerms"), word("off")));

	Foam::string time_target = diskSourceTermsDictionary.lookup("targetTime");

	scalar correction_coefficient = readScalar(diskSourceTermsDictionary.lookup("correctionCoefficient"));

	if (diskPreprocess == true)
	{
		disk_topology_indices.resize(list_disks.size());
		disk_source_terms.resize(list_disks.size());
		std::vector<double> disk_topology_volume(list_disks.size());

		// Read topology
		{
			rapidxml::xml_document<> doc;
			std::vector<char> xml_string;
			OpenSMOKE::OpenInputFileXML(doc, xml_string, path_folder_disks / "DiskTopology.xml");

			rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
			rapidxml::xml_node<>* number_of_disks_node = opensmoke_node->first_node("Disks");

			try
			{
				const int ndisks = boost::lexical_cast<int>(boost::trim_copy(std::string(number_of_disks_node->value())));
			}
			catch(...)
			{
				Info << "Wrong number of disks in DiskTopology.xml file" << endl;
				abort();
			}

			try
			{
				for(unsigned int i=0;i<list_disks.size();i++)
				{
					Info << "Reading topology for Disk " << list_disks[i] << endl;

					std::string disk_name = "Disk." + std::to_string(list_disks[i]);
					rapidxml::xml_node<>* disk_node = opensmoke_node->first_node(disk_name.c_str());
					rapidxml::xml_node<>* number_points_node = disk_node->first_node("NumberOfPoints");
					rapidxml::xml_node<>* volume_node = disk_node->first_node("Volume");
					rapidxml::xml_node<>* cells_node = disk_node->first_node("Cells");

					const int np = boost::lexical_cast<int>(boost::trim_copy(std::string(number_points_node->value())));
					disk_topology_indices[i].resize(np);


					disk_topology_volume[i] = boost::lexical_cast<double>(boost::trim_copy(std::string(volume_node->value())));

					std::cout << " * Cells: " << np << " Volume: " << disk_topology_volume[i] << std::endl;

					std::stringstream data;
					data << cells_node->value();
					for(unsigned int j=0;j<np;j++)
					{
						std::string dummy;
						data >> dummy; disk_topology_indices[i](j) = std::stoi(dummy);
						data >> dummy; 
					}
				}
			}
			catch(...)
			{
				Info << "Wrong input data in DiskTopology.xml file" << endl;
				abort();
			}
		}

		// Read disks
		for(unsigned int i=0;i<list_disks.size();i++)
		{
			Info << "Reading source terms for Disk " << list_disks[i] << endl;

			std::string disk_name = "Disk." + std::to_string(list_disks[i]) + ".source." + time_target + ".xml";

			rapidxml::xml_document<> doc;
			std::vector<char> xml_string;
			OpenSMOKE::OpenInputFileXML(doc, xml_string, path_folder_disks / disk_name);

			rapidxml::xml_node<>* opensmoke_node = doc.first_node("opensmoke");
			rapidxml::xml_node<>* number_of_species_node = opensmoke_node->first_node("number-species");
			rapidxml::xml_node<>* slice_volume_node = opensmoke_node->first_node("slice-volume");
			rapidxml::xml_node<>* source_terms_node = opensmoke_node->first_node("source-terms");

			// Source terms
			disk_source_terms[i].resize(thermodynamicsMapXML->NumberOfSpecies());
			disk_source_terms[i].setZero();

			try
			{
				const int number_of_species = boost::lexical_cast<int>(boost::trim_copy(std::string(number_of_species_node->value())));
				const double slice_volume = boost::lexical_cast<double>(boost::trim_copy(std::string(slice_volume_node->value())));

				std::stringstream data;
				data << source_terms_node->value();
				for(unsigned int j=0;j<number_of_species;j++)
				{
					std::string name_species;
					data >> name_species;
					const int j_species = thermodynamicsMapXML->IndexOfSpecies(name_species)-1;
					if (j_species < 0)
						Info << "Warning: " << name_species << " is not available in the current kinetic mechanism" << endl;

					std::string dummy;				
					data >> dummy; 
					data >> dummy; 				
					data >> dummy; 

					disk_source_terms[i](j_species) = correction_coefficient*std::stod(dummy)*slice_volume/disk_topology_volume[i];	// [kg/s/m3]
					if (exclude_negative_disk_source_terms == true)
						if (disk_source_terms[i](j_species) < 0.)
							disk_source_terms[i](j_species) = 0.;
				}
			}
			catch(...)
			{
				Info << "Wrong input data in " << disk_name << " file" << endl;
				abort();
			}
		}

		// Creation source terms
		{
			Info<< "Source terms from disks creation... " << endl;

			const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
			sourceFromDisk.resize(NC);

			for (int i=0;i<NC;i++)
			{
				sourceFromDisk.set
				(
					i,
					new volScalarField
					(
						IOobject
						(
							"SDisk_" + thermodynamicsMapXML->NamesOfSpecies()[i],
							mesh.time().timeName(),
							mesh,
							IOobject::NO_READ,
							IOobject::AUTO_WRITE
						),
						mesh,
						dimensionedScalar("SDisk", dimensionSet(1, -3, -1, 0, 0), 0.0)
					)
				);
			}

			for(unsigned int i=0;i<list_disks.size();i++)
				for(unsigned int j=0;j<disk_topology_indices[i].size();j++)
				{
					const int celli = disk_topology_indices[i](j);

					#if OPENFOAM_VERSION >= 40
					for(int k=0;k<NC;k++)
						sourceFromDisk[k].ref()[celli] = disk_source_terms[i](k);
					#else
					for(int i=0;i<NC+1;i++)
						sourceFromDisk[k].internalField()[celli] = = disk_source_terms[i](k);
					#endif	
				}

			for(int k=0;k<NC;k++)
				sourceFromDisk[k].write();

			Info << "Source terms from disks correctly written!" << endl;
			abort();
		}
	}

	if (diskPreprocess == false)
	{
		Info<< "Source terms from disks creation... " << endl;

		const unsigned int NC = thermodynamicsMapXML->NumberOfSpecies();
		sourceFromDisk.resize(NC);

		for (int i=0;i<NC;i++)
		{
			sourceFromDisk.set
			(
				i,
				new volScalarField
				(
					IOobject
					(
						"SDisk_" + thermodynamicsMapXML->NamesOfSpecies()[i],
						mesh.time().timeName(),
						mesh,
						IOobject::MUST_READ,
						IOobject::AUTO_WRITE
					),
					mesh
				)
			);
		}
	}
}

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

if (iSOM == true)
{
	Info<< "Classification of points (SOM)... " ;

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	{
		Eigen::VectorXd original(nf+1);
		Eigen::VectorXd transformed(npca);
		Eigen::VectorXd yCluster(nclusters);

		#if OPENFOAM_VERSION >= 40
		const scalarField& TCells = T.internalField();
		scalarField& clusterIndexCells = clusterIndex.ref();
		scalarField& clusterWeightCells =  clusterWeight.ref();
		scalarField& clusterSpeciesCells =  clusterSpecies.ref();
		scalarField& clusterReactionsCells = clusterReactions.ref();
		#else
		const scalarField& TCells = T.internalField();
		scalarField& clusterIndexCells = clusterIndex.internalField();
		scalarField& clusterWeightCells =  clusterWeight.internalField();
		scalarField& clusterSpeciesCells =  clusterSpecies.internalField();
		scalarField& clusterReactionsCells = clusterReactions.internalField();
		#endif

		forAll(TCells, celli)
		{
			if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
			{
				// PCA transformation
				original(0) = (TCells[celli]-pca_mu(0))/pca_sigma(0);
				if (nf == thermodynamicsMapXML->NumberOfSpecies())
				{
					for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
						original(i+1) = (Y[i].internalField()[celli]-pca_mu(i+1))/pca_sigma(i+1);
				}
				else
				{
					for(int i=0;i<nf;i++)
						original(i+1) = (Y[listretspecies(i)-1].internalField()[celli]-pca_mu(i+1))/pca_sigma(i+1);
				}

				// Transformation
				transformed = original.transpose()*pca_weights;

				// Recognizing the cluster
				int index = -1;
				myNeuralNetworkBasedOnPCA(transformed.data(), yCluster.data()); 
				clusterWeightCells[celli] = yCluster.maxCoeff(&index);
				clusterIndexCells[celli] = index;
				clusterSpeciesCells[celli] = vector_thermodynamicsMapXML[index]->NumberOfSpecies();
				clusterReactionsCells[celli] = vector_kineticsMapXML[index]->NumberOfReactions();
			}
			else
			{
				// Default values (no classification is applied)
				clusterWeightCells[celli] = 1;
				clusterIndexCells[celli] = 0;		// never used
				clusterSpeciesCells[celli] = 3;		// fuel, oxidizer, inert
				clusterReactionsCells[celli] = 0;	// no reactions
			}
		}
   	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}

if (iVQ2 == true || iSOFTMAX == true)
{
	
	if (iVQ2 == true) 	Info<< "Classification of points (VQ2)... " ;
	else if (iSOFTMAX== true) 	Info<< "Classification of points (SOFTMAX)... " ;

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	{
		#if OPENFOAM_VERSION >= 40
		const scalarField& TCells = T.internalField();
		scalarField& clusterIndexCells = clusterIndex.ref();
		scalarField& clusterWeightCells =  clusterWeight.ref();
		scalarField& clusterSpeciesCells =  clusterSpecies.ref();
		scalarField& clusterReactionsCells = clusterReactions.ref();
		#else
		const scalarField& TCells = T.internalField();
		scalarField& clusterIndexCells = clusterIndex.internalField();
		scalarField& clusterWeightCells =  clusterWeight.internalField();
		scalarField& clusterSpeciesCells =  clusterSpecies.internalField();
		scalarField& clusterReactionsCells = clusterReactions.internalField();
		#endif

		int n_cells = TCells.size();

		arma::Mat<double> observations(n_cells, thermodynamicsMapXML->NumberOfSpecies()+1);
		forAll(TCells, celli)
		{
			observations(celli,0) = TCells[celli];
			for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				observations(celli,i+1) = Y[i].internalField()[celli];
		}

		// Classification
		std::vector<int> idx(n_cells);
		if (iVQ2 == true) 		VQ2.Classify(observations, idx);
		else if (iSOFTMAX == true) 	SOFTMAX.Classify(observations, idx);

		forAll(TCells, celli)
		{
			if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
			{
				int index = idx[celli];
				clusterIndexCells[celli] = index;
				clusterSpeciesCells[celli] = vector_thermodynamicsMapXML[index]->NumberOfSpecies();
				clusterReactionsCells[celli] = vector_kineticsMapXML[index]->NumberOfReactions();
			}
			else
			{
				// Default values (no classification is applied)
				clusterWeightCells[celli] = 1;
				clusterIndexCells[celli] = 0;		// never used
				clusterSpeciesCells[celli] = 3;		// fuel, oxidizer, inert
				clusterReactionsCells[celli] = 0;	// no reactions
			}
		}
   	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}

if (iFitCTree == true)
{
	Info<< "Classification of points (FitCTree)... " ;

	double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

	{
		const int nf = thermodynamicsMapXML->NumberOfSpecies();

		Eigen::VectorXd original(nf+1);
		Eigen::VectorXd transformed(fitctree_npcs_retained);
		Eigen::VectorXd yCluster(nclusters);

		#if OPENFOAM_VERSION >= 40
		const scalarField& TCells = T.internalField();
		scalarField& clusterIndexCells = clusterIndex.ref();
		scalarField& clusterWeightCells =  clusterWeight.ref();
		scalarField& clusterSpeciesCells =  clusterSpecies.ref();
		scalarField& clusterReactionsCells = clusterReactions.ref();
		#else
		const scalarField& TCells = T.internalField();
		scalarField& clusterIndexCells = clusterIndex.internalField();
		scalarField& clusterWeightCells =  clusterWeight.internalField();
		scalarField& clusterSpeciesCells =  clusterSpecies.internalField();
		scalarField& clusterReactionsCells = clusterReactions.internalField();
		#endif

		classifyPoint_initialize();

		forAll(TCells, celli)
		{
			if (TCells[celli] > direct_integration_minimum_temperature_for_chemistry)
			{
				// PCA transformation
				original(0) = (TCells[celli]-fitctree_mu(0))/fitctree_sigma(0);
				if (nf == thermodynamicsMapXML->NumberOfSpecies())
				{
					for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
						original(i+1) = (Y[i].internalField()[celli]-fitctree_mu(i+1))/fitctree_sigma(i+1);
				}

				// Recognizing the cluster
				int index = -1;
				if (fitctree_npcs_retained == 0)
				{
					index = static_cast<int>( classifyPoint(original.data()) ) - 1;
				}
				else
				{
					transformed = original.transpose()*fitctree_pca_coeff;
					index = static_cast<int>( classifyPoint(transformed.data()) ) - 1;
				}

				// DEBUG
				// std::cout << celli << " " << TCells[celli] << " " << index << " " << vector_thermodynamicsMapXML[index]->NumberOfSpecies() << " " << vector_kineticsMapXML[index]->NumberOfReactions() << std::endl;

				clusterWeightCells[celli] = 1;
				clusterIndexCells[celli] = index;
				clusterSpeciesCells[celli] = vector_thermodynamicsMapXML[index]->NumberOfSpecies();
				clusterReactionsCells[celli] = vector_kineticsMapXML[index]->NumberOfReactions();
			}
			else
			{
				// Default values (no classification is applied)
				clusterWeightCells[celli] = 1;
				clusterIndexCells[celli] = 0;		// never used
				clusterSpeciesCells[celli] = 3;		// fuel, oxidizer, inert
				clusterReactionsCells[celli] = 0;	// no reactions
			}
		}

		classifyPoint_terminate();
   	}
	
	double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

	Info << "done in " << tEnd - tStart << " s  (" << (tEnd-tStart)/double(mesh.nCells())*1000. << " ms per cell)" << endl;
}

#endif

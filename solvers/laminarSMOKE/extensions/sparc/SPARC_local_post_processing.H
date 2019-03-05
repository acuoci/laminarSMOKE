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

// SOM variables
if (iSOM == true)
{

	const double som_mean_chemistryCpu =  sum(cpuChemistry).value()/double(mesh.globalData().nTotalCells());
	const double som_max_chemistryCpu = max(cpuChemistry).value();
	const double som_mean_species = sum(clusterSpecies).value()/double(mesh.globalData().nTotalCells());
	const double som_max_species = max(clusterSpecies).value();
        const double som_mean_reactions = sum(clusterReactions).value()/double(mesh.globalData().nTotalCells());
        const double som_max_reactions = max(clusterReactions).value();

	if (somIntegralCounter == somIntegralUpdate)
	{	
		if (Pstream::master())
		{
			fSom() << std::setw(20) << runTime.value();
			fSom() << std::setw(20) << som_mean_species;
			fSom() << std::setw(20) << som_max_species;
			fSom() << std::setw(20) << som_mean_reactions;
			fSom() << std::setw(20) << som_max_reactions;
			fSom() << std::setw(20) << som_mean_chemistryCpu;
			fSom() << std::setw(20) << som_max_chemistryCpu;
			fSom() << std::endl;
		}
		
		somIntegralCounter = 0;
	}

	somIntegralCounter++;
}

// VQ2 variables
if (iVQ2 == true)
{

	const double vq2_mean_chemistryCpu =  sum(cpuChemistry).value()/double(mesh.globalData().nTotalCells());
	const double vq2_max_chemistryCpu = max(cpuChemistry).value();
	const double vq2_mean_species = sum(clusterSpecies).value()/double(mesh.globalData().nTotalCells());
	const double vq2_max_species = max(clusterSpecies).value();
        const double vq2_mean_reactions = sum(clusterReactions).value()/double(mesh.globalData().nTotalCells());
        const double vq2_max_reactions = max(clusterReactions).value();

	if (vq2IntegralCounter == vq2IntegralUpdate)
	{	
		if (Pstream::master())
		{
			fVQ2() << std::setw(20) << runTime.value();
			fVQ2() << std::setw(20) << vq2_mean_species;
			fVQ2() << std::setw(20) << vq2_max_species;
			fVQ2() << std::setw(20) << vq2_mean_reactions;
			fVQ2() << std::setw(20) << vq2_max_reactions;
			fVQ2() << std::setw(20) << vq2_mean_chemistryCpu;
			fVQ2() << std::setw(20) << vq2_max_chemistryCpu;
			fVQ2() << std::endl;
		}
		
		vq2IntegralCounter = 0;
	}

	vq2IntegralCounter++;
}

// SOFTMAX variables
if (iSOFTMAX == true)
{

	const double softmax_mean_chemistryCpu =  sum(cpuChemistry).value()/double(mesh.globalData().nTotalCells());
	const double softmax_max_chemistryCpu = max(cpuChemistry).value();
	const double softmax_mean_species = sum(clusterSpecies).value()/double(mesh.globalData().nTotalCells());
	const double softmax_max_species = max(clusterSpecies).value();
        const double softmax_mean_reactions = sum(clusterReactions).value()/double(mesh.globalData().nTotalCells());
        const double softmax_max_reactions = max(clusterReactions).value();

	if (softmaxIntegralCounter == softmaxIntegralUpdate)
	{	
		if (Pstream::master())
		{
			fSOFTMAX() << std::setw(20) << runTime.value();
			fSOFTMAX() << std::setw(20) << softmax_mean_species;
			fSOFTMAX() << std::setw(20) << softmax_max_species;
			fSOFTMAX() << std::setw(20) << softmax_mean_reactions;
			fSOFTMAX() << std::setw(20) << softmax_max_reactions;
			fSOFTMAX() << std::setw(20) << softmax_mean_chemistryCpu;
			fSOFTMAX() << std::setw(20) << softmax_max_chemistryCpu;
			fSOFTMAX() << std::endl;
		}
		
		softmaxIntegralCounter = 0;
	}

	softmaxIntegralCounter++;
}

#endif

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
|   Authors: A. Cuoci                                                     |
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
|   Copyright(C) 2016, 2015, 2014 A. Cuoci                                |
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

enum species_order_policy_enum {SPECIES_ORDER_POLICY_CONSTANT, SPECIES_ORDER_POLICY_SWEEP, SPECIES_ORDER_POLICY_RANDOM_SHUFFLE, SPECIES_ORDER_POLICY_ROTATE };

label jacobianUpdate = 1;
label propertiesUpdate = 1;
Switch implicitSourceTerm = true;
Switch sparseJacobian  = true;
species_order_policy_enum species_order_policy = SPECIES_ORDER_POLICY_CONSTANT;
std::vector<std::string> exceptional_species;

const dictionary& steadyStateDictionary = solverOptions.subDict("SteadyState");
{
	sparseJacobian     = Switch(steadyStateDictionary.lookup(word("sparseJacobian")));
	jacobianUpdate     = readLabel(steadyStateDictionary.lookup("jacobianUpdate"));
	propertiesUpdate   = readLabel(steadyStateDictionary.lookup("propertiesUpdate"));
	implicitSourceTerm = Switch(steadyStateDictionary.lookup(word("implicitSourceTerm")));

	implicitSourceTerm = Switch(steadyStateDictionary.lookup(word("implicitSourceTerm")));

	word order_policy(steadyStateDictionary.lookup("orderSpecies"));
	if (order_policy == "constant")		species_order_policy = SPECIES_ORDER_POLICY_CONSTANT;
	else if (order_policy == "sweep")	species_order_policy = SPECIES_ORDER_POLICY_SWEEP;
	else if (order_policy == "random")	species_order_policy = SPECIES_ORDER_POLICY_RANDOM_SHUFFLE;
	else if (order_policy == "rotate")	species_order_policy = SPECIES_ORDER_POLICY_ROTATE;
	else
	{
		Info << "Wrong orderPolicy option: constant | sweep | random | rotate" << endl;
		abort();
	}

	// Exceptional species
	exceptional_species.resize(thermodynamicsMapXML->NumberOfSpecies());
	for(unsigned int k=0;k<thermodynamicsMapXML->NumberOfSpecies();k++)
		exceptional_species[k] = "Yi";
	
	List<word> list_exceptional_species = List<word>(steadyStateDictionary.lookup("exceptionalSpecies"));
	if (list_exceptional_species[0] != "none")
	{
		for(unsigned int k=0;k<list_exceptional_species.size();k++)
		{
			const int index = thermodynamicsMapXML->IndexOfSpecies(list_exceptional_species[k])-1;
			exceptional_species[index] = "Yi_" + list_exceptional_species[k];
		}
	}
}

label jacobianCounter = jacobianUpdate;
label propertiesCounter = propertiesUpdate;

std::vector<Eigen::VectorXd> Jstore;
if (jacobianUpdate != 1)
{
	Jstore.resize(mesh.nCells());
	for(int i=0;i<mesh.nCells();i++)
		Jstore[i].resize(thermodynamicsMapXML->NumberOfSpecies()+1);
}

std::vector<int> species_order(thermodynamicsMapXML->NumberOfSpecies());
for(int i=0;i<species_order.size();i++)
{
	species_order[i] = i;
}


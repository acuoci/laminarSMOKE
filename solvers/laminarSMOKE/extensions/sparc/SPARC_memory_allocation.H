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

// Batch reactor homogeneous
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressureSPARC_ODE_OpenSMOKE> denseOdeConstantPressureSPARC;
typedef OdeSMOKE::MethodGear<denseOdeConstantPressureSPARC> methodGearConstantPressureSPARC;

PtrList< BatchReactorHomogeneousConstantPressureSPARC > 		vector_batchReactorHomogeneousConstantPressure(nclusters);
PtrList< BatchReactorHomogeneousConstantVolume >   			vector_batchReactorHomogeneousConstantVolume(nclusters);
PtrList< OdeSMOKE::MultiValueSolver<methodGearConstantPressureSPARC> >  vector_odeSolverConstantPressure(nclusters);
PtrList< OdeSMOKE::MultiValueSolver<methodGearConstantVolume>   >  	vector_odeSolverConstantVolume(nclusters);
for(unsigned int i=0;i<nclusters;i++)
{
	Info << "Creating constant pressure batch reactor: " << i << " - "; 
	vector_batchReactorHomogeneousConstantPressure.set(i, new BatchReactorHomogeneousConstantPressureSPARC(*vector_thermodynamicsMapXML[i], *vector_kineticsMapXML[i]));
	//vector_batchReactorHomogeneousConstantPressure.set(i, new BatchReactorHomogeneousConstantPressureSPARC(*thermodynamicsMapXML, *kineticsMapXML));
	Info << "Number of species: " << vector_batchReactorHomogeneousConstantPressure[i].NumberOfEquations() << endl; 

	Info << "Creating constant volume batch reactor:   " << i << " - "; 
	vector_batchReactorHomogeneousConstantVolume.set(i, new BatchReactorHomogeneousConstantVolume(*vector_thermodynamicsMapXML[i], *vector_kineticsMapXML[i]));
	//vector_batchReactorHomogeneousConstantVolume.set(i, new BatchReactorHomogeneousConstantVolume(*thermodynamicsMapXML, *kineticsMapXML));
	Info << "Number of species: " << vector_batchReactorHomogeneousConstantPressure[i].NumberOfEquations() << endl; 

	Info << "Creating constant pressure ODE solver:   " << i << endl; 
	vector_odeSolverConstantPressure.set(i, new OdeSMOKE::MultiValueSolver<methodGearConstantPressureSPARC>);
	vector_odeSolverConstantPressure[i].SetReactor(&vector_batchReactorHomogeneousConstantPressure[i]);

	Info << "Creating constant volume ODE solver:   " << i << endl; 
	vector_odeSolverConstantVolume.set(i, new OdeSMOKE::MultiValueSolver<methodGearConstantVolume>);
	vector_odeSolverConstantVolume[i].SetReactor(&vector_batchReactorHomogeneousConstantVolume[i]);
}

#endif

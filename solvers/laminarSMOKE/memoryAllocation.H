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

//- Memory allocation: gas-phase properties
OpenSMOKE::OpenSMOKEVectorDouble massFractions(thermodynamicsMapXML->NumberOfSpecies());
OpenSMOKE::OpenSMOKEVectorDouble moleFractions(thermodynamicsMapXML->NumberOfSpecies());

//- Memory allocation: gas-phase chemistry
OpenSMOKE::OpenSMOKEVectorDouble omega(thermodynamicsMapXML->NumberOfSpecies());

#if STEADYSTATE != 1

// Batch reactor homogeneous
BatchReactorHomogeneousConstantPressure batchReactorHomogeneousConstantPressure(*thermodynamicsMapXML, *kineticsMapXML);
BatchReactorHomogeneousConstantVolume   batchReactorHomogeneousConstantVolume(*thermodynamicsMapXML, *kineticsMapXML);

// ODE Solver (constant pressure)
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_OpenSMOKE> denseOdeConstantPressure;
typedef OdeSMOKE::MethodGear<denseOdeConstantPressure> methodGearConstantPressure;
autoPtr< OdeSMOKE::MultiValueSolver<methodGearConstantPressure> > odeSolverConstantPressure;
odeSolverConstantPressure.reset(new OdeSMOKE::MultiValueSolver<methodGearConstantPressure>);
odeSolverConstantPressure().SetReactor(&batchReactorHomogeneousConstantPressure);

// ODE Solver (constant volume)
typedef OdeSMOKE::KernelDense<OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_OpenSMOKE> denseOdeConstantVolume;
typedef OdeSMOKE::MethodGear<denseOdeConstantVolume> methodGearConstantVolume;
autoPtr< OdeSMOKE::MultiValueSolver<methodGearConstantVolume> > odeSolverConstantVolume;
odeSolverConstantVolume.reset(new OdeSMOKE::MultiValueSolver<methodGearConstantVolume>);
odeSolverConstantVolume().SetReactor(&batchReactorHomogeneousConstantVolume);

OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CHEMEQ2 *chemeq2SolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CHEMEQ2)
{
	chemeq2SolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CHEMEQ2::GetInstance();
	chemeq2SolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CHEMEQ2 *chemeq2SolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CHEMEQ2)
{
	chemeq2SolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CHEMEQ2::GetInstance();
	chemeq2SolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}

#if OPENSMOKE_USE_DVODE == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE *dvodeSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DVODE::GetInstance();
	dvodeSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE *dvodeSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
{
	dvodeSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DVODE::GetInstance();
	dvodeSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_ODEPACK == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE *dlsodeSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODE::GetInstance();
	dlsodeSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE *dlsodeSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodeSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODE::GetInstance();
	dlsodeSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA *dlsodaSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DLSODA::GetInstance();
	dlsodaSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA *dlsodaSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE)
{
	dlsodaSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DLSODA::GetInstance();
	dlsodaSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_SUNDIALS == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE *cvodeSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_CVODE::GetInstance();
	cvodeSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE *cvodeSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
{
	cvodeSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_CVODE::GetInstance();
	cvodeSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_DASPK == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK *daspkSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK *daspkSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
{
	daspkSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_DASPK::GetInstance(thermodynamicsMapXML->NumberOfSpecies()+1);
	daspkSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_RADAU == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5 *radauSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radauSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_RADAU5::GetInstance();
	radauSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5 *radauSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
{
	radauSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_RADAU5::GetInstance();
	radauSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#if OPENSMOKE_USE_MEBDF == 1
OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF *mebdfSolverConstantPressure;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfSolverConstantPressure = OpenSMOKE::BatchReactorHomogeneousConstantPressure_ODE_MEBDF::GetInstance();
	mebdfSolverConstantPressure->SetBatchReactorHomogeneousConstantPressure(&batchReactorHomogeneousConstantPressure);
}
OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF *mebdfSolverConstantVolume;
if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
{
	mebdfSolverConstantVolume = OpenSMOKE::BatchReactorHomogeneousConstantVolume_ODE_MEBDF::GetInstance();
	mebdfSolverConstantVolume->SetBatchReactorHomogeneousConstantVolume(&batchReactorHomogeneousConstantVolume);
}
#endif

#endif

#if STEADYSTATE != 1
#if OPENSMOKE_USE_ISAT == 1
    #include "memoryAllocation_ISAT.H"
#endif
#endif



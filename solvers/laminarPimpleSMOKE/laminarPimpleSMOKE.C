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
|                                                                         |
|   Application: laminarPimpleSMOKE                                       |
|                                                                         |
|   Description: unsteady solver for reactive flows with detailed kinetic |
|                mechanisms.                                              |
|                                                                         |
\*-----------------------------------------------------------------------*/

// This is not a steady state simulation
#define STEADYSTATE 0

// OpenSMOKE++ Definitions
#include "OpenSMOKEpp"

// CHEMKIN maps
#include "maps/Maps_CHEMKIN"
#include "utilities/soot/polimi/OpenSMOKE_PolimiSoot_Analyzer.h"

// OpenSMOKE++ Dictionaries
#include "dictionary/OpenSMOKE_Dictionary"

// ODE solvers
#include "math/native-ode-solvers/MultiValueSolver"
#include "math/external-ode-solvers/ODE_Parameters.h"

// OpenFOAM
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "interpolation.H"

// Customized radiation model
#include "OpenSMOKEradiationModel.H"

#if OPENFOAM_VERSION <= 30
#include "fvIOoptionList.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#endif

#if OPENFOAM_VERSION >= 40
#include "fvOptions.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "pressureControl.H"
#endif

// Additional include files
#include "sparkModel.H"
#include "utilities.H"
#include "laminarSMOKEthermoClass.H"

// Homogeneous reactors
#include "DRG.h"
#include "BatchReactorHomogeneousConstantPressure.H"
#include "BatchReactorHomogeneousConstantPressure_ODE_Interface.H"
#include "BatchReactorHomogeneousConstantVolume.H"
#include "BatchReactorHomogeneousConstantVolume_ODE_Interface.H"

// ISAT
#if OPENSMOKE_USE_ISAT == 1
    #include "ISAT.h"
    #include "numericalJacobian4ISAT.H"
    #include "mappingGradients/mappingGradient4OpenFOAM.h"
#endif

// Soot
#include "sootUtilities.H"

template<typename Solver, typename OdeBatch>
void SolveOpenSourceSolvers(OdeBatch& ode, const double t0, const double tf, const OpenSMOKE::OpenSMOKEVectorDouble& y0, OpenSMOKE::OpenSMOKEVectorDouble& yf, const OpenSMOKE::ODE_Parameters& parameters)
{
	Solver o(ode);
	o.SetDimensions(y0.Size());
	o.SetAbsoluteTolerance(parameters.absolute_tolerance());
	o.SetRelativeTolerance(parameters.relative_tolerance());
	o.SetAnalyticalJacobian(false);
	o.SetInitialValues(t0, y0.GetHandle());
	o.Solve(tf);
	o.Solution(yf.GetHandle());
}		

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	#if OPENFOAM_VERSION == 40
		#include "laminarPimpleSMOKE.4x.H"
	#elif OPENFOAM_VERSION == 30
		#include "laminarPimpleSMOKE.3x.H"
	#else
		#include "laminarPimpleSMOKE.2x.H"
	#endif
}


// ************************************************************************* //

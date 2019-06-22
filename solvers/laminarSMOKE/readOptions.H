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

Info<< "Reading Solver options\n" << endl;
IOdictionary solverOptions
(
	IOobject
	(
		"solverOptions",
		U.time().constant(),
		U.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

label inertIndex;
OpenSMOKE::ThermodynamicsMap_CHEMKIN*			thermodynamicsMapXML; 
OpenSMOKE::KineticsMap_CHEMKIN* 			kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN* 		transportMapXML;

const dictionary& kineticsDictionary = solverOptions.subDict("Kinetics");
{
	//- Kinetics folder
	Foam::string kinetics_folder = kineticsDictionary.lookup("folder");

	// Reading homogeneous data
	Info<< " * loading homogeneous kinetics...\n" << endl;
	{	
		word kinetics("opensmoke");
		boost::filesystem::path path_kinetics = kinetics_folder;

		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc, xml_string,path_kinetics / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN(doc); 
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN(doc); 
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN(*thermodynamicsMapXML, doc); 
							
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
	}

	//- Inert species
	word inertSpecies(kineticsDictionary.lookup("inertSpecies"));
	inertIndex = thermodynamicsMapXML->IndexOfSpecies(inertSpecies)-1;
}

//- Mass fractions tolerance
scalar massFractionsTol = 1e-10;

//- Internal models
enum { STRANG_MOMENTUM_TRANSPORT_REACTION, STRANG_MOMENTUM_REACTION_TRANSPORT, STRANG_COMPACT } strangAlgorithm;

// Physical model
Switch energyEquation;
Switch speciesEquations;
Switch momentumEquations;
Switch zMixEquation;
Switch tauEquation;
Switch homogeneousReactions;
Switch iCorrectiveFluxes;
Switch constPressureBatchReactor;
Switch includeDpDt;
Eigen::VectorXd LewisNumbers;
scalar thermodynamicPressure = 101325.;
scalar direct_integration_minimum_temperature_for_chemistry = 200.;

const dictionary& physicalModelDictionary = solverOptions.subDict("PhysicalModel");
{
	momentumEquations = Switch(physicalModelDictionary.lookup(word("momentumEquations")));
	energyEquation = Switch(physicalModelDictionary.lookup(word("energyEquation")));
	speciesEquations = Switch(physicalModelDictionary.lookup(word("speciesEquations")));
	zMixEquation = Switch(physicalModelDictionary.lookup(word("zMixEquation")));
	tauEquation = Switch(physicalModelDictionary.lookup(word("tauEquation")));
	homogeneousReactions = Switch(physicalModelDictionary.lookup(word("homogeneousReactions")));
	includeDpDt = Switch(physicalModelDictionary.lookup(word("includeDpDt")));
	iCorrectiveFluxes = Switch(physicalModelDictionary.lookup(word("correctiveFluxes")));
	direct_integration_minimum_temperature_for_chemistry  = readScalar(physicalModelDictionary.lookup("minTemperature4Chemistry"));
	
	// Algorithm
	{
		word strang(physicalModelDictionary.lookup("strangAlgorithm"));
		if (strang == "MomentumTransportReaction")	strangAlgorithm = STRANG_MOMENTUM_TRANSPORT_REACTION;
		else if (strang == "MomentumReactionTransport")	strangAlgorithm = STRANG_MOMENTUM_REACTION_TRANSPORT;
		else if (strang == "Compact")			strangAlgorithm = STRANG_COMPACT;
		else
		{
			Info << "Wrong strangAlgorithm option: MomentumTransportReaction || MomentumReactionTransport || Compact" << endl;
			abort();
		}
	}
	
	// Diffusivity of species
	{
		LewisNumbers.resize(thermodynamicsMapXML->NumberOfSpecies());
	
		const dictionary& lewisNumberDictionary = solverOptions.subDict("LewisNumbers");
		for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			LewisNumbers(i) = readScalar(lewisNumberDictionary.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
	
		Info << endl;
		Info << "Diffusivities of species will be calculated using the following Lewis numbers: " << endl;
		for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			Info << " " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t" << LewisNumbers(i) << endl;
		Info << endl;
	}
	
	//- Constant pressure reactors
	constPressureBatchReactor = Switch(physicalModelDictionary.lookup(word("constPressureBatchReactor")));

	// Themodynamic pressure [Pa]
	thermodynamicPressure = readScalar(physicalModelDictionary.lookup("thermodynamicPressure"));
}

//- Detect spark
const dictionary& sparkDictionary = solverOptions.subDict("Spark");
sparkModel spark;
spark.Read(sparkDictionary);

Eigen::VectorXd outputFormationRatesIndices;
Eigen::VectorXd outputDiffusivitiesIndices;
IOobject::writeOption outputGasProperties = IOobject::NO_WRITE;
IOobject::writeOption outputGasRadiation = IOobject::NO_WRITE;
const dictionary& outputDictionary = solverOptions.subDict("Output"); 
{
	// Write additional fields containing properties of gaseous phase
	// In order to reduce the IO operation, by default it is set to false
	{
		Switch flag = Switch(outputDictionary.lookup(word("gasProperties")));
		if (flag == true)	outputGasProperties = IOobject::AUTO_WRITE;
		else			outputGasProperties = IOobject::AUTO_WRITE;
	}

	// Write additional fields relevant to the radiative heat transfer
	// In order to reduce the IO operation, by default it is set to false
	{
		Switch flag = Switch(outputDictionary.lookup(word("gasRadiation")));
		if (flag == true)	outputGasRadiation = IOobject::AUTO_WRITE;
		else			outputGasRadiation = IOobject::AUTO_WRITE;
	}

	// Write formation rates of species
	{
		Switch outputFormationRates = Switch(outputDictionary.lookup(word("formationRates")));
		if (outputFormationRates == true)
		{
			List<word>  listFormationRates(outputDictionary.lookup("listFormationRates"));
			outputFormationRatesIndices.resize(listFormationRates.size());
			for (int i=0;i<listFormationRates.size();i++)
				outputFormationRatesIndices(i) = thermodynamicsMapXML->IndexOfSpecies(listFormationRates[i])-1;
		}
	}

	// Write diffusivities of species
	{
		Switch outputDiffusivities = Switch(outputDictionary.lookup(word("diffusivities")));
		if (outputDiffusivities == true)
		{
			List<word>  listDiffusivities(outputDictionary.lookup("listDiffusivities"));
			outputDiffusivitiesIndices.resize(listDiffusivities.size());
			for (int i=0;i<listDiffusivities.size();i++)
				outputDiffusivitiesIndices(i) = thermodynamicsMapXML->IndexOfSpecies(listDiffusivities[i])-1;
		}
	}
}

#if STEADYSTATE != 1

OpenSMOKE::DRG* drg;
Switch drg_analysis = false;
double drg_minimum_temperature_for_chemistry = 300.;
List<double>  drg_epsilon;
List<double>  drg_temperature;
IOobject::writeOption outputDRG = IOobject::NO_WRITE;
const dictionary& drgDictionary = solverOptions.subDict("DRG");
{
	drg_analysis = Switch(drgDictionary.lookup(word("drg")));
	
	if ( drg_analysis == true)
	{
		outputDRG = IOobject::AUTO_WRITE;

		drg_minimum_temperature_for_chemistry  = readScalar(drgDictionary.lookup("minTemperature"));

		std::vector<unsigned int> drgListSpecies;
		List<word>  listSpecies(drgDictionary.lookup("species"));
		drgListSpecies.resize(listSpecies.size());
		for (int i=0;i<drgListSpecies.size();i++)
			drgListSpecies[i] = thermodynamicsMapXML->IndexOfSpecies(listSpecies[i])-1;

		drg_epsilon = readList<double>(drgDictionary.lookup("epsilon"));
		drg_temperature = readList<double>(drgDictionary.lookup("temperature"));

		drg = new OpenSMOKE::DRG(thermodynamicsMapXML, kineticsMapXML);
		drg->SetKeySpecies(drgListSpecies);
	}
}

#endif

/*********************************************************************************************************************/ 

#if STEADYSTATE != 1

scalar chemeq2_epsilonMinimum 	= 0.01;
scalar chemeq2_epsilonMaximum 	= 100.;
scalar chemeq2_yMinimum 	= 1.e-32;
scalar chemeq2_dtMinimum 	= 1.e-20;
label  chemeq2_subIterations 	= 1;

// Batch reactor homogeneous: ode parameters
const dictionary& odeHomogeneousDictionary = solverOptions.subDict("OdeHomogeneous");
OpenSMOKE::ODE_Parameters odeParameterBatchReactorHomogeneous;
{
	//- Mass fractions tolerance
	scalar relTolerance = readScalar(odeHomogeneousDictionary.lookup("relTolerance"));
	scalar absTolerance = readScalar(odeHomogeneousDictionary.lookup("absTolerance"));
	odeParameterBatchReactorHomogeneous.SetRelativeTolerance(relTolerance);
	odeParameterBatchReactorHomogeneous.SetAbsoluteTolerance(absTolerance);
	
	//- Full pivoting (only for OpenSMOKE solver)
	Switch fullPivoting(odeHomogeneousDictionary.lookup("fullPivoting"));
	odeParameterBatchReactorHomogeneous.SetFullPivoting(fullPivoting);
	
	//- Maximum order of integration (only for OpenSMOKE solver)
	label maximumOrder = readLabel(odeHomogeneousDictionary.lookup("maximumOrder"));
	odeParameterBatchReactorHomogeneous.SetMaximumOrder(maximumOrder);
	
	//- 
	
	// Type
	word homogeneousODESolverString(odeHomogeneousDictionary.lookup("odeSolver"));
	if (	homogeneousODESolverString != "OpenSMOKE" 	&& homogeneousODESolverString != "DVODE"  && 
		homogeneousODESolverString != "DLSODE" 		&& homogeneousODESolverString != "DLSODA" && 
		homogeneousODESolverString != "CVODE" 		&& homogeneousODESolverString != "DASPK"  &&
		homogeneousODESolverString != "MEBDF" 		&& homogeneousODESolverString != "RADAU5"  &&
		homogeneousODESolverString != "CHEMEQ2"
	   )
	{
		Info << "Wrong homogeneous ODE Solver: OpenSMOKE || DVODE || DLSODE || DLSODA || CVODE || DASPK || MEBDF || RADAU5" << endl;
		abort();
	}

	if (homogeneousODESolverString == "OpenSMOKE")	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_OPENSMOKE);
	if (homogeneousODESolverString == "DVODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE);
	if (homogeneousODESolverString == "DLSODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE);
	if (homogeneousODESolverString == "DLSODA") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA);
	if (homogeneousODESolverString == "CVODE") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE);
	if (homogeneousODESolverString == "DASPK") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK);
	if (homogeneousODESolverString == "MEBDF") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF);
	if (homogeneousODESolverString == "RADAU5") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5);	
	if (homogeneousODESolverString == "CHEMEQ2") 	odeParameterBatchReactorHomogeneous.SetType(OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CHEMEQ2);	

	if (homogeneousODESolverString == "CHEMEQ2")
	{
		const dictionary& chemeq2Dictionary = odeHomogeneousDictionary.subDict("CHEMEQ2");

		chemeq2_epsilonMinimum 	= readScalar(chemeq2Dictionary.lookup("epsilonMinimum"));
		chemeq2_epsilonMaximum 	= readScalar(chemeq2Dictionary.lookup("epsilonMaximum"));
		chemeq2_yMinimum 	= readScalar(chemeq2Dictionary.lookup("yMinimum"));
		chemeq2_dtMinimum 	= readScalar(chemeq2Dictionary.lookup("dtMinimum"));
		chemeq2_subIterations 	= readLabel(chemeq2Dictionary.lookup("subIterations"));
	}	
}


// Check type
{
	if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DVODE)
	{
		#if OPENSMOKE_USE_DVODE != 1
		{
			Info << "The solver was compiled without the DVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_CVODE)
	{
		#if OPENSMOKE_USE_SUNDIALS != 1
		{
			Info << "The solver was compiled without the CVODE support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if ( odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODE ||
	          odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DLSODA )
	{
		#if OPENSMOKE_USE_ODEPACK != 1
		{
			Info << "The solver was compiled without the ODEPACK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_RADAU5)
	{
		#if OPENSMOKE_USE_RADAU != 1
		{
			Info << "The solver was compiled without the RADAU support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}	
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_DASPK)
	{
		#if OPENSMOKE_USE_DASPK != 1
		{
			Info << "The solver was compiled without the DASPK support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}
	else if (odeParameterBatchReactorHomogeneous.type() == OpenSMOKE::ODE_Parameters::ODE_INTEGRATOR_MEBDF)
	{
		#if OPENSMOKE_USE_MEBDF != 1
		{
			Info << "The solver was compiled without the MEBDF support. Please choose a different ODE solver." << endl;
			abort();
		}
		#endif
	}			
}
#endif

#if STEADYSTATE != 1
#if OPENSMOKE_USE_ISAT == 1
    #include "readOptions_ISAT.H"
#endif
#endif

Info<< "Solver options correctly read\n" << endl;

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
OpenSMOKE::ThermodynamicsMap_CHEMKIN*		thermodynamicsMapXML; 
OpenSMOKE::KineticsMap_CHEMKIN* 		kineticsMapXML;
OpenSMOKE::TransportPropertiesMap_CHEMKIN* 	transportMapXML;

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
enum { STRANG_MOMENTUM_TRANSPORT_REACTION, STRANG_MOMENTUM_REACTION_TRANSPORT} strangAlgorithm;
enum { DIFFUSIVITY_MODEL_MULTICOMPONENT, DIFFUSIVITY_MODEL_LEWIS_NUMBERS} diffusivityModel;

// Physical model
Switch energyEquation;
Switch speciesEquations;
Switch momentumEquations;
Switch zMixEquation;
Switch tauEquation;
Switch homogeneousReactions;
Switch iMassDiffusionInEnergyEquation;
Switch constPressureBatchReactor;
Switch includeDpDt;
Switch soretEffect;
Switch thermophoreticEffect;
Eigen::VectorXd LewisNumbers;
std::vector<bool> soretEffectList;
std::vector<int> indexSoret;
scalar thermodynamicPressure = 101325.;
label zMixDmixIndex = inertIndex;
scalar direct_integration_minimum_temperature_for_chemistry = 200.;
autoPtr<std::ofstream> fSootIntegrals;

const dictionary& physicalModelDictionary = solverOptions.subDict("PhysicalModel");
{
	momentumEquations = Switch(physicalModelDictionary.lookup(word("momentumEquations")));
	energyEquation = Switch(physicalModelDictionary.lookup(word("energyEquation")));
	speciesEquations = Switch(physicalModelDictionary.lookup(word("speciesEquations")));
	zMixEquation = Switch(physicalModelDictionary.lookup(word("zMixEquation")));
	tauEquation = Switch(physicalModelDictionary.lookup(word("tauEquation")));
	homogeneousReactions = Switch(physicalModelDictionary.lookup(word("homogeneousReactions")));
	includeDpDt = Switch(physicalModelDictionary.lookup(word("includeDpDt")));
	iMassDiffusionInEnergyEquation = Switch(physicalModelDictionary.lookup(word("massDiffusionInEnergyEquation")));
	direct_integration_minimum_temperature_for_chemistry  = physicalModelDictionary.lookupOrDefault<scalar>(word("minTemperature4Chemistry"),  scalar(200.));
	
	// Thermophoretic effect
	thermophoreticEffect = Switch(physicalModelDictionary.lookup(word("thermophoreticEffect")));

	// Soret effect
	soretEffect = Switch(physicalModelDictionary.lookup(word("soretEffect")));
	soretEffectList.resize(thermodynamicsMapXML->NumberOfSpecies());
	indexSoret.resize(thermodynamicsMapXML->NumberOfSpecies());

	for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		soretEffectList[i] = false;	
	for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		indexSoret[i] = -1;
	
	// Algorithm
	{
		word strang(physicalModelDictionary.lookup("strangAlgorithm"));
		if (strang == "MomentumTransportReaction")	strangAlgorithm = STRANG_MOMENTUM_TRANSPORT_REACTION;
		else if (strang == "MomentumReactionTransport")	strangAlgorithm = STRANG_MOMENTUM_REACTION_TRANSPORT;
		else
		{
			Info << "Wrong strangAlgorithm option: MomentumTransportReaction || MomentumReactionTransport" << endl;
			abort();
		}
	}
	
	// Diffusivity of species
	{
		word diffusivity(physicalModelDictionary.lookup("diffusivityModel"));
		if (diffusivity == "multi-component")		diffusivityModel = DIFFUSIVITY_MODEL_MULTICOMPONENT;
		else if (diffusivity == "lewis-numbers")	diffusivityModel = DIFFUSIVITY_MODEL_LEWIS_NUMBERS;
		else
		{
			Info << "Wrong diffusivityModel option: multi-component || lewis-numbers" << endl;
			abort();
		}
		
		if (diffusivityModel == DIFFUSIVITY_MODEL_LEWIS_NUMBERS)
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
	}
	
	//- Constant pressure reactors
	constPressureBatchReactor = Switch(physicalModelDictionary.lookup(word("constPressureBatchReactor")));

	// Themodynamic pressure [Pa]
	thermodynamicPressure = readScalar(physicalModelDictionary.lookup("thermodynamicPressure"));

	//- Inert species
	if (zMixEquation == true)
	{
		word referenceSpecies(physicalModelDictionary.lookup("zMixDmixSpecies"));
		zMixDmixIndex = thermodynamicsMapXML->IndexOfSpecies(referenceSpecies)-1;
	}
}

Switch hmom_analysis = false;
Switch hmom_analysis_onthefly = true;
Switch hmom_test = false;
Switch hmom_pah_gas_consumption = true;
int hmom_index_H = -1;
int hmom_index_OH = -1;
int hmom_index_H2O = -1;
int hmom_index_H2 = -1;
int hmom_index_O2 = -1;
int hmom_index_C2H2 = -1;
std::vector<int> hmom_pah_species_indices;
std::vector<std::string> hmom_pah_species_names;

int hmom_pah_nc = 16;
int hmom_number_moments = 4;
int hmom_collision_diameter_model = 2;
int hmom_fractal_diameter_model = 1;

int hmom_nucleation = 1;
int hmom_surface_growth = 1;
int hmom_oxidation = 1;
int hmom_condensation = 1;
int hmom_coagulation = 1;
int hmom_coagulation_continous = 1;
double hmom_schmidt_number = 50.;

scalar hmom_minimum_normalized_value = 1.e-16;
scalar hmom_minimum_fv_for_integration = 1.e-16;


const dictionary& hmomDictionary = solverOptions.subDict("HMOM");
{
	hmom_analysis = Switch(hmomDictionary.lookup(word("hmom")));
	hmom_test = Switch(hmomDictionary.lookup(word("test")));

	if ( hmom_analysis == true)
	{
		hmom_index_H = thermodynamicsMapXML->IndexOfSpecies("H")-1;
		hmom_index_OH = thermodynamicsMapXML->IndexOfSpecies("OH")-1;
		hmom_index_H2O = thermodynamicsMapXML->IndexOfSpecies("H2O")-1;
		hmom_index_H2 = thermodynamicsMapXML->IndexOfSpecies("H2")-1;
		hmom_index_O2 = thermodynamicsMapXML->IndexOfSpecies("O2")-1;
		hmom_index_C2H2 = thermodynamicsMapXML->IndexOfSpecies("C2H2")-1;
		
		List<word>  listSpecies(hmomDictionary.lookup("pah_species"));
		hmom_pah_species_indices.resize(listSpecies.size());
		hmom_pah_species_names.resize(listSpecies.size());
		for (int i=0;i<listSpecies.size();i++)
		{
			hmom_pah_species_indices[i] = thermodynamicsMapXML->IndexOfSpecies(listSpecies[i])-1;
			hmom_pah_species_names[i] = listSpecies[i];
		}

		hmom_pah_nc = readLabel(hmomDictionary.lookup("pah_nc"));
		hmom_collision_diameter_model = readLabel(hmomDictionary.lookup("collision_diameter_model"));
		hmom_fractal_diameter_model = readLabel(hmomDictionary.lookup("fractal_diameter_model"));

		hmom_nucleation = readLabel(hmomDictionary.lookup("nucleation"));
		hmom_surface_growth = readLabel(hmomDictionary.lookup("surface_growth"));
		hmom_oxidation = readLabel(hmomDictionary.lookup("oxidation"));
		hmom_condensation = readLabel(hmomDictionary.lookup("condensation"));
		hmom_coagulation = readLabel(hmomDictionary.lookup("coagulation"));
		hmom_coagulation_continous = readLabel(hmomDictionary.lookup("coagulation_continous"));

		hmom_pah_gas_consumption = Switch(hmomDictionary.lookup(word("pah_gas_consumption")));

		hmom_schmidt_number = hmomDictionary.lookupOrDefault<scalar>(word("schmidt_number"),  scalar(50.));

		hmom_minimum_normalized_value   = hmomDictionary.lookupOrDefault<scalar>(word("mininum_normalized_value"),  scalar(1.e-16));
		hmom_minimum_fv_for_integration = hmomDictionary.lookupOrDefault<scalar>(word("mininum_fv_for_integration"),  scalar(1.e-16));

		Info << "HMOM: minimum normalized moments: " << hmom_minimum_normalized_value << endl;
		Info << "HMOM: minimum fv for integration: " << hmom_minimum_fv_for_integration << endl;

		// Soot integrals
		if ( !boost::filesystem::exists("soot.integrals") )
		{
			fSootIntegrals.reset(new std::ofstream("soot.integrals", std::ios::out));
			fSootIntegrals().setf(std::ios::scientific);

			fSootIntegrals() << std::setw(20) << "time";
			fSootIntegrals() << std::setw(20) << "soot_mass[kg]";
			fSootIntegrals() << std::setw(20) << "soot_volume[m3]";
			fSootIntegrals() << std::setw(20) << "soot_particles[#]";		
			fSootIntegrals() << std::setw(20) << "pah_mass[kg]";
			fSootIntegrals() << std::endl;
		}
		else
		{
			fSootIntegrals.reset(new std::ofstream("soot.integrals", std::ios::app));
			fSootIntegrals().setf(std::ios::scientific);
		}
	}
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

// Thermophoretic effect
std::vector<bool> thermophoreticEffectList;
Switch physicalSootDiffusivity = false;
label  physicalSootDiffusivityReferenceIndex = -1;
std::vector<double> physicalSootDiffusivityCorrection;
std::vector<int> physicalSootDiffusivityCorrectionIndex;
{ 
	const dictionary& polimiSootDictionary = solverOptions.subDict("PolimiSoot");

	if (thermophoreticEffect == true)
	{
		thermophoreticEffectList.resize(thermodynamicsMapXML->NumberOfSpecies());
		for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			thermophoreticEffectList[i] = false;

		// Populate the list of thermophoretic species
		if (thermophoreticEffect == true)
		{
			word tag = polimiSootDictionary.lookup(word("thermophoreticEffectTag"));	
	
			unsigned int count = 0;
			for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			{
				std::string name = thermodynamicsMapXML->NamesOfSpecies()[i];
				std::string prefix = name.substr(0,tag.size());
				if (prefix == tag)
				{
					thermophoreticEffectList[i] = true;
					count++;
				}
			}
				
			if (count == 0)
				std::cout << "Warning: No " << tag << " species were found in the kinetic scheme! Are you sure you want to account for the thermophoretic effect?" << std::endl;
			else
				std::cout << "Thermophoretic effect: Number of " << tag << " species found equal to " << count << std::endl;
		}
	}

	physicalSootDiffusivity = Switch(polimiSootDictionary.lookup(word("physicalSootDiffusivity")));

	if (physicalSootDiffusivity == true)
	{
		//scalar physicalSootDiffusivityExponent = readScalar(polimiSootDictionary.lookup("exponent"));
		scalar physicalSootDiffusivityExponent = -0.681;
		label physicalSootDiffusivityIndex = readLabel(polimiSootDictionary.lookup("physicalSootDiffusivityIndex"));

		for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		{
			if (thermodynamicsMapXML->NamesOfSpecies()[i] == "BIN1A")
			{
				physicalSootDiffusivityReferenceIndex = i;
				break;
			}
		}
		if (physicalSootDiffusivityReferenceIndex == -1)
		{
			Info << "Species BIN1A is not included in the kinetic mechanism" << endl;
			abort();
		}
		else
		{
			int count = -1;
			for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
			{
				if (thermodynamicsMapXML->NamesOfSpecies()[i].compare(0,3,"BIN") == 0)
				{
					count++;
					double MWratio = std::min(	thermodynamicsMapXML->MW()[i+1]/thermodynamicsMapXML->MW()[physicalSootDiffusivityReferenceIndex+1], 
									std::pow(2., double(physicalSootDiffusivityIndex)) );
					physicalSootDiffusivityCorrection.push_back( std::pow(MWratio, physicalSootDiffusivityExponent) );
					physicalSootDiffusivityCorrectionIndex.push_back(i);
			
					Info << "Physical diffusion coefficient: " << thermodynamicsMapXML->NamesOfSpecies()[i] << " : " << physicalSootDiffusivityCorrection[count] << endl;
				}
			}
		}
	}
} 

// Soot post-processing
OpenSMOKE::PolimiSoot_Analyzer* sootAnalyzer;
Switch iPolimiSoot = false;
label sootIntegralUpdate  = 1;
label sootIntegralCounter = sootIntegralUpdate;
Switch sootOnTheFlyFill = false;


//List<vector> pnts_soot_psdf;
{	
		const dictionary& postProcessingDictionary = solverOptions.subDict("PostProcessing");		
		iPolimiSoot = Switch(postProcessingDictionary.lookup(word("soot")));

		if (iPolimiSoot == true)
		{
			const dictionary& postProcessingPolimiSootDictionary = postProcessingDictionary.subDict("PolimiSoot");
		/*	List<word> list_soot_precursors = List<word>(postProcessingPolimiSootDictionary.lookup("sootPrecursors"));

			for(unsigned int k=0;k<list_soot_precursors.size();k++)
			{
				bool found = false;
				for(unsigned int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
				{
					if (list_soot_precursors[k] == thermodynamicsMapXML->NamesOfSpecies()[i])
					{
						soot_precursors_indices.push_back(i);
						found = true;
					}	
				}
				if (found == false)
				{
					Info << "The following soot precursor is not included in the kinetic mechanism: " << list_soot_precursors[k] << endl;
					abort();
				}
				else
				{
					Info << "Soot precursor " << k+1 << " " << list_soot_precursors[k] << " : index " << soot_precursors_indices[k] << endl;
				}
			}
		*/
			
			sootOnTheFlyFill    = Switch(postProcessingPolimiSootDictionary.lookup(word("onTheFlyFill")));
			sootIntegralUpdate  = readLabel(postProcessingPolimiSootDictionary.lookup("integralUpdate"));
			sootIntegralCounter = sootIntegralUpdate;

			Foam::string minimum_bin = postProcessingPolimiSootDictionary.lookup("binMinimum");
			label bin_index_zero     = readLabel(postProcessingPolimiSootDictionary.lookup("binIndexZero"));
			label bin_index_final    = readLabel(postProcessingPolimiSootDictionary.lookup("binIndexFinal"));
			scalar bin_density_zero  = readScalar(postProcessingPolimiSootDictionary.lookup("binDensityZero"));
			scalar bin_density_final = readScalar(postProcessingPolimiSootDictionary.lookup("binDensityFinal"));
			scalar fractal_diameter  = readScalar(postProcessingPolimiSootDictionary.lookup("fractalDiameter"));
			Foam::string planck_coefficient = postProcessingPolimiSootDictionary.lookup("planckCoefficient");



			sootAnalyzer = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);
			sootAnalyzer->SetFractalDiameter(fractal_diameter);
			sootAnalyzer->SetMinimumSection(minimum_bin);
			sootAnalyzer->SetDensity(bin_index_zero, bin_index_final, bin_density_zero, bin_density_final);
			sootAnalyzer->SetPlanckAbsorptionCoefficient(planck_coefficient);
			sootAnalyzer->Setup();

			// Particle size distribution function
			// List<vector> pnts_soot_psdf_dummy(postProcessingPolimiSootDictionary.lookup("PSDF"));
			// pnts_soot_psdf = pnts_soot_psdf_dummy;

			// Soot integrals
			if ( !boost::filesystem::exists("soot.integrals") )
			{
  				fSootIntegrals.reset(new std::ofstream("soot.integrals", std::ios::out));
				fSootIntegrals().setf(std::ios::scientific);

				fSootIntegrals() << std::setw(20) << "time";
				fSootIntegrals() << std::setw(20) << "soot_mass[kg]";
				fSootIntegrals() << std::setw(20) << "soot_volume[m3]";
				fSootIntegrals() << std::setw(20) << "soot_particles[#]";		
				fSootIntegrals() << std::setw(20) << "prec_mass[kg]";
				fSootIntegrals() << std::setw(20) << "prec_volume[m3]";
				fSootIntegrals() << std::setw(20) << "prec_partciles[#]";
				fSootIntegrals() << std::setw(20) << "dimer_mass[kg]";		
				fSootIntegrals() << std::setw(20) << "pah12_340nm[kg]";
				fSootIntegrals() << std::setw(20) << "pah34_400nm[kg]";
				fSootIntegrals() << std::setw(20) << "pahlarge_500nm[kg]";
				fSootIntegrals() << std::endl;
			}
			else
			{
				fSootIntegrals.reset(new std::ofstream("soot.integrals", std::ios::app));
				fSootIntegrals().setf(std::ios::scientific);
			}
		}
}

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

//label minMaxUpdate  = 25;
//label minMaxCounter = minMaxUpdate;

unsigned int vc_main_species = 0;

Info<< "Solver options correctly read\n" << endl;

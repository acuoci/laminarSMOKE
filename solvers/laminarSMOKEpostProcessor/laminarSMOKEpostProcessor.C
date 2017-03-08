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
|   Application: laminarSMOKEpostProcessor                                |
|                                                                         |
|   Description: sample field data with a choice of interpolation schemes |
|                sampling options and write formats.                      |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|    Keywords:                                                            |
|                                                                         |
|    \param setFormat : set output format, choice of \n                   |
|      - xmgr                                                             |
|      - jplot                                                            |
|      - gnuplot                                                          |
|      - raw                                                              |
|                                                                         |
|    \param surfaceFormat : surface output format, choice of \n           |
|      - null        : suppress output                                    |
|      - foamFile    : separate points, faces and values file             |
|      - dx          : DX scalar or vector format                         |
|      - vtk         : VTK ascii format                                   |
|      - raw         : xyz value format for use with e.g. gnuplot splot   |
|      - obj         : Wavefron stl. Does not contain values!             |
|      - stl         : ascii stl. Does not contain values!                |
|                                                                         |
|    \param interpolationScheme : interpolation scheme, choice of \n      |
|      - cell          : use cell-centre value, constant over cells (def) |
|      - cellPoint     : use cell-centre and vertex values                |
|      - cellPointFace : use cell-centre, vertex and face values. \n      |
|        -# vertex values determined from neighbouring cell-centre values |
|        -# face values determined using the current face interpolation   |
|           scheme for the field (linear, limitedLinear, etc.)            |
|                                                                         |
|    \param fields : list of fields to sample                             |
|                                                                         |
|    \param sets : list of sets to sample, choice of \n                   |
|      - uniform             evenly distributed points on line            |
|      - face                one point per face intersection              |
|      - midPoint            one point per cell, inbetween two            |
|                            face intersections                           |
|      - midPointAndFace     combination of face and midPoint             |
|                                                                         |
|      - curve               specified points, not nessecary on line,     |
|                            uses tracking                                |
|      - cloud               specified points, uses findCell              |
|                                                                         |
|        Option axis: how to write point coordinate. Choice of            |
|          - x/y/z: x/y/z coordinate only                                 |
|          - xyz: three columns                                           |
|            (probably does not make sense for anything but raw)          |
|          - distance: from start of sampling line (if uses line)         |
|            or distance from first specified sampling point              |
|                                                                         |
|        Type specific options:                                           |
|            uniform, face, midPoint, midPointAndFace : start and end     |
|            uniform: extra number of sampling points                     |
|            curve, cloud: list of coordinates                            |
|                                                                         |
|    \param surfaces : list of surfaces to sample, choice of \n           |
|      - plane : values on plane defined by point, normal.                |
|      - patch : values on patch.                                         |
|                                                                         |
\*-----------------------------------------------------------------------*/

// OpenSMOKE
#include "OpenSMOKE_Definitions.h"
#include <string>
#include <iostream>
#include <numeric>
#include <Eigen/Dense>

// Base classes
#include "thermo/ThermoPolicy_CHEMKIN.h"
#include "kinetics/ReactionPolicy_CHEMKIN.h"
#include "math/PhysicalConstants.h"
#include "math/OpenSMOKEUtilities.h"

// Maps
#include "maps/ThermodynamicsMap_CHEMKIN.h"
#include "maps/TransportPropertiesMap_CHEMKIN.h"
#include "maps/KineticsMap_CHEMKIN.h"

// OpenFOAM
#include "argList.H"
#include "timeSelector.H"
#include "IOsampledSets.H"
#include "IOsampledSurfaces.H"
#include "fvCFD.H"
#include "multivariateScheme.H"
#include "interpolation.H"

// Soot
#include "sootUtilities.H"
#include "soot/OpenSMOKE_PolimiSoot_Analyzer.h"
#include "SootClassesReader.h"

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    	timeSelector::addOptions();
    	#include "addRegionOption.H"
    	#include "addDictOption.H"
    	#include "setRootCase.H"
    	#include "createTime.H"
    	instantList timeDirs = timeSelector::select0(runTime, args);
    	#include "createNamedMesh.H"

	const word solverOptionsDictionaryName("solverOptions");

	Info<< "Reading field U\n" << endl;
	volVectorField U
	(
  		IOobject
   	 	(
        	"U",
        	runTime.timeName(),
        	mesh,
        	IOobject::MUST_READ,
        	IOobject::AUTO_WRITE
    		),
    		mesh
	);

	Info<< "Reading solverOptions dictionary\n" << endl;
	IOdictionary solverOptionsDictionary
	(
		IOobject
		(
			solverOptionsDictionaryName,
			U.time().constant(),
			U.db(),
			IOobject::MUST_READ,
			IOobject::NO_WRITE
		)
	);

	// Read the kinetic scheme in XML format
	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>* thermodynamicsMapXML; 
	OpenSMOKE::KineticsMap_CHEMKIN<double>* kineticsMapXML;
	OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>* transportMapXML;
	OpenSMOKE::PolimiSoot_Analyzer* sootAnalyzer;
	
	Foam::string kinetics_folder;

	{	
		const dictionary& kineticsDictionary = solverOptionsDictionary.subDict("Kinetics");
		kinetics_folder = kineticsDictionary.lookup("folder");
		boost::filesystem::path path_kinetics = kinetics_folder;

		rapidxml::xml_document<> doc;
		std::vector<char> xml_string;
		OpenSMOKE::OpenInputFileXML(doc,xml_string,path_kinetics / "kinetics.xml");

		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();
		thermodynamicsMapXML = new OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>(doc); 
		transportMapXML = new OpenSMOKE::TransportPropertiesMap_CHEMKIN<double>(doc); 
		kineticsMapXML = new OpenSMOKE::KineticsMap_CHEMKIN<double>(*thermodynamicsMapXML, doc); 					
		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();
		std::cout << " * Time to read XML file: " << tEnd-tStart << std::endl;
	}

	bool calculateThermophoreticVelocity = false;
	bool calculateRatesAcrossBoundaries = false;
	bool calculateLocalStrainRate = false;
	bool calculateMoleFractions = false;
	bool calculateConcentrations = false;
	bool reconstructMixtureFraction = false;
	bool xmlProbeLocations = false;
	bool exportDisks = false;

	std::vector<std::string> 	fuel_names;
	std::vector<std::string> 	oxidizer_names;
	std::vector<double> 		mass_fuel;
	std::vector<double> 		mass_oxidizer;

	List<vector> pnts_xml;

	bool postProcessingPolimiSoot = false;
	List<word> polimiSootBoundaries;
	std::vector<int> soot_precursors_indices;
	SootClassesReader soot_classes_reader;
	bool calculateSootClasses = false;

	boost::filesystem::path soot_folder("sootpp");
	boost::filesystem::create_directory(soot_folder);

	OFstream fSootFvLarge( (soot_folder / "soot_fv_large").string() );
	OFstream fSootFvSmall( (soot_folder / "soot_fv_small").string() );
	OFstream fSootRhoLarge( (soot_folder / "soot_rho_large").string() );
	OFstream fSootRhoSmall( (soot_folder / "soot_rho_small").string() );
	OFstream fSootNLarge( (soot_folder / "soot_N_large").string() );
	OFstream fSootNSmall( (soot_folder / "soot_N_small").string() );	
	OFstream fSootOmegaLarge( (soot_folder / "soot_omega_large").string() );
	OFstream fSootOmegaSmall( (soot_folder / "soot_omega_small").string() );
	OFstream fSootRLarge( (soot_folder / "soot_R_large").string() );
	OFstream fSootRSmall( (soot_folder / "soot_R_small").string() );	
	OFstream fSootR12( (soot_folder / "soot_R_pah_1_2").string() );
	OFstream fSootR34( (soot_folder / "soot_R_pah_3_4").string() );
	OFstream fSootRmore4( (soot_folder / "soot_R_pah_more_4").string() );

	PtrList<OFstream> fsootClassesIntegrals;

	List<vector> pnts_soot_psdf;
	{	
		const dictionary& postProcessingDictionary = solverOptionsDictionary.subDict("PostProcessing");
		
		calculateMoleFractions = Switch(postProcessingDictionary.lookup(word("moleFractions")));
		calculateConcentrations = Switch(postProcessingDictionary.lookup(word("concentrations")));
		calculateThermophoreticVelocity = Switch(postProcessingDictionary.lookup(word("thermophoreticVelocity")));
		calculateRatesAcrossBoundaries = Switch(postProcessingDictionary.lookup(word("flowRates")));
		calculateLocalStrainRate = Switch(postProcessingDictionary.lookup(word("localStrainRate")));

		reconstructMixtureFraction = Switch(postProcessingDictionary.lookup(word("reconstructMixtureFraction")));
		postProcessingPolimiSoot = Switch(postProcessingDictionary.lookup(word("soot")));
		xmlProbeLocations = Switch(postProcessingDictionary.lookup(word("xmlProbeLocations")));

		exportDisks = Switch(postProcessingDictionary.lookup(word("exportDisks")));

		if (xmlProbeLocations == true)
		{
			const dictionary& xmlProbeLocationsDictionary = postProcessingDictionary.subDict("XMLProbeLocations");

			List<vector> pnts_xml_dummy(xmlProbeLocationsDictionary.lookup("Points"));
			pnts_xml = pnts_xml_dummy;
		}

		if (postProcessingPolimiSoot == true)
		{
			const dictionary& postProcessingPolimiSootDictionary = postProcessingDictionary.subDict("PolimiSoot");
			polimiSootBoundaries = List<word>(postProcessingPolimiSootDictionary.lookup("boundaries"));
			Foam::string minimum_bin = postProcessingPolimiSootDictionary.lookup("binMinimum");
			label bin_index_zero     = readLabel(postProcessingPolimiSootDictionary.lookup("binIndexZero"));
			label bin_index_final    = readLabel(postProcessingPolimiSootDictionary.lookup("binIndexFinal"));
			scalar bin_density_zero  = readScalar(postProcessingPolimiSootDictionary.lookup("binDensityZero"));
			scalar bin_density_final = readScalar(postProcessingPolimiSootDictionary.lookup("binDensityFinal"));
			scalar fractal_diameter = readScalar(postProcessingPolimiSootDictionary.lookup("fractalDiameter"));
			calculateSootClasses = Switch(postProcessingPolimiSootDictionary.lookup(word("sootClasses")));
			

			sootAnalyzer = new OpenSMOKE::PolimiSoot_Analyzer(thermodynamicsMapXML);
			sootAnalyzer->SetFractalDiameter(fractal_diameter);
			sootAnalyzer->SetMinimumSection(minimum_bin);
			sootAnalyzer->SetDensity(bin_index_zero, bin_index_final, bin_density_zero, bin_density_final);
			sootAnalyzer->Setup();

			// Particle size distribution function
			List<vector> pnts_soot_psdf_dummy(postProcessingPolimiSootDictionary.lookup("PSDF"));
			pnts_soot_psdf = pnts_soot_psdf_dummy;

			// Soot classes reader
			Info << "Calculate soot classes (new): " << calculateSootClasses << endl;
			if (calculateSootClasses == true)
			{
				boost::filesystem::path path_to_classes_definition = kinetics_folder + "ReacOfProcess.xml";

				if (boost::filesystem::exists(path_to_classes_definition) == false)
				{
					Info << "XML File " << path_to_classes_definition.string() << " does not exist!" << endl;
					abort();
				}

				rapidxml::xml_document<> doc;
				std::vector<char> xml_string;
				OpenSMOKE::OpenInputFileXML(doc, xml_string, path_to_classes_definition);

				soot_classes_reader.ReadFromFile(doc);

				fsootClassesIntegrals.setSize(soot_classes_reader.number_of_classes());
				for (int i=0;i<soot_classes_reader.number_of_classes();i++)
					fsootClassesIntegrals.set( i, new OFstream( (soot_folder / soot_classes_reader.class_name(i)).string()) ) ;
			}
		}

		if (reconstructMixtureFraction == true)
		{
			const dictionary& reconstructMixtureFractionDictionary = postProcessingDictionary.subDict("MixtureFraction");
			List<word>  listFuelNames(reconstructMixtureFractionDictionary.lookup("fuelNames"));
			List<word>  listOxidizerNames(reconstructMixtureFractionDictionary.lookup("oxidizerNames"));
			List<scalar>  listFuelMasses(reconstructMixtureFractionDictionary.lookup("fuelMassFractions"));
			List<scalar>  listOxidizerMasses(reconstructMixtureFractionDictionary.lookup("oxidizerMassFractions"));

			fuel_names.resize(listFuelNames.size());
			oxidizer_names.resize(listOxidizerNames.size());
			mass_fuel.resize(listFuelNames.size());
			mass_oxidizer.resize(listOxidizerNames.size());

			for(unsigned int j=0;j<listFuelNames.size();j++)
			{
				fuel_names[j] = listFuelNames[j];
				mass_fuel[j] = listFuelMasses[j];
			}

			for(unsigned int j=0;j<listOxidizerNames.size();j++)
			{
				oxidizer_names[j] = listOxidizerNames[j];
				mass_oxidizer[j] = listOxidizerMasses[j];
			}
		}
	}
	
	// Formation rates of selected species
	Eigen::VectorXd outputFormationRatesIndices;
	{
		const dictionary& outputDictionary = solverOptionsDictionary.subDict("Output"); 
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
	}

    	forAll(timeDirs, timeI)
    	{
       		runTime.setTime(timeDirs[timeI], timeI);
        	Info<< "Time = " << runTime.timeName() << endl;

        	// Handle geometry/topology changes
        	polyMesh::readUpdateState state = mesh.readUpdate();

		// Read basic fields
		#include "readBasicFields.H"
		#include "readSpecies.H"
		#include "calculateDensity.H"
		#include "compressibleCreatePhi.H"
		#include "calculateEnthalpy.H"
		#include "calculateElementsMassFractions.H"

		// Calculation of rates
		#include "calculateMassFlowRates.H"
		#include "calculateEnthalpyRates.H"
		#include "calculateElementsRates.H"

		// Gas phase
		#include "calculateMoleFractions.H"
		#include "calculateConcentrations.H"
		#include "calculateFormationRates.H"

		// XML probe locations
		#include "calculateProbeLocationsXML.H"

		// Soot
		#include "calculateThermophoreticVelocity.H"
		#include "postProcessingPolimiSoot.H"

		// Reconstructions
	 	#include "postProcessingMixtureFraction.H"

		// Local strain rate
		#include "calculateLocalStrainRate.H"

		// Export Disks
		#include "exportDisks.H"

        	Info<< endl;
    	}

    	Info<< "End\n" << endl;

    	return 0;
}


// ************************************************************************* //

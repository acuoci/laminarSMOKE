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

//- Reading ISAT parameters
Switch isatCheck(solverOptions.subDict("ISAT").lookup("ISAT"));
ISAT *isat_HOM;
Eigen::VectorXd scalingFactors_ISAT;
int luSolver_ISAT = 1;
label numberSubSteps_ISAT = 1;

if (isatCheck == true)
{
	scalar epsilon_ISAT = solverOptions.subDict("ISAT").lookupOrDefault<double>("tolerance", 1e-4);
	       numberSubSteps_ISAT = solverOptions.subDict("ISAT").lookupOrDefault<int>("numberSubSteps", 1);

	scalar maxSizeMRU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSizeMRU", 100);
	scalar maxSizeMFU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSizeMFU", 100);
	scalar maxSizeBT = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSizeBT", 100000);
	scalar maxSearchMRU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSearchMRU", 10);
	scalar maxSearchMFU = solverOptions.subDict("ISAT").lookupOrDefault<int>("maxSearchMFU", 30);
	Switch searchMRU(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("searchMRU", true)); 
	Switch searchMFU(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("searchMFU", true));
	Switch searchBF(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("searchBruteForce", false));
	Switch clearIfFull(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("clearingIfFull", false));
	Switch cleanAndBalance(solverOptions.subDict("ISAT").lookupOrDefault<Switch>("cleanAndBalance", true));
	scalar maxTimeOldCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxTimeOldCoeff", 0.7);
	scalar maxGrowCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxGrowCoeff", 0.5);
	scalar maxHeightCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("maxHeightCoeff", 20.);
	scalar minUsedCoeff = solverOptions.subDict("ISAT").lookupOrDefault<double>("minUsedCoeff", 0.01);
	scalar balanceFactorRetrieve = solverOptions.subDict("ISAT").lookupOrDefault<double>("balanceFactorRetrieve", 2.);
	scalar balanceFactorAddition = solverOptions.subDict("ISAT").lookupOrDefault<double>("balanceFactorAddition", 0.1);
	word   luFactorization = solverOptions.subDict("ISAT").lookupOrDefault<word>("luFactorization","Partial");
	word   qrFactorization = solverOptions.subDict("ISAT").lookupOrDefault<word>("qrFactorization","Full");

	if (luFactorization != "Partial" && luFactorization != "Full")
	{
		Info << "Wrong luFactorization options: Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (luFactorization == "Full") 
			luSolver_ISAT = 0;
		else
			luSolver_ISAT = 1;		
	}

	int qrSolver_ISAT = 0;
	if (qrFactorization != "NoPivoting" && qrFactorization != "Partial" && qrFactorization != "Full")
	{
		Info << "Wrong qrFactorization options: NoPivoting || Partial || Full" << endl;
		abort();
	} 
	else 
	{
		if (qrFactorization == "Full") 
			qrSolver_ISAT = 0;
		else if (qrFactorization == "Partial") 
			qrSolver_ISAT = 1;		
		else if (qrFactorization == "NoPivoting")
			qrSolver_ISAT = 2;
	}


	//- ISAT scale factor 
	scalingFactors_ISAT.setConstant(thermodynamicsMapXML->NumberOfSpecies()+2,1.);
	if(isatCheck == true) 
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();

		massFractionsTol = epsilon_ISAT;
	
		dictionary scalingFactors(solverOptions.subDict("ISAT").subDict("scalingFactors"));
		for(unsigned int i=0;i<NC;i++) 
		{	
			if(scalingFactors.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingFactors_ISAT(i) = 1./readScalar(scalingFactors.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingFactors_ISAT(i) = 1./readScalar(scalingFactors.lookup("others"));
		}
		scalingFactors_ISAT(NC) = 1./readScalar(scalingFactors.lookup("T"));
		scalingFactors_ISAT(NC+1) = 1./readScalar(scalingFactors.lookup("tau"));
	}
	//- ISAT scaling error factors
	Eigen::VectorXd scalingErrors_ISAT;
	scalingErrors_ISAT.setConstant(thermodynamicsMapXML->NumberOfSpecies()+2,1.);
	if(isatCheck == true) 
	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
	
		dictionary scalingErrors(solverOptions.subDict("ISAT").subDict("scalingErrors"));
	
		for(unsigned int i=0;i<NC;i++) 
		{
			if(scalingErrors.found(thermodynamicsMapXML->NamesOfSpecies()[i])) 
				scalingErrors_ISAT(i) = readScalar(scalingErrors.lookup(thermodynamicsMapXML->NamesOfSpecies()[i]));
			else 
				scalingErrors_ISAT(i) = readScalar(scalingErrors.lookup("others"));
		}
		scalingErrors_ISAT(NC) = readScalar(scalingErrors.lookup("T"));
		scalingErrors_ISAT(NC+1) = readScalar(scalingErrors.lookup("tau"));
	}

	//- ISAT HOM
	isat_HOM = new ISAT(scalingErrors_ISAT, epsilon_ISAT, thermodynamicsMapXML->NumberOfSpecies()+2);	

	// - Setting ISAT_HOM param
	isat_HOM->setMaxSizeBT(maxSizeBT);
	isat_HOM->setMaxSizeMRU(maxSizeMRU);
	isat_HOM->setMaxSizeMFU(maxSizeMFU);
	isat_HOM->setMaxSearchMRU(maxSearchMRU);
	isat_HOM->setMaxSearchMFU(maxSearchMFU);
	isat_HOM->setFlagSearchMRU(searchMRU);
	isat_HOM->setFlagSearchMFU(searchMFU);
	isat_HOM->setFlagClearingIfFull(clearIfFull);
	isat_HOM->setMaxGrowCoeff(maxGrowCoeff);
	isat_HOM->setMaxHeightCoeff(maxHeightCoeff);
	isat_HOM->setMaxTimeOldCoeff(maxTimeOldCoeff);
	isat_HOM->setMinUsedCoeff(minUsedCoeff);
	isat_HOM->setBalanceFactorRet(balanceFactorRetrieve);
	isat_HOM->setBalanceFactorAdd(balanceFactorAddition);
	isat_HOM->setQRType(qrSolver_ISAT);
	isat_HOM->setFlagCleanAndBalance(cleanAndBalance);

	{
		const unsigned int NC  = thermodynamicsMapXML->NumberOfSpecies();
	
		Info << endl << "ISAT parameters " << endl;
		Info << "   tolerance           : " << epsilon_ISAT << endl;
		Info << "   luFactorization     : " << luFactorization << endl; 	
		Info << "   qrFactorization     : " << qrFactorization << endl; 	
	
		Info << "   scalingFactors      : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingFactors_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingFactors_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(solverOptions.subDict("ISAT").subDict("scalingFactors").lookup("others")) << endl;
		Info << "       T               : " << scalingFactors_ISAT(NC) << endl;
		Info << "       tau             : " << scalingFactors_ISAT(NC+1) << endl;


		Info << "   clear if BT full    : " << clearIfFull << endl;

		Info << "   search in MRU       : " << searchMRU << endl;
		Info << "   search in MFU       : " << searchMFU << endl;
		Info << "   search BF           : " << searchBF << endl;

		Info << "   dimension parameters: " << endl;
		Info << "      max size BT      : " << maxSizeBT << endl;	
		if(searchMRU == true) 
		{
			Info << "      max size MRU     : " << maxSizeMRU << endl;
			Info << "      max search MRU   : " << maxSearchMRU << endl;
		}
		if(searchMFU == true) 
		{
			Info << "      max size MFU     : " << maxSizeMRU << endl;
			Info << "      max search MFU   : " << maxSearchMFU << endl;  
		}
	
		if (cleanAndBalance == true)
		{
			Info << "   balance parameters  : " << endl;
			Info << "      balanceFactorRetrieve   : " << balanceFactorRetrieve << endl;	
			Info << "      balanceFactorAddition   : " << balanceFactorAddition << endl;	
			Info << "      maxHeightCoeff          : " << maxHeightCoeff << endl;	
			Info << "      maxGrowCoeff            : " << maxGrowCoeff << endl;
			Info << "      minUsedCoeff            : " << minUsedCoeff << endl;
			Info << "      maxTimeOldCoeff         : " << maxTimeOldCoeff << endl; 
		}

		Info << "   scaling error       : " << endl;
		for(unsigned int i=0;i<NC;i++)  
		{
			if(scalingErrors_ISAT(i) != 1.) 
				Info << "      " << thermodynamicsMapXML->NamesOfSpecies()[i] << "\t       : " << scalingErrors_ISAT(i) << endl;
		}
		Info << "       others          : " << readScalar(solverOptions.subDict("ISAT").subDict("scalingErrors").lookup("others")) << endl;
		Info << "       T               : " << scalingErrors_ISAT(NC) << endl;
		Info << "       tau             : " << scalingErrors_ISAT(NC+1) << endl;
	
		Info << endl;
	}
}


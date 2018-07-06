for (unsigned int i=0;i<VCT.rows();i++)
{
	const double P_Pa = VCT(i,0);
	const double T = VCT(i,1);
	const double YF = VCT(i,4);
	const double YO = VCT(i,5);
	const double YI = VCT(i,6);
	const double YP1 = VCT(i,7);
	const double YP2 = VCT(i,8);
	const double YP3 = VCT(i,9);
	const double YP4 = VCT(i,10);
	const double YN2 = VCT(i,11);

	massFractions[1] = YF;
	massFractions[2] = YO;
	massFractions[3] = YI;
	massFractions[4] = YP1;
	massFractions[5] = YP2;
	massFractions[6] = YP3;
	massFractions[7] = YP4;
	massFractions[8] = YN2;


	const double CpMix = VCT(i,2)*4186;	// [J/kg/K]
	const double HMix = VCT(i,3)*4186;	// [J/kg/K]
	const double MWmix = VCT(i,20);

	const double OmegaF = VCT(i,12);
	const double OmegaO = VCT(i,13);
	const double OmegaI = VCT(i,14);
	const double OmegaP1 = VCT(i,15);
	const double OmegaP2 = VCT(i,16);
	const double OmegaP3 = VCT(i,17);
	const double OmegaP4 = VCT(i,18);
	const double OmegaN2 = VCT(i,19);

	Info << "T:    " << T << endl;
	Info << "FUEL: " << YF << endl;
	Info << "OX:   " << YO << endl;
	Info << "N2:   " << YN2 << endl;
	Info << "I:    " << YI << endl;
	Info << "P1:   " << YP1 << endl;
	Info << "P2:   " << YP2 << endl;
	Info << "P3:   " << YP3 << endl;
	Info << "P4:   " << YP4 << endl;
	Info << endl;

	const double cTot = P_Pa/PhysicalConstants::R_J_kmol/T;

	const double Calc_MWmix = virtualChemistryTable->MWMix(massFractions.GetHandle());
	const double Calc_CpMix = virtualChemistryTable->CpMix(T,P_Pa,massFractions.GetHandle());
	const double Calc_HMix = virtualChemistryTable->HMix(T,P_Pa,massFractions.GetHandle());
	virtualChemistryTable->FormationRates(cTot, Calc_MWmix, T, massFractions.GetHandle(), Omega.GetHandle());
	const double Calc_Qdot = virtualChemistryTable->Qdot(T,P_Pa,Omega.GetHandle());

	Info << MWmix << " " << Calc_MWmix << endl;
	Info << CpMix << " " << Calc_CpMix << endl;
	Info << HMix << " " << Calc_HMix << endl;
	Info << 0 << " " << Calc_Qdot << endl;
	
	Info << "OmegaF: " << OmegaF << " " << Omega[1] << endl;
	Info << "OmegaO: " << OmegaO << " " << Omega[2] << endl;
	Info << "OmegaI: " << OmegaI << " " << Omega[3] << endl;
	Info << "Omega1: " << OmegaP1 << " " << Omega[4] << endl;
	Info << "Omega2: " << OmegaP2 << " " << Omega[5] << endl;
	Info << "Omega3: " << OmegaP3 << " " << Omega[6] << endl;
	Info << "Omega4: " << OmegaP4 << " " << Omega[7] << endl;

	double sum = 0.;
	for(int i=0;i<thermodynamicsMapXML->NumberOfSpecies();i++)
		sum += Omega[i];
	Info << "Sum: " << sum << endl;

	getchar();
}

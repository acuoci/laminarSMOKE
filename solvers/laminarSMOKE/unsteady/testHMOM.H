if (hmom_test_ == true)
{
	double Y_H = 1e-5;
	double Y_OH = 1e-5; 
	double C_OH = 1e-5;
	double C_H = 1e-5;
	double C_H2O = 1e-5;
	double C_H2 = 1e-5;
	double C_C2H2 = 1e-5;
	double C_O2 = 1e-5;
	double C_PAH = 1e-5; 
	double temp=1000.;
	double press=1.;
	double vislam = 1.8e-5;
	double M00_norm = 1e-9;
	double M10_norm = 1e-9;
	double M01_norm = 1e-9;
	double N0_norm = 1e-9; 
	double M00_source = 0;
	double M10_source = 0;
	double M01_source = 0;
	double N0_source = 0;

	double in_fv;
	double in_nd;
	double in_dp;
	double in_dc;
	double in_np;
	double in_M00_source;
	double in_M10_source;
	double in_M01_source;
	double in_N0_source;

	std::ifstream fInput("test_solution.dat", std::ios::in);

	std::string dummy;
	for(int k=1;k<=25;k++)
		fInput >> dummy;

	for(int k=1;k<=195;k++)
	{
		Info << "Point " << k << endl;
		{
			fInput >> Y_H; 
			fInput >> Y_OH;
			fInput >> C_OH;
			fInput >> C_H; 
			fInput >> C_H2O;
			fInput >> C_H2;
			fInput >> C_C2H2;
			fInput >> C_O2;
			fInput >> C_PAH;
			fInput >> temp;
			fInput >> press;
			fInput >> vislam;
			fInput >> M00_norm;
			fInput >> M10_norm;
			fInput >> M01_norm;
			fInput >> N0_norm;
			fInput >> in_M00_source;
			fInput >> in_M10_source;
			fInput >> in_M01_source;
			fInput >> in_N0_source;
			fInput >> in_fv;
			fInput >> in_nd;
			fInput >> in_dp;
			fInput >> in_dc;
			fInput >> in_np;

			/*
			Info << "YH " << Y_H << endl; 
			Info << "YO " << Y_OH << endl;
			Info << "OH " << C_OH << endl;
			Info << "CH " << C_H << endl; 
			Info << "AC " << C_H2O << endl;
			Info << "H2 " << C_H2 << endl;
			Info << "C2 " << C_C2H2 << endl;
			Info << "O2 " << C_O2 << endl;
			Info << "PA " << C_PAH << endl;
			Info << "T  " << temp << endl;
			Info << "P  " << press << endl;
			Info << "mu " << vislam << endl;
			Info << "00 " << M00_norm << endl;
			Info << "10 " << M10_norm << endl;
			Info << "01 " << M01_norm << endl;
			Info << "N0 " << N0_norm << endl;
			Info << "00 " << in_M00_source << endl;
			Info << "10 " << in_M10_source << endl;
			Info << "01 " << in_M01_source << endl;
			Info << "N0 " << in_N0_source << endl;
			Info << "fv " << in_fv << endl;
			Info << "nd " << in_nd << endl;
			Info << "dp " << in_dp << endl;
			Info << "dc " << in_dc << endl;
			Info << "np " << in_np << endl;
			*/

		}

		Info << "Calculating source terms... " << endl;

		// Mass fractions [-]
		// Concentrations [mol/cm3]
		// Temperature [K]
		// Pressure [atm]
		// Moments (normalized) [mol/m3] ?
		// Moments sources (normalized) [mol/m3/s] ?
/*
		source_terms_soot_hmom_ (	&Y_H, &Y_OH, 
						&C_OH, &C_H, &C_H2O, &C_H2, &C_C2H2, &C_O2,
						&C_PAH, &temp, &press, &vislam,
						&M00_norm, &M10_norm, &M01_norm, &N0_norm, 
						&M00_source, &M10_source, &M01_source, &N0_source);

		Info << "M00_source: " << in_M00_source << "\t" << M00_source << "\t" << std::fabs(in_M00_source - M00_source)/std::max(std::fabs(M00_source), 1.e-32)*100. << endl; 
		Info << "M10_source: " << in_M10_source << "\t" << M10_source << "\t" << std::fabs(in_M10_source - M10_source)/std::max(std::fabs(M10_source), 1.e-32)*100. << endl; 
		Info << "M01_source: " << in_M01_source << "\t" << M01_source << "\t" << std::fabs(in_M01_source - M01_source)/std::max(std::fabs(M01_source), 1.e-32)*100. << endl; 
		Info << "N0_source:  " << in_N0_source << "\t"  << N0_source  << "\t" << std::fabs(in_N0_source - N0_source)/std::max(std::fabs(N0_source), 1.e-32)*100. << endl; 

		double one  = 1.;
		double zero = 0.;
		double minus_one = -1.;
		double minus_two = -2.;
		double three = 3.;

		double Avcoll = 1.-2./1.8;
		double Ascoll = 3./1.8-1.;
		double K_coll = 6./std::pow(36.*3.14159, 1./1.8);
		double K_sphe = std::pow(36.*3.14159, 1./3.);

		Info << "Calculating properties... " << endl;
		__mod_soot_hmom_MOD_calc_mom_var (&M00_norm, &M10_norm, &M01_norm, &N0_norm);
		const double mom0 = __mod_soot_hmom_MOD_soot_mom (&zero, &zero);

		const double fv = __mod_soot_hmom_MOD_soot_mom (&one, &zero);
		const double nd = mom0*1e-6;
		const double dp = 6.*__mod_soot_hmom_MOD_soot_mom (&one, &minus_one)/mom0 * 1.e9;
		const double dc = K_coll*__mod_soot_hmom_MOD_soot_mom(&Avcoll, &Ascoll)/mom0 * 1.e9;
		const double np = std::pow(K_sphe, -3.)*__mod_soot_hmom_MOD_soot_mom (&minus_two, &three)/mom0;

		Info << "fv:     " << in_fv << "\t" << fv << "\t" << std::fabs(in_fv - fv)/std::max(std::fabs(fv), 1.e-32)*100. << endl; 
		Info << "nd:     " << in_nd << "\t" << nd << "\t" << std::fabs(in_nd - nd)/std::max(std::fabs(nd), 1.e-32)*100. << endl;
		Info << "dp:     " << in_dp << "\t" << dp << "\t" << std::fabs(in_dp - dp)/std::max(std::fabs(dp), 1.e-32)*100. << endl;
		Info << "dc:     " << in_dc << "\t" << dc << "\t" << std::fabs(in_dc - dc)/std::max(std::fabs(dc), 1.e-32)*100. << endl;
		Info << "np:     " << in_np << "\t" << np << "\t" << std::fabs(in_np - np)/std::max(std::fabs(np), 1.e-32)*100. << endl;
*/
		getchar();
	}
}

#include <queue>

namespace OpenSMOKE
{
	VirtualChemistry::VirtualChemistry(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
						boost::filesystem::path file_name) :
	thermodynamicsMap_(thermodynamicsMap)
	{
		ns_ = thermodynamicsMap_.NumberOfSpecies();
		NSVector_.resize(ns_);
		Le_.resize(ns_);
		Le_.setConstant(1.);

		iOptimization_ = 170911;

		// Import table
		nv_ = 7;
		np_ = 16;

		// Open the table file
		std::ifstream fInput(file_name.c_str(), std::ios::in);

		// Read the first line
		std::string first_line;
		std::getline (fInput, first_line);

		// Memory allocation
		v_.resize(np_,nv_);
		x_.resize(np_);
		ratios_.resize(np_,nv_);
		table_.resize(nv_);

		// Fill the table
		for (unsigned int i=0;i<np_;i++)
		{
			fInput >> x_(i); 
			for (unsigned int k=0;k<nv_;k++)
				fInput >> v_(i,k);
		}

		// Check closure
		std::string dummy;
		fInput >> dummy;
		if (dummy != "END")
		{
			Info << "Virtual Chemistry: error in reading the table. Expected END, Found " << dummy << endl;
			abort();
		}

		// Close input file
		fInput.close();

		// Min and max values
		min_x_ = x_.minCoeff();
		max_x_ = x_.maxCoeff();

		Info << min_x_ << " " << max_x_ << endl; 

		// Precalculation of ratios to be used in interpolation
		for (unsigned int i=1;i<np_-1;i++)
			for (unsigned int k=0;k<nv_;k++)
				ratios_(i-1,k) = (v_(i,k)-v_(i-1,k))/(x_(i)-x_(i-1));

		Info << "Precalculation" << endl;

		// Indices of products
		I_index_  = thermodynamicsMap_.IndexOfSpecies("I")-1;
		P1_index_ = thermodynamicsMap_.IndexOfSpecies("P1")-1;
		P2_index_ = thermodynamicsMap_.IndexOfSpecies("P2")-1;
		P3_index_ = thermodynamicsMap_.IndexOfSpecies("P3")-1;
		P4_index_ = thermodynamicsMap_.IndexOfSpecies("P4")-1; 
	}

	void VirtualChemistry::SetOptimization(const unsigned int value)
	{
		if (value != 170911 && value != 171013)
		{
			Info << "Virtual Chemistry: error in choosing the optimization type: " << value << endl;
			Info << "                   available: 170911 | 171013" << endl;
			abort();
		}

		iOptimization_ = value;
	}

	void VirtualChemistry::SetFuel(const std::string name, const double mw)
	{
		fuel_name_  = name;
		fuel_mw_    = mw;
		fuel_index_ = thermodynamicsMap_.IndexOfSpecies(fuel_name_)-1; 
	}

	void VirtualChemistry::SetOxidizer(const std::string name, const double mw)
	{
		oxidizer_name_  = name;
		oxidizer_mw_    = mw;
		oxidizer_index_ = thermodynamicsMap_.IndexOfSpecies(oxidizer_name_)-1; 
	}

	void VirtualChemistry::SetInert(const std::string name, const double mw)
	{
		inert_name_  = name;
		inert_mw_    = mw;
		inert_index_ = thermodynamicsMap_.IndexOfSpecies(inert_name_)-1; 
	}

	void VirtualChemistry::SetTransportProperties(const double mu0, const double T0, const double Beta0, const double Pr0)
	{
		mu0_   = mu0;
		T0_    = T0;
		Beta0_ = Beta0;
		Pr0_   = Pr0;
	}

	void VirtualChemistry::SetLewisNumber(const unsigned int i, const double Le)
	{
		Le_(i) = Le;
	}

	double VirtualChemistry::MWMix(double* Y)
	{
		double mw = 0.;
		if (iOptimization_ == 170911)
		{
			const double mw_table = InterpolationSingleValue(Y[inert_index_], 0);
			const double sum = 1.-(Y[fuel_index_]+Y[oxidizer_index_]+Y[inert_index_]);
			mw = 1./(Y[fuel_index_]/fuel_mw_+Y[oxidizer_index_]/oxidizer_mw_+Y[inert_index_]/inert_mw_+sum/mw_table);
		}
		else if (iOptimization_ == 171013)
		{
			const double mw_table_ = InterpolationSingleValue(Y[inert_index_], 0);
			mw = 1./( Y[fuel_index_]/fuel_mw_ + Y[oxidizer_index_]/oxidizer_mw_ + Y[inert_index_]/inert_mw_ + Y[I_index_]/mw_table_ +
				  Y[P1_index_]/15.336 + Y[P2_index_]/16.345 + Y[P3_index_]/93.850 + Y[P4_index_]/37.960 );
		}
		return mw;     
	}

	double VirtualChemistry::CpMix(const double T, const double P_Pa, double* Y)
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the specific heats of species in mass units [J/kg/K]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.cpMolar_Species(NSVector_.data());

		// Return the specific heat of the mixture [J/kg]
		double sum = 0;
		for (unsigned int i=0;i<ns_;i++)
			sum += Y[i]*NSVector_(i);
		return sum;
	}

	double VirtualChemistry::HMix(const double T, const double P_Pa, double* Y)
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the enthalpies of species in mass units [J/kg]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.hMolar_Species(NSVector_.data());

		// Return the enthalpy of the mixture [J/kg]
		double sum = 0;
		for (unsigned int i=0;i<ns_;i++)
			sum += Y[i]*NSVector_(i);
		return sum;
	}

	double VirtualChemistry::Qdot(const double T, const double P_Pa, double* Omega) 
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the enthalpies of species in mass units [J/kg]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.hMolar_Species(NSVector_.data());

		// Returns the reaction heat [J/m3/s]
		double sum = 0;
		for (unsigned int i=0;i<ns_;i++)
			sum -= Omega[i]*NSVector_(i);
		return sum;
	}

	double VirtualChemistry::DynamicViscosity(const double T)
	{
		return mu0_*std::pow(T/T0_, Beta0_);
	}

	double VirtualChemistry::ThermalConductivity(const double mu, const double cp)
	{
		return mu*cp/Pr0_;
	}

	void VirtualChemistry::MassDiffusionCoefficients(const double lambda, const double rho, const double cp, double* Dmix)
	{
		const double alpha = lambda/rho/cp;
		for (unsigned int i=0;i<ns_;i++)
			Dmix[i] = alpha/Le_(i);
	}

	void VirtualChemistry::FormationRates(const double cTot, const double MW, const double T, double* Y, double* Omega)
	{
		// Mass fraction of N2
		const double YN2 = Y[inert_index_];
		Interpolation(YN2, table_.data());

		// Mass fractions
		const double YF = Y[fuel_index_];
		const double YO = Y[oxidizer_index_];
		const double YI = Y[I_index_];

		// Concentrations [mol/cm3]
		const double CF = cTot * YF*MW/1. / 1000.;
		const double CO = cTot * YO*MW/1. / 1000.;
		const double CI = cTot * YI*MW/1. / 1000.;

		double r1 = 0.;
		double r2 = 0.;

		if (iOptimization_ == 170911)
		{
			// Kinetic constants [mol, cm3, s]
			const double k1 = 1.5384796696859927E+18*std::pow(T,0.)*std::exp(-3.5075174514637045E+04/1.987/T) * table_(1);	
			const double k2 = 3.9225000838284247E+18*std::pow(T,0.)*std::exp(-8.5680603523860715E+04/1.987/T);

			// Reaction rates [kmol/m3/s]
			r1 = k1*std::pow(CF, 1.7099829697856470)*std::pow(CO, 0.8686294702424591)*1000;
			r2 = k2*std::pow(CI, 2.3399221348980621*table_(2))*1000;

			// 0.2FUEL + 0.8OX => I 		1.5384796696859927E+18	0.	3.5075174514637045E+04
			// FORD /FUEL 1.7099829697856470/
			// FORD /OX   0.8686294702424591/

			// I => 0.25P1+0.25P2+0.25P3+0.25P4 	3.9225000838284247E+18	0.	8.5680603523860715E+04
			// FORD /I 2.3399221348980621/
		}
		else if (iOptimization_ == 171013)
		{
			// Kinetic constants [mol, cm3, s]
			const double k1 = 4.99463851138e+18*std::pow(T,0.)*std::exp(-50292.9081962/1.987/T) * table_(1);	
			const double k2 = 5.62191853964e+18*std::pow(T,0.)*std::exp(-111552.26047/1.987/T);

			// Reaction rates [kmol/m3/s]
			r1 = k1*std::pow(CF, 0.000241938688303)*std::pow(CO, 2.49626471112)*1000;
			r2 = k2*std::pow(CI, 1.83822628643*table_(2))*1000;

			// 0.2FUEL + 0.80OX => I 	 		4.99463851138e+18    0.0    50292.9081962
			// FORD /FUEL 0.000241938688303/
			// FORD /OX 2.49626471112/

			// I => 0.25P1+0.25P2+0.25P3+0.25P4 	 5.62191853964e+18   0.0    111552.26047
			// FORD /I 1.83822628643/
		}
		
		// Formation rates [kg/m3/s]
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the formation rates of species in mass units [J/kg]
		Omega[fuel_index_] = -0.2*r1;		// FUEL
		Omega[oxidizer_index_] = -0.8*r1;	// OX
		Omega[2] = r1-r2;			// I
		Omega[3] = table_(3)*r2;		// P1
		Omega[4] = table_(4)*r2;		// P2
		Omega[5] = table_(5)*r2;		// P3
		Omega[6] = table_(6)*r2;		// P4
		Omega[7] = 0.;				// N2
	}

	void VirtualChemistry::Interpolation(const double x, double* v)
	{
		if (x<=min_x_)
		{
			for (unsigned int k=0;k<nv_;k++)
				v[k] = v_(0,k);
			return;
		}
		else if (x>=max_x_)
		{
			for (unsigned int k=0;k<nv_;k++)
				v[k] = v_(np_-1,k);
			return;
		}		
		else
		{
			for (unsigned int i=1;i<np_-1;i++)
				if (x_(i) >= x)
				{
					for (unsigned int k=0;k<nv_;k++)
						v[k] = v_(i-1,k) +  ratios_(i-1,k)*(x-x_(i-1));
					return;
				}
		}
	}

	double VirtualChemistry::InterpolationSingleValue(const double x, const unsigned int k)
	{
		if (x<=min_x_)
		{
			return v_(0,k);
		}
		else if (x>=max_x_)
		{
			return v_(np_-1,k);
		}		
		else
		{
			for (unsigned int i=1;i<np_-1;i++)
				if (x_(i) >= x)
				{
					return (v_(i-1,k) +  ratios_(i-1,k)*(x-x_(i-1)));
				}
		}
	}
}

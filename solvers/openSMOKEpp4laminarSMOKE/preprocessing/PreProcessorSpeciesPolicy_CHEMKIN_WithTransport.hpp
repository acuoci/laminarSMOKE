/*----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Author: Alberto Cuoci <alberto.cuoci@polimi.it>                       |
|   CRECK Modeling Group <http://creckmodeling.chem.polimi.it>            |
|   Department of Chemistry, Materials and Chemical Engineering           |
|   Politecnico di Milano                                                 |
|   P.zza Leonardo da Vinci 32, 20133 Milano                              |
|                                                                         |
|-------------------------------------------------------------------------|
|                                                                         |
|   This file is part of OpenSMOKE++ framework.                           |
|                                                                         |
|	License                                                               |
|                                                                         |
|   Copyright(C) 2014, 2013, 2012  Alberto Cuoci                          |
|   OpenSMOKE++ is free software: you can redistribute it and/or modify   |
|   it under the terms of the GNU General Public License as published by  |
|   the Free Software Foundation, either version 3 of the License, or     |
|   (at your option) any later version.                                   |
|                                                                         |
|   OpenSMOKE++ is distributed in the hope that it will be useful,        |
|   but WITHOUT ANY WARRANTY; without even the implied warranty of        |
|   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         |
|   GNU General Public License for more details.                          |
|                                                                         |
|   You should have received a copy of the GNU General Public License     |
|   along with OpenSMOKE++. If not, see <http://www.gnu.org/licenses/>.   |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "math/PhysicalConstants.h"
#include "CollisionIntegralMatrices.hpp"

namespace OpenSMOKE
{
	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::BOLTZMANN		= PhysicalConstants::kBoltzmann*1.e7;	// [erg/K]

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::BOLTZMANN3		= boost::math::pow<3>(BOLTZMANN);		// [erg3/K3]

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::MU_MIN			= 1.e-38;								// [sqrt(A3.erg)]

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ZROTA			= 0.5*std::pow(boost::math::constants::pi<double>(), 1.5);	

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ZROTB			= 2.0+0.25*boost::math::pow<2>(boost::math::constants::pi<double>());

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ZROTC			= std::pow(boost::math::constants::pi<double>(), 1.5);

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::CONST_2SUPI	= 2./boost::math::constants::pi<double>();

	template<typename Species>
	const double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::CONST_5SU3R	= 5./(3.*PhysicalConstants::R_J_mol);

	template<typename Species>
	PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::PreProcessorSpeciesPolicy_CHEMKIN_WithTransport() 
	{
	}

	template<typename Species>
	PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::PreProcessorSpeciesPolicy_CHEMKIN_WithTransport(const PreProcessorSpeciesPolicy_CHEMKIN_WithTransport& orig) {
	}

	template<typename Species>
	PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::~PreProcessorSpeciesPolicy_CHEMKIN_WithTransport() 
	{
		this->species_.clear();
		this->names_.clear();
		delete [] fittingBinaryDiffusivities;
	}


	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::Allocate()
	{
		this->NC = boost::lexical_cast<unsigned int>(this->species_.size());

		ChangeDimensions(this->NC,&shape_factor,true);
		ChangeDimensions(this->NC,&epsylon_over_kb,true);
		ChangeDimensions(this->NC,&sigma,true);
		ChangeDimensions(this->NC,&mu,true);
		ChangeDimensions(this->NC,&alfa,true);
		ChangeDimensions(this->NC,&zRot298,true);
		ChangeDimensions(this->NC,&kb_over_epsylon,true);
		ChangeDimensions(this->NC,&muStar,true);
		ChangeDimensions(this->NC,&epsylon,true);

		ChangeDimensions(this->NC,&this->MW,true);
		ChangeDimensions(this->NC,&this->uMW,true);

		// Mass Diffusivity
		ChangeDimensions(this->NC,this->NC,&deltajkStar,true);
		ChangeDimensions(this->NC,this->NC,&coeff_Djk,true);
		ChangeDimensions(this->NC,this->NC,&Djk,true);

		// Viscosity
		ChangeDimensions(this->NC,&deltakStar,true);
		ChangeDimensions(this->NC,&omega22k,true);
		ChangeDimensions(this->NC,&TkStar,true);
		ChangeDimensions(this->NC,&eta,true);
		ChangeDimensions(this->NC,&coeff_eta,true);

		// Thermal conductivity
		ChangeDimensions(this->NC,&fVib,true);
		ChangeDimensions(this->NC,&fRot,true);
		ChangeDimensions(this->NC,&fTrans,true);
		ChangeDimensions(this->NC,&cVtrans,true);
		ChangeDimensions(this->NC,&cVrot,true);
		ChangeDimensions(this->NC,&cVvib,true);
		ChangeDimensions(this->NC,&f298,true);
		ChangeDimensions(this->NC,&fT,true);
		ChangeDimensions(this->NC,&zRot,true);
		ChangeDimensions(this->NC,&rho,true);
		ChangeDimensions(this->NC,&Dkk,true);
		ChangeDimensions(this->NC,&coeff_Dkk,true);
		ChangeDimensions(this->NC,&omega11kk,true);
		ChangeDimensions(this->NC,&lambda,true);

		// Thermodynamic properties
		ChangeDimensions(this->NC,&this->Cp,true);
		ChangeDimensions(this->NC,&this->Cv,true);

		// Transport data
		for(unsigned int i=1;i<=this->NC;i++)
		{
			shape_factor[i]=this->species_[i-1].shape_factor();		// 0=monoatomic 1=linear 2=c.non linear
			epsylon_over_kb[i]=this->species_[i-1].epsylon_over_kb();	// [K]
			sigma[i]=this->species_[i-1].sigma();						// [A]
			mu[i]=this->species_[i-1].mu() * 1.e-6;					// [Debye --> sqrt(A3.erg)]
			alfa[i]=this->species_[i-1].alfa();						// [A3]
			zRot298[i]=this->species_[i-1].zRot298();					// [-]
		}

		// Molecular wights [kg/kmol]
		for(unsigned int i=1;i<=this->NC;i++)
		{
			this->MW[i] = this->species_[i-1].MolecularWeight();
			this->uMW[i] = 1./this->MW[i];
		}

		// Lennard-Jones Potentials [erg]
		Product(BOLTZMANN, epsylon_over_kb, &epsylon);
		for(unsigned int i=1;i<=this->NC;i++)
			kb_over_epsylon[i] = 1./epsylon_over_kb[i];
	
		// Dipolar moments (check)
		for(unsigned int i=1;i<=this->NC;i++)
			if(mu[i] < MU_MIN)  mu[i] = 0.;  // dipolar moment [sqrt(A3.erg)]
	
		// Reduced dipolar moment for polar molecules
		for(unsigned int i=1;i<=this->NC;i++)
			muStar[i] = mu[i] / sqrt(epsylon[i]*boost::math::pow<3>(sigma[i])); // [-]
	
		// The thermal diffusivities are calculated only the light species
		for(unsigned int i=1;i<=this->NC;i++)
			if (this->MW[i] <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)	iThermalDiffusionRatios.Append(i);
		ChangeDimensions(iThermalDiffusionRatios.Size(), this->NC, &TetaStar, true);
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::Setup()
	{
		Allocate();
		InitializeMassDiffusivity();
		InitializeViscosity();
		InitializeThermalConductivity();

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteFittingCoefficientsOnASCIIFile(const std::string file_name) const
	{
		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);

		WriteViscosityFittingCoefficients(fOutput);
		WriteThermalConductivityFittingCoefficients(fOutput);
		WriteBinaryDiffusivityFittingCoefficients(fOutput);
		WriteThermalDiffusionRatiosFittingCoefficients(fOutput);

		fOutput.close();

		return true;
	}

	// Reduced molecular weight [kg/kmol]
	double Mjk(const double MWj, const double MWk)
	{
		return 2.*MWj*MWk/(MWj+MWk);
	}

	// Dipolar collision moment [sqrt(A3.erg)]
	double mujk(const double mu_j, const double mu_k)
	{
		return sqrt(mu_j * mu_k);	// 
	}

	// Dipole effective collision moment (5.14)
	double mujkStar(const double mu_j, const double mu_k, const double epsjk_over_kb, const double sigmajk)
	{
		return mujk(mu_j, mu_k)  / sqrt( (epsjk_over_kb*PhysicalConstants::kBoltzmann*1.e7) * boost::math::pow<3>(sigmajk)); // [-]
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::LennardJonesPotentialAndCollisionDiameter(const int j, const int k, double& epsjk_over_kb, double& sigmajk)
	{
		// Lennard-Jones potentials
		epsjk_over_kb = sqrt(epsylon_over_kb[j] * epsylon_over_kb[k]);	// [K]

		// Lennard-Jones collision diameters
		sigmajk = .5 * (sigma[j] + sigma[k]);	

		// Case 2: polar and non-polar species
		if(mu[j] == 0. && mu[k] != 0.)
		{
			// a. non-polar species: reduced polarizability (5.13)
			//    alfa in [A3] - sigma in [A]
			double alfaStar = alfa[j] / boost::math::pow<3>(sigma[j]);	//[-]

			// c. Correction coefficient (5.12)
			double xsi = 1. + .25 * alfaStar * boost::math::pow<2>(muStar[k]) * sqrt(epsylon[k] / epsylon[j]); //[-]

			// d. corrections
			epsjk_over_kb *= (xsi * xsi);		// (5.9)
			sigmajk *= std::pow(xsi, -1./6.);			// (5.10)
		}

		else if(mu[k] == 0. && mu[j] != 0.)
		{
			// a. polar species: reduced polarizability (5.13)
			double alfaStar = alfa[k] / boost::math::pow<3>(sigma[k]);	//[-]

			// c. Correction coefficient (5.12)
			double xsi = 1. + .25 * alfaStar * boost::math::pow<2>(muStar[j]) * sqrt(epsylon[j] / epsylon[k]); //[-]

			// d. corrections
			epsjk_over_kb *= (xsi * xsi);		// (5.9)
			sigmajk *= std::pow(xsi,-1./6.);			// (5.10)
		}
	}

	template<typename Species>
	double PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::UncorrectedLennardJonesPotential(const int j, const int k) const
	{
		return sqrt(epsylon_over_kb[j] * epsylon_over_kb[k]);	// [K]
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::InitializeMassDiffusivity()
	{
		for(unsigned int j=1;j<=this->NC;j++)
			for(unsigned int k=1;k<=this->NC;k++)
			{
				double epsylonjk_over_kb_, sigmajk_;
				LennardJonesPotentialAndCollisionDiameter(j, k, epsylonjk_over_kb_, sigmajk_);

				deltajkStar[j][k] = 0.5*boost::math::pow<2>( mujkStar( mu[j], mu[k], epsylonjk_over_kb_, sigmajk_) );
				coeff_Djk[j][k] = 2.6693e-7  / ( sqrt(Mjk(this->MW[j],this->MW[k]))*boost::math::pow<2>(sigmajk_) );
			}
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::InitializeThermalConductivity()
	{
		for(unsigned int k=1;k<=this->NC;k++)
		{
			// Parker-Brau-Jonkman at 298K (5.33)
			double aux = sqrt(epsylon_over_kb[k] / 298.); //[-]
			f298[k] = 1. + aux * ( ZROTA + aux * (ZROTB + aux * ZROTC)); // [-]
		}

		// Translational and rotational parts of specific heat at constant volume [J/mol/K]
		cVtrans = 1.5 * PhysicalConstants::R_J_mol; 
		for(unsigned int k=1;k<=this->NC;k++)
		{
			if(shape_factor[k] == 0)		// monoatomic species
				cVrot[k] = 0.;				

			else if(shape_factor[k] == 1)				// linear chain molecule
				cVrot[k] = PhysicalConstants::R_J_mol;	

			else							// non-linear chain molecule
				cVrot[k] = cVtrans[k];		
		}

		// Coefficient for the evaluation of binary coefficients
		for(unsigned int k=1;k<=this->NC;k++)
				coeff_Dkk[k] = 2.6693e-7  / ( sqrt(this->MW[k])*boost::math::pow<2>(sigma[k]) );
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::InitializeViscosity()
	{
		// Coefficient for the evaluation of viscosity of single species
		for(unsigned int k=1;k<=this->NC;k++)
			coeff_eta[k] = 26.693e-7 * sqrt(this->MW[k]) / boost::math::pow<2>(sigma[k]);
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SpeciesViscosities(const double T)
	{
		// Evaluates the viscosities of species [kg/m/s]
		// The temperature must be provided in K

		ReducedTemperatureSingleSpecies(T);
		Omega22k(TkStar, deltakStar, omega22k, this->names_);
		SingleSpeciesViscosity(T);
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SpeciesThermalConductivities(const double T)
	{
		// Evaluates the thermal conductivities of species [W/m/K]
		// The temperature must be provided in K

		ComputeFT(T);
		ComputeZrotational();
		ComputeCvVibrational();
		DensitySingleSpecies(T);
		Omega11kk();
		SelfDiffusionCoefficient(T);
		ThermalConductivitySingleSpecies();
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SpeciesBinaryDiffusivities(const double T, const double P_bar)
	{
		// Evaluates the binary diffusion coefficients [m2/s]
		// The temperature must be provided in K and the pressure in bar

//		ReducedTemperatures(T);
//		Omega11jk(T);
		BinaryDiffusionCoefficients(T,P_bar);
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SpeciesThermalDiffusionRatios(const double T)
	{
	
		for(int i=1;i<=iThermalDiffusionRatios.Size();i++)
			for(unsigned int j=1;j<=this->NC;j++)
			{
				int k = iThermalDiffusionRatios[i];

				double epsylonjk_over_kb_, sigmajk_;
				LennardJonesPotentialAndCollisionDiameter(k, j, epsylonjk_over_kb_, sigmajk_);
				double TSLOG = log(T/epsylonjk_over_kb_);

				double T1 = TSLOG;
				double T2 = TSLOG*T1;
				double T3 = TSLOG*T2;
				double T4 = TSLOG*T3;
				double T5 = TSLOG*T4;
				double T6 = TSLOG*T5;

				// Evaluates A*, B* and C* (ratios of collision integrals) (5.55, 5.56, 5.57)
				double ASTAR =	CollisionIntegralMatrices::FITASTAR_1 + CollisionIntegralMatrices::FITASTAR_2*T1 + 
								CollisionIntegralMatrices::FITASTAR_3*T2 + CollisionIntegralMatrices::FITASTAR_4*T3 + 
								CollisionIntegralMatrices::FITASTAR_5*T4 + CollisionIntegralMatrices::FITASTAR_6*T5 + 
								CollisionIntegralMatrices::FITASTAR_7*T6;
				double BSTAR =	CollisionIntegralMatrices::FITBSTAR_1 + CollisionIntegralMatrices::FITBSTAR_2*T1 + 
								CollisionIntegralMatrices::FITBSTAR_3*T2 + CollisionIntegralMatrices::FITBSTAR_4*T3 + 
								CollisionIntegralMatrices::FITBSTAR_5*T4 + CollisionIntegralMatrices::FITBSTAR_6*T5 + 
								CollisionIntegralMatrices::FITBSTAR_7*T6;
				double CSTAR =	CollisionIntegralMatrices::FITCSTAR_1 + CollisionIntegralMatrices::FITCSTAR_2*T1 + 
								CollisionIntegralMatrices::FITCSTAR_3*T2 + CollisionIntegralMatrices::FITCSTAR_4*T3 + 
								CollisionIntegralMatrices::FITCSTAR_5*T4 + CollisionIntegralMatrices::FITCSTAR_6*T5 + 
								CollisionIntegralMatrices::FITCSTAR_7*T6;

				// Evaluates the binary thermal diffusion ratios (5.54)
				TetaStar[i][j] = 7.5*(2.*ASTAR+5.)*(6.*CSTAR-5.)/(ASTAR*(16.*ASTAR-12.*BSTAR+55.)) *
								 (this->MW[k]-this->MW[j])/(this->MW[k]+this->MW[j]);
			}
	}


	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::BinaryDiffusionCoefficients(const double T, const double P_bar)
	{
		// The matrix is symmetric and only theupper half is calculated
		// The main diagonal is not calculated here (it corresponds to the self diffusion
		// coefficients), becuase it is calculated by a different function

		double T_P = std::pow(T,1.5) / P_bar;

		for(unsigned int j=1;j<=this->NC;j++)
			for(unsigned int k =j+1;k<=this->NC;k++)
			{
				Djk[k][j] = Djk[j][k] = coeff_Djk[j][k] * T_P / CollisionIntegral11(T/UncorrectedLennardJonesPotential(j,k),deltajkStar[j][k]);

				/*
				// Comment (Alberto Cuoci): in my opinion this version is more correct
				// However the results obtained with the version reported above are much closer to what is calculated by CHEMKIN
				// I have no explanation about this
				
				double epsylonjk_over_kb_, sigmajk_;
				LennardJonesPotentialAndCollisionDiameter(j, k, epsylonjk_over_kb_, sigmajk_);
				Djk[k][j] = Djk[j][k] = coeff_Djk[j][k] * T_P / CollisionIntegral11(T/epsylonjk_over_kb_,deltajkStar[j][k]);
				*/
			}

		for (unsigned k = 1; k <= this->NC; k++)
			Djk[k][k] = coeff_Dkk[k] * T_P / CollisionIntegral11(T*kb_over_epsylon[k], deltakStar[k]);
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ComputeFT(const double T)
	{
		// Parker-Brau-Jonkman function (5.34)
		// It is used for calculating the thermal conductivites of single species

		double uT = 1./T;
		for(unsigned int k=1;k<=this->NC;k++)
		{
			double aux = sqrt(epsylon_over_kb[k] * uT);
			fT[k] = 1. + aux * ( ZROTA + aux * ( ZROTB + aux * ZROTC));
		}
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ComputeZrotational(void)
	{
		// -Rotational relaxation collision number (5.33)
		// It is used for calculating the thermal conductivites of single species

		for(unsigned int k=1;k<=this->NC;k++)
			zRot[k] = zRot298[k] * f298[k] / fT[k];
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ComputeCvVibrational(void)
	{
		// Vibrational part of Cp at constant volume [J/kmol/K]
		// It is used for calculating the thermal conductivites of single species

		for(unsigned int k=1;k<=this->NC;k++)
			if(shape_factor[k] != 0)
				cVvib[k] = this->Cv[k]*this->MW[k]*1.e-3 - cVtrans[k] - cVrot[k];
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::DensitySingleSpecies(const double T)
	{
		// Evaluation of densities of single species [kg/m3]
		// It is used for calculating the thermal conductivites of single speciesnot multiplied
		// by the pressure

		double coeff= 100. / (PhysicalConstants::R_J_mol * T);
		for(unsigned int k=1;k<=this->NC;k++)
			rho[k] = coeff*this->MW[k];
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ThermalConductivitySingleSpecies(void)
	{
		// Evaluation of thermal conductivities of single species [W/m/K]
		// The viscosities eta must be calculate before calling this function
		// The following data are necessary: fT, zRot, rho, CvVib e Dkk
		for(unsigned int k=1;k<=this->NC;k++)
		{
			fVib[k] = rho[k] * Dkk[k] / eta[k];		// (5.20)

			double A = 2.5-fVib[k];												// (5.21)
			double B = zRot[k] + CONST_2SUPI*(CONST_5SU3R*cVrot[k] + fVib[k]);	// (5.22)
			fTrans[k] = 2.50*(1.-CONST_2SUPI*cVrot[k]/cVtrans[k]*A/B);			// (5.18)

			if(shape_factor[k] == 0)
				lambda[k] = eta[k] * this->uMW[k] * (fTrans[k]*cVtrans[k]) ;			// (5.30)

			else
			{
				fRot[k] = fVib[k]*(1.+CONST_2SUPI*A/B);										// (5.19)
				double rdf=(fTrans[k]*cVtrans[k] + fRot[k]*cVrot[k] + fVib[k]*cVvib[k]);	// (5.17)

				lambda[k] = eta[k] * this->uMW[k] * (fTrans[k]*cVtrans[k] + fRot[k]*cVrot[k] + fVib[k]*cVvib[k]);	// (5.17)
			}
		}

		lambda *= 1000.;		// [W/m/K]
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::Omega11kk(void)
	{
		for(unsigned int k=1;k<=this->NC;k++)
			omega11kk[k] = CollisionIntegral11(TkStar[k],deltakStar[k]);
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SelfDiffusionCoefficient(const double T)
	{
		// The main diagonal of this matrix is needed only for evaluating the 
		// thermal conductivities of single species and it is not necessary for
		// the mass diffusivities
		// The coefficients of the matrix are not divided by the pressure in bar

		double T_P = std::pow(T,1.5);
		for (unsigned k = 1; k <= this->NC; k++)
			Dkk[k] = coeff_Dkk[k] * T_P / omega11kk[k];		// (5.31)
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::ReducedTemperatureSingleSpecies(const double T)
	{
		// This vector is used only as input coordinate for the evaluation of the
		// collision integral Omega22

		for(unsigned int k=1;k<=this->NC;k++)
			TkStar[k] = T*kb_over_epsylon[k];									// (5.2)
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SingleSpeciesViscosity(const double T)
	{
		// Viscosity of single species [kg/m/s]
		// Standard kinetic theory expression (5.1)

		double sqrT = sqrt(T);
		for(unsigned int k=1;k<=this->NC;k++)
			eta[k] = coeff_eta[k] * sqrT / omega22k[k];
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::Fitting()
	{
		double tStart = OpenSMOKE::OpenSMOKEGetCpuTime();

		const int nPoints = 10;
		const double TMIN =  300.;
		const double TMAX = 3600.;
		const double dT=(TMAX-TMIN)/double(nPoints-1);

		Eigen::Matrix4d XTX;
		Eigen::MatrixXd XT(4, nPoints);

		// Assembling X and XT Matrices
		{
			Eigen::MatrixXd X(nPoints, 4);
		
			double T=TMIN;
			for (int i=0;i<nPoints;i++)
			{
				double logT = log(T);

				X(i,0)=1.;
				X(i,1)=logT;
				X(i,2)=logT*logT;
				X(i,3)=logT*logT*logT;
				T+=dT;
			}

			XT = X.transpose();
			XTX=XT*X;
		}
	
		// Fitting Viscosity
		{
			std::cout << " * Fitting viscosity..." << std::endl;

			Eigen::MatrixXd y(nPoints, this->NC);
			Eigen::MatrixXd Y(4, this->NC);

			double T=TMIN;
			for (int i=0;i<nPoints;i++)
			{
				SpeciesViscosities(T);
				for (unsigned int j=1;j<=this->NC;j++)
					y(i,j-1) = log(eta[j]);
				T+=dT;
			}
	
			Y=XT*y;
			fittingEta = XTX.fullPivLu().solve(Y);
		}

		// Fitting Conductivities
		{
			std::cout << " * Fitting thermal conductivity..." << std::endl;

			Eigen::MatrixXd y(nPoints, this->NC);
			Eigen::MatrixXd Y(4, this->NC);

			double T=TMIN;
			for (int i=0;i<nPoints;i++)
			{
				this->SpeciesCp(T);
				this->SpeciesCv();
				SpeciesViscosities(T);
				SpeciesThermalConductivities(T);
				for (unsigned int j=1;j<=this->NC;j++)
					y(i,j-1) = log(lambda[j]);
				T+=dT;
			}
	
			Y=XT*y;
			fittingLambda.resize(4, this->NC);
			fittingLambda = XTX.fullPivLu().solve(Y);
		}

		// Fitting Mass Diffusivities
		{
			std::cout << " * Fitting mass diffusivity..." << std::endl;

			fittingBinaryDiffusivities = new Eigen::MatrixXd[this->NC];
			Eigen::MatrixXd* y = new Eigen::MatrixXd[this->NC];
			for (unsigned int j=1;j<=this->NC;j++)
				y[j-1].resize(nPoints, this->NC);

			double T=TMIN;
			for (int i=0;i<nPoints;i++)
			{
				SpeciesBinaryDiffusivities(T,1.);

				for (unsigned int j=1;j<=this->NC;j++)
					for (unsigned int k=1;k<=this->NC;k++)
						y[j-1](i,k-1) = log(Djk[j][k]);
				T+=dT;
			}

			for (unsigned int j=1;j<=this->NC;j++)
			{
				Eigen::MatrixXd Y(4, this->NC);
	
				Y=XT*y[j-1];
				y[j-1].resize(0,0);

				fittingBinaryDiffusivities[j-1].resize(4, this->NC);
				fittingBinaryDiffusivities[j-1] = XTX.fullPivLu().solve(Y);
			}

			delete [] y;
		}

		// Thermal Diffusion Ratios
		{
			std::cout << " * Fitting thermal diffusion ratios..." << std::endl;

			const unsigned int NT = iThermalDiffusionRatios.Size();
			if (NT>0)
			{
				Eigen::MatrixXd y(nPoints, this->NC*NT);
				Eigen::MatrixXd Y(4, this->NC*NT);
				Eigen::MatrixXd X(nPoints, 4);

				double T=TMIN;
				for (int i=0;i<nPoints;i++)
				{
					SpeciesThermalDiffusionRatios(T);

					for (unsigned int k=1;k<=NT;k++)
						for (unsigned int j=1;j<=this->NC;j++)
							y(i,j+this->NC*(k-1)-1) = TetaStar[k][j];

					X(i,0)=1.;
					X(i,1)=T;
					X(i,2)=T*T;
					X(i,3)=T*T*T;

					T+=dT;
				}

				XT=X.transpose();
				XTX=XT*X;

				Y=XT*y;
				fittingTetaBinary.resize(4, this->NC*NT);
				fittingTetaBinary = XTX.fullPivLu().solve(Y);
			}
		}

		double tEnd = OpenSMOKE::OpenSMOKEGetCpuTime();

		std::cout << " * Transport properties fitted in: " << tEnd-tStart << " s" << std::endl;

		return true;
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteViscosityFittingCoefficients(std::ostream &fOutput) const
	{
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                       VISCOSITY FITTING COEFFICIENTS    " << std::endl;
		fOutput << std::endl;
		fOutput << "             mu = exp( A + B*logT + C*(logT)^2 + D*(logT)^3 )  [kg/m/s]      " << std::endl;
		fOutput << std::endl;
		fOutput << "    Species                         A                 B              C                 D              mu(298K)        mu(1000K)    " << std::endl;
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
   
		for(unsigned int i=0;i<this->NC;i++)
		{
			double logT  = log(298.);
			double mu298 = std::exp(fittingEta(0,i)+logT*(fittingEta(1,i)+logT*(fittingEta(2,i)+logT*fittingEta(3,i))));
		
				   logT	  = log(1000.);
			double mu1000 = std::exp(fittingEta(0,i)+logT*(fittingEta(1,i)+logT*(fittingEta(2,i)+logT*fittingEta(3,i))));

			fOutput << std::right << std::setw(5) << i+1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << this->names_[i+1];

			for(int j=0;j<=3;j++)
				fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << fittingEta(j,i);
		
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << mu298;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << mu1000;
		
			fOutput << std::endl;
		}
		fOutput << std::endl;
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl << std::endl;
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteThermalConductivityFittingCoefficients(std::ostream &fOutput) const
	{
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                       THERMAL CONDUCTIVITY FITTING COEFFICIENTS    " << std::endl;
		fOutput << std::endl;
		fOutput << "           lambda = exp( A + B*logT + C*(logT)^2 + D*(logT)^3 )  [W/m/K]      " << std::endl;
		fOutput << std::endl;
		fOutput << "    Species                         A                 B              C                 D           lambda(298K)     lambda(1000K)  " << std::endl;
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
   
		for(unsigned int i=0;i<this->NC;i++)
		{
			double logT  = log(298.);
			double mu298 = std::exp(fittingLambda(0,i)+logT*(fittingLambda(1,i)+logT*(fittingLambda(2,i)+logT*fittingLambda(3,i))));
		
				   logT	  = log(1000.);
			double mu1000 = std::exp(fittingLambda(0,i)+logT*(fittingLambda(1,i)+logT*(fittingLambda(2,i)+logT*fittingLambda(3,i))));

			fOutput << std::right << std::setw(5) << i+1;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << this->names_[i+1];

			for(int j=0;j<=3;j++)
				fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << fittingLambda(j,i);
		
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << mu298;
			fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << mu1000;
		
			fOutput << std::endl;
		}
		fOutput << std::endl;
		fOutput << " ----------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl << std::endl;
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteBinaryDiffusivityFittingCoefficients(std::ostream &fOutput) const
	{
		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                     BINARY DIFFUSIVITY FITTING COEFFICIENTS                                                                                                  " << std::endl;
		fOutput << std::endl;
		fOutput << "             Djk = ( exp( A + B*logT + C*(logT)^2 + D*(logT)^3 ) ) / P    [m2/s]  (P in bar)                                                                  " << std::endl;
		fOutput << std::endl;
		fOutput << "             Species             Species                      A                 B              C                 D             Djk(298K)        Djk(1000K)    " << std::endl;
		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
 
		OpenSMOKE::OpenSMOKEVectorInt iD(this->NC*(this->NC-1)/2);
		int count=1;
		for (unsigned int i=1;i<=this->NC;i++)
			for (unsigned int k=i+1;k<=this->NC;k++)
				iD[count++]=k+(i-1)*this->NC-1;

		const double P_bar=1.;
		int *s=iD.GetHandle();
		for (unsigned int i=0;i<this->NC;i++)
		{
			for (unsigned int k=i+1;k<this->NC;k++)
			{
				double logT   = log(298.);
				double Djk298 = std::exp(fittingBinaryDiffusivities[i](0,k)+logT*(fittingBinaryDiffusivities[i](1,k)+logT*(fittingBinaryDiffusivities[i](2,k)+logT*fittingBinaryDiffusivities[i](3,k)))) / P_bar;
	
					   logT	    = log(1000.);
				double Djk1000	= std::exp(fittingBinaryDiffusivities[i](0,k)+logT*(fittingBinaryDiffusivities[i](1,k)+logT*(fittingBinaryDiffusivities[i](2,k)+logT*fittingBinaryDiffusivities[i](3,k)))) / P_bar;


				fOutput << std::right << std::setw(5) << i+1;
				fOutput << std::right << std::setw(5) << k+1;
				fOutput << std::right << std::setw(3) << "";
				fOutput << std::setw(20) << std::left  << this->names_[i+1];
				fOutput << std::setw(20) << std::left  << this->names_[k+1];

				for(int j=0;j<4;j++)
					fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << fittingBinaryDiffusivities[i](j,k);
			
				fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << Djk298;
				fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << Djk1000;
			
				fOutput << std::endl;

				s++;
			}

			fOutput << std::endl;
		}

		fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl << std::endl;
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteThermalDiffusionRatiosFittingCoefficients(std::ostream &fOutput) const
	{
		if (iThermalDiffusionRatios.Size() > 0)
		{
			fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOutput << "                   THERMAL DIFFUSION RATIOS FITTING COEFFICIENTS                                                                                                  " << std::endl;
			fOutput << std::endl;
			fOutput << "                       Teta_jk = A + B*T + C*T^2 + D*T^3   [-]                                                                  " << std::endl;
			fOutput << std::endl;
			fOutput << "             Species             Species                      A                 B              C                 D          Tetajk(298K)      Tetajk(1000K)   " << std::endl;
			fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOutput << std::endl;
 
			for (int i=1;i<=iThermalDiffusionRatios.Size();i++)
			{
				for (unsigned int k=1;k<=this->NC;k++)
				{
					unsigned int j=k+(i-1)*this->NC-1;

					double T = 298.;
					double Tetajk298	= fittingTetaBinary(0,j)+T*(fittingTetaBinary(1,j)+T*(fittingTetaBinary(2,j)+T*fittingTetaBinary(3,j)));

					T = 1000.;
					double Tetajk1000	= fittingTetaBinary(0,j)+T*(fittingTetaBinary(1,j)+T*(fittingTetaBinary(2,j)+T*fittingTetaBinary(3,j)));

					fOutput << std::right << std::setw(5) << i;
					fOutput << std::right << std::setw(5) << k;
					fOutput << std::right << std::setw(3) << "";
					fOutput << std::setw(20) << std::left  << this->names_[iThermalDiffusionRatios[i]];
					fOutput << std::setw(20) << std::left  << this->names_[k];

					for(unsigned int m=0;m<4;m++)
						fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << fittingTetaBinary(m,j);
			
					fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << Tetajk298;
					fOutput << std::setw(17) << std::scientific << std::right << std::setprecision(6) << Tetajk1000;
			
					fOutput << std::endl;
				}

				fOutput << std::endl;
			}

			fOutput << " -------------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
			fOutput << std::endl << std::endl;
		}
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteTransportDataOnXMLFile(std::stringstream& xml_string) const
	{
		xml_string << "<Transport type=\"CHEMKIN\">" << std::endl;
		WriteTransportDataOnASCIIFile(xml_string);
		xml_string << "</Transport>" << std::endl;

		xml_string << "<Lennard-Jones>" << std::endl;
		for (unsigned int j = 1; j <= this->NC; j++)
			xml_string << this->MW[j]/PhysicalConstants::Nav_kmol << " " << sigma[j] * 1.e-10 << " " << epsylon_over_kb[j] << std::endl;
		xml_string << "</Lennard-Jones>" << std::endl;

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteSpeciesBundlingOnXMLFile(std::stringstream& xml_string, const double epsilon) const
	{
		xml_string << "<SpeciesBundling>" << std::endl;
		SpeciesBundling(xml_string, epsilon);
		xml_string << "</SpeciesBundling>" << std::endl;
		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteTransportDataOnASCIIFileOldStyle(std::ofstream &fOutput) const
	{
		const unsigned int tag_transport = 1;
		const double minValue = -1;
		const double maxValue =  1;

		// Write tag for transport file
		fOutput << tag_transport << std::endl;

		// Viscosity fitting coefficients
		{
			fOutput << minValue << std::endl;
			fOutput << maxValue << std::endl;

			fOutput << fittingEta.cols() << " " << fittingEta.rows() << std::endl;
			for (unsigned int j=0;j<this->NC;j++)
			{	
				for (unsigned int i=0;i<4;i++)
					fOutput << fittingEta(i,j) << " ";
				fOutput << std::endl;
			}
		}

		// Thermal conductivity fitting coefficients
		{
			fOutput << minValue << std::endl;
			fOutput << maxValue << std::endl;

			fOutput << fittingLambda.cols() << " " << fittingLambda.rows() << std::endl;
			for (unsigned int j=0;j<this->NC;j++)
			{	
				for (unsigned int i=0;i<4;i++)
					fOutput << fittingLambda(i,j) << " ";
				fOutput << std::endl;
			}
		}

		// Mass diffusion fitting coefficients
		{
			fOutput << minValue << std::endl;
			fOutput << maxValue << std::endl;

			fOutput << fittingBinaryDiffusivities[0].cols()*this->NC << " " << fittingBinaryDiffusivities[0].rows() << std::endl;
			for (unsigned int j=0;j<this->NC;j++)
				for (unsigned int k=0;k<this->NC;k++)
				{	
					if (j==k)
					{
						for (unsigned int i=0;i<4;i++)
							fOutput << 0 << " ";
					}
					else
					{
						for (unsigned int i=0;i<4;i++)
							fOutput << fittingBinaryDiffusivities[j](i,k) << " ";
					}
					fOutput << std::endl;
				}
		}
		
		// Thermal diffusion coefficients
		{
			fOutput << minValue << std::endl;
			fOutput << maxValue << std::endl;
		
			fOutput << iThermalDiffusionRatios.Size() << std::endl;
			for (int j=1;j<=iThermalDiffusionRatios.Size();j++)
				fOutput << iThermalDiffusionRatios[j] << " ";
			fOutput << std::endl;
		
			fOutput << fittingTetaBinary.cols() << " " << fittingTetaBinary.rows() << std::endl;
			for (int j=0;j<iThermalDiffusionRatios.Size();j++)
				for (unsigned int k=0;k<this->NC;k++)
				{	
					for (unsigned int i=0;i<4;i++)
						fOutput << fittingTetaBinary(i,k+j*this->NC) << " ";
					fOutput << std::endl;
				}
		}

		return true;
	}

	template<typename Species>
	void PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::SpeciesBundling(std::ostream &fOutput, const double epsilon) const
	{
		fOutput << "<SelfDiffusion>" << std::endl;
		for (unsigned int k = 0; k<this->NC; k++)
		{
			for (unsigned int i = 0; i<4; i++)
				fOutput << fittingBinaryDiffusivities[k](i, k) << " ";
			fOutput << std::endl;
		}
		fOutput << "</SelfDiffusion>" << std::endl;

		std::cout << " * Species bundling..." << std::endl;

		std::cout << "   Calculating the error matrix" << std::endl;
		Eigen::MatrixXd eps(this->NC, this->NC);
		Eigen::MatrixXi A(this->NC, this->NC);

		// Populating the error matrix (must be done only once)
		{
			eps.setZero();

			const double Tmin = 300.;
			const double Tmax = 3000.;
			const double deltaT = 300.;
			const unsigned int n = static_cast<unsigned int>((Tmax - Tmin) / deltaT) + 1;

			for (unsigned int i = 0; i < this->NC; i++)
			for (unsigned int j = i + 1; j < this->NC; j++)
			{
				for (unsigned int m = 0; m < n; m++)
				{
					const double T = Tmin + deltaT*m;
					const double lnT = std::log(T);

					for (unsigned int k = 0; k < this->NC; k++)
					{
						double epsilon_local = 0.;
						for (unsigned int l = 0; l < 4; l++)
							epsilon_local += (fittingBinaryDiffusivities[i](l, k) - fittingBinaryDiffusivities[j](l, k)) * std::pow(lnT, double(l));
						epsilon_local = std::fabs(epsilon_local);

						if (epsilon_local > eps(i, j))
							eps(i, j) = epsilon_local;
					}
				}
			}

			for (unsigned int i = 0; i < this->NC; i++)
			for (unsigned int j = i + 1; j < this->NC; j++)
				eps(j, i) = eps(i, j);
		}

		std::vector<double> list_of_epsilon(7);
		list_of_epsilon[0] =  0.010;
		list_of_epsilon[1] =  0.025; 
		list_of_epsilon[2] =  0.050; 
		list_of_epsilon[3] =  0.075; 
		list_of_epsilon[4] =  0.100; 
		list_of_epsilon[5] =  0.250; 
		list_of_epsilon[6] =  0.500; 

		for (unsigned int s = 0; s < list_of_epsilon.size(); s++)
		{
			const double epsilon = list_of_epsilon[s];

			std::cout << "   Calculating the adiancency matrix A for epsilon " << epsilon << std::endl;
			A.setZero();
			for (unsigned int i = 0; i<this->NC; i++)
				for (unsigned int j = 0; j<this->NC; j++)
					if (eps(i, j) < epsilon)	A(i, j) = 1;

			
			std::vector<unsigned int> d(this->NC);
			std::vector<double> d_sum(this->NC);
			for (unsigned int i = 0; i < this->NC; i++)
			{
				d[i] = 0;
				d_sum[i] = 0.;

				for (unsigned int j = 0; j < this->NC; j++)
				{
					d[i] += A(i, j);
					d_sum[i] += A(i, j)*eps(i, j);
				}
			}

			// Sort the degree in descending order
			std::vector<size_t> d_sorted_indices = OpenSMOKE::SortAndTrackIndicesDecreasing(d); 

			// Additional data
			std::vector<unsigned int> d_transfer(this->NC);
			std::vector<bool> available_species(this->NC);
			for (unsigned int i = 0; i < this->NC; i++)
			{
				available_species[i] = true;
				for (unsigned int j = 0; j < this->NC; j++)
				if (i == d_sorted_indices[j])
				{
					d_transfer[i] = j;
					break;
				}
			}


			std::cout << "   Applying the Coarse Grained Algorithm (GCA)" << std::endl;

			std::vector<unsigned int> reference_species;
			std::vector< std::vector<unsigned int> > group_species;
			for (unsigned int k = 0; k < this->NC; k++)
			{
				if (available_species[k] == true)
				{
					unsigned int index = d_sorted_indices[k];

					std::vector<unsigned int> group_species_local;
					for (unsigned int i = 0; i < this->NC; i++)
					{
						if (A(index, i) == 1 && available_species[d_transfer[i]] == true)
						{
							group_species_local.push_back(i);
							available_species[d_transfer[i]] = false;
						}
					}

					reference_species.push_back(index);
					group_species.push_back(group_species_local);
				}
				else
				{
					continue;
				}
			}

			std::vector<int> group(this->NC);
			for (unsigned int k = 0; k < this->NC; k++)
			{
				group[k] = -1;
				for (unsigned int i = 0; i < reference_species.size(); i++)
				{
					for (unsigned int j = 0; j < group_species[i].size(); j++)
					if (group_species[i][j] == k)
					{
						group[k] = i;
						break;
					}
				}
			}

			// Check if every species has a group
			for (unsigned int k = 0; k < this->NC; k++)
			if (group[k] == -1)
				OpenSMOKE::FatalErrorMessage("Species bundling failed!");

			std::cout << "   Number of groups: " << reference_species.size() << "/" << this->NC << std::endl;

			// Write on XML file
			{
				fOutput << "<Bundling epsilon=\"" << epsilon << "\">" << std::endl;

				fOutput << "<NumberOfGroups>" << std::endl;
				fOutput << reference_species.size() << std::endl;
				fOutput << "</NumberOfGroups>" << std::endl;

				fOutput << "<ReferenceSpecies>" << std::endl;
				fOutput << reference_species.size() << std::endl;
				for (unsigned int i = 0; i < reference_species.size(); i++)
					fOutput << reference_species[i] << " ";
				fOutput << std::endl;
				fOutput << "</ReferenceSpecies>" << std::endl;

				fOutput << "<SpeciesInGroups>" << std::endl;
				for (unsigned int i = 0; i < reference_species.size(); i++)
				{
					fOutput << group_species[i].size() << std::endl;
					for (unsigned int j = 0; j < group_species[i].size(); j++)
						fOutput << group_species[i][j] << " ";
					fOutput << std::endl;
				}
				fOutput << "</SpeciesInGroups>" << std::endl;

				fOutput << "<GroupOfSpecies>" << std::endl;
				fOutput << this->NC << std::endl;
				for (unsigned int k = 0; k < this->NC; k++)
					fOutput << group[k] << " ";
				fOutput << std::endl;
				fOutput << "</GroupOfSpecies>" << std::endl;

				fOutput << "</Bundling>" << std::endl;
			}
		}

		std::cout << "Done" << std::endl;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteTransportDataOnASCIIFile(std::ostream &fOutput) const
	{
		unsigned int number_species_thermal_diffusion_ratios = 0;
		for (unsigned int j=0;j<this->NC;j++)
			if (this->species_[j].MolecularWeight() <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)
				number_species_thermal_diffusion_ratios++;
		fOutput << number_species_thermal_diffusion_ratios << std::endl;
		
		unsigned int count_species_thermal_diffusion_ratios = 0;
		for (unsigned int j=0;j<this->NC;j++)
		{
			fOutput << this->species_[j].MolecularWeight() << std::endl;

			// Thermal conductivity fitting coefficients
			for (unsigned int i=0;i<4;i++)
				fOutput << fittingLambda(i,j) << " ";
			fOutput << std::endl;

			// Dynamic viscosity fitting coefficients
			for (unsigned int i=0;i<4;i++)
				fOutput << fittingEta(i,j) << " ";
			fOutput << std::endl;

			// Mass diffusion coefficients
			for (unsigned int k=j+1;k<this->NC;k++)
			{	
				for (unsigned int i=0;i<4;i++)
						fOutput << fittingBinaryDiffusivities[j](i,k) << " ";
				fOutput << std::endl;
			}

			// Thermal diffusion ratios
			if (this->species_[j].MolecularWeight() <= PhysicalConstants::MAX_MW_THERMALDIFFUSION_RATIOS)
			{
				count_species_thermal_diffusion_ratios++;
				for (unsigned int k=0;k<this->NC;k++)
				{	
					for (unsigned int i=0;i<4;i++)
						fOutput << fittingTetaBinary(i,k+(count_species_thermal_diffusion_ratios-1)*this->NC) << " ";
					fOutput << std::endl;
				}
			}
		}

		fOutput << "E" << std::endl;

		return true;
	}

	template<typename Species>
	bool PreProcessorSpeciesPolicy_CHEMKIN_WithTransport<Species>::WriteTransportTableOnASCIIFile(std::ostream& fOutput) const
	{
		fOutput << "------------------------------------------------------------------------------------------------------------------" << std::endl;
		fOutput << "  Species                  Shape       eps/kb          sigma             mu           alfa        zRot298         " << std::endl;
		fOutput << "------------------------------------------------------------------------------------------------------------------" << std::endl;
		for(int i=1;i<=shape_factor.Size();i++)
		{
			fOutput << std::right << std::setw(5) << i;
			fOutput << ". ";
			fOutput << std::setw(20) << std::left  << this->names_[i];
			fOutput << std::setw(3)  << std::right << shape_factor[i];
			fOutput << std::setw(15) << std::fixed << std::right << std::setprecision(4)  << epsylon_over_kb[i];
			fOutput << std::setw(15) << std::fixed << std::right << std::setprecision(4)  << sigma[i];
			fOutput << std::setw(15) << std::fixed << std::right << std::setprecision(4)  << mu[i]/1.e-6;
			fOutput << std::setw(15) << std::fixed << std::right << std::setprecision(4)  << alfa[i];
			fOutput << std::setw(15) << std::fixed << std::right << std::setprecision(4)  << zRot298[i];
			fOutput << std::endl;
		}
		fOutput << "------------------------------------------------------------------------------------------------------------------";
		fOutput << std::endl << std::endl << std::endl;
	
		return true;
	}

}

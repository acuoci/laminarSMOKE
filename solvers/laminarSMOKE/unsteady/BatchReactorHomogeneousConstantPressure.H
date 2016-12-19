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

#ifndef BatchReactorHomogeneousConstantPressure_H
#define	BatchReactorHomogeneousConstantPressure_H

class BatchReactorHomogeneousConstantPressure
{
public:

	BatchReactorHomogeneousConstantPressure(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
												OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap);

	void SetReactor( const double P0 );
	void SetTemperature( const double T, double* y );
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };
	void SetISAT( const bool flag ) { isat_ = flag; };

	void SetHMOMTest(const bool flag);
	void SetHMOMAnalysis(const bool flag);
	void SetHMOMPrimaryParticles(const double HMOM_S0, const double HMOM_V0);
	void SetHMOMOptions(const int index_H, const int index_OH, const int index_H2, const int index_H2O, const int index_C2H2, const int index_O2, const std::vector<int> index_PAHs);


	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() const;

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	void Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	double GetTemperature() const;

	void SetDRG(OpenSMOKE::DRG* drg) { drg_ = drg; drgAnalysis_ = true; }
	void SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega );

	void SetNormalizedHMOM( const double& HMOM_Norm_N0, const double& HMOM_Norm_M00, const double& HMOM_Norm_M10, const double& HMOM_Norm_M01);

private:

	double T_;
	double P0_;
	double rho_;

	double cTot_;
	double MW_;
	double CpMixMass_;
	double QR_;

	double HMOM_Norm_N0_;
	double HMOM_Norm_M00_;
	double HMOM_Norm_M10_;
	double HMOM_Norm_M01_;
	double AvogadroNumber_;
	double HMOM_S0_;
	double HMOM_V0_;

	double HMOM_N0_;
	double HMOM_M00_;
	double HMOM_M10_;
	double HMOM_M01_;

	double HMOM_N_;
	double HMOM_fv_;
	double HMOM_S_;
	double HMOM_V_;
	double HMOM_epsilon_M00_;

	unsigned int NC_;
	unsigned int NR_;
	unsigned int NHMOM_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN<double>&			kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble Rf_;
	OpenSMOKE::OpenSMOKEVectorDouble Rb_;
	OpenSMOKE::OpenSMOKEVectorDouble r_;

	bool isat_;
	bool checkMassFractions_;
	bool energyEquation_;

	double enthalpy_;

	OpenSMOKE::DRG* drg_;
	bool drgAnalysis_;

	void ProcessNormalizedHMOM();

	bool hmom_analysis_;
	bool hmom_test_;
	int hmom_index_H_;
	int hmom_index_OH_;
	int hmom_index_H2_;
	int hmom_index_H2O_;
	int hmom_index_C2H2_;
	int hmom_index_O2_;
	std::vector<int> hmom_index_PAHs_;
};

void BatchReactorHomogeneousConstantPressure::SetHMOMAnalysis(const bool flag)
{
	hmom_analysis_ = flag;
}

void BatchReactorHomogeneousConstantPressure::SetHMOMTest(const bool flag)
{
	hmom_test_ = flag;
}

void BatchReactorHomogeneousConstantPressure::SetHMOMPrimaryParticles(const double HMOM_S0, const double HMOM_V0)
{
	AvogadroNumber_ = 6.022e23;	// [#/mol]
	HMOM_S0_ = HMOM_S0;
	HMOM_V0_ = HMOM_V0;
}

void BatchReactorHomogeneousConstantPressure::SetHMOMOptions(const int index_H, const int index_OH, const int index_H2, const int index_H2O, const int index_C2H2, const int index_O2, const std::vector<int> index_PAHs)
{
	hmom_index_H_ = index_H;
	hmom_index_OH_ = index_OH;
	hmom_index_H2_ = index_H2;
	hmom_index_H2O_ = index_H2O;
	hmom_index_C2H2_ = index_C2H2;
	hmom_index_O2_ = index_O2;
	hmom_index_PAHs_ = index_PAHs;
}

unsigned int BatchReactorHomogeneousConstantPressure::NumberOfEquations() const
{
	if (isat_ == false)	return NC_+1+NHMOM_;
	else			return NC_+2+NHMOM_;
}

BatchReactorHomogeneousConstantPressure::BatchReactorHomogeneousConstantPressure
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NHMOM_ = 4;
		NC_    = thermodynamicsMap_.NumberOfSpecies();
		NR_    = kineticsMap_.NumberOfReactions();
		QR_    = 0.;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		ChangeDimensions(NC_, &Rf_, true);
		ChangeDimensions(NC_, &Rb_, true);
		ChangeDimensions(NR_, &r_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;
		isat_ = false;
		drgAnalysis_ = false;

		HMOM_epsilon_M00_ = 1.e-16;
	}

void BatchReactorHomogeneousConstantPressure::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressure::SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega )
{
	omega_ = omega;
}

void BatchReactorHomogeneousConstantPressure::SetNormalizedHMOM( const double& HMOM_Norm_N0, const double& HMOM_Norm_M00, const double& HMOM_Norm_M10, const double& HMOM_Norm_M01)
{
	HMOM_Norm_N0_  = HMOM_Norm_N0;
	HMOM_Norm_M00_ = HMOM_Norm_M00;
	HMOM_Norm_M10_ = HMOM_Norm_M10;
	HMOM_Norm_M01_ = HMOM_Norm_M01;
}

void BatchReactorHomogeneousConstantPressure::ProcessNormalizedHMOM()
{
	// Denormalize
	HMOM_N0_  = HMOM_Norm_N0_*AvogadroNumber_;
	HMOM_M00_ = HMOM_Norm_M00_*AvogadroNumber_ + HMOM_epsilon_M00_;
	HMOM_M10_ = HMOM_Norm_M10_*AvogadroNumber_*HMOM_V0_;
	HMOM_M01_ = HMOM_Norm_M01_*AvogadroNumber_*HMOM_S0_;

	// Soot variables
	HMOM_N_  = HMOM_M00_;			// [#/m3]
	HMOM_fv_ = HMOM_M10_+1.e-16;		// [#/m3]
	HMOM_S_  = HMOM_M01_/HMOM_M00_;		// [m2]
	HMOM_V_  = HMOM_M10_/HMOM_M00_;		// [m3]
}

void BatchReactorHomogeneousConstantPressure::SetTemperature( const double T, double* y)
{
	// Temperature
	T_    = T;

	// Recover mass fractions
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = y[i-1];

	// Recover normalized moments
	HMOM_Norm_N0_  = y[NC_+0];
	HMOM_Norm_M00_ = y[NC_+1];
	HMOM_Norm_M10_ = y[NC_+2];
	HMOM_Norm_M01_ = y[NC_+3];

	// Process normalized moments 
	ProcessNormalizedHMOM();

	// Mole fractions
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);

	// Enthalpy
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);
	thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(enthalpy_, x_);
	enthalpy_ /= MW_;
}

double BatchReactorHomogeneousConstantPressure::GetTemperature() const
{
	return T_;
}

int BatchReactorHomogeneousConstantPressure::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	if (drgAnalysis_ == false)
	{
		// Recover mass fractions
		if (checkMassFractions_ == true)
		{	for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = max(y[i], 0.);
		}
		else
		{
			for(unsigned int i=1;i<=NC_;++i)
				omega_[i] = y[i];
		}

		// Recover temperature
		T_ = y[NC_+1];

		// Recover normalized moments	
		HMOM_Norm_N0_  = y[NC_+2];
		HMOM_Norm_M00_ = y[NC_+3];
		HMOM_Norm_M10_ = y[NC_+4];
		HMOM_Norm_M01_ = y[NC_+5];

		// Process normalized moments 
		ProcessNormalizedHMOM();
		
		// Test HMOM library
		#include "testHMOM.H"

		// Calculates the pressure and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
		cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
		Product(cTot_, x_, &c_);
	    	rho_ = cTot_*MW_;

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P0_);

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P0_);
		kineticsMap_.KineticConstants();
		kineticsMap_.ReactionRates(c_);
		kineticsMap_.FormationRates(&R_);
	
		// Species equations
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = thermodynamicsMap_.MW()[i]*R_[i]/rho_;
		   
	    	// Energy equation
	    	dy[NC_+1] = 0.;     
	    	if (energyEquation_ == true)
	    	{
			double CpMixMolar; 
			thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
			CpMixMass_ = CpMixMolar / MW_;
			QR_ = kineticsMap_.HeatRelease(R_);
		
			dy[NC_+1]  = QR_ / (rho_*CpMixMass_);
		}

		// HMOM Equations
		if (hmom_analysis_ == true)
		{
			double Y_H    = omega_[hmom_index_H_+1];	// [-]
			double Y_OH   = omega_[hmom_index_OH_+1];	// [-]
			double C_H    = c_[hmom_index_H_+1]/1000.;	// [mol/cm3]
			double C_OH   = c_[hmom_index_OH_+1]/1000.;	// [mol/cm3]
			double C_H2O  = c_[hmom_index_H2O_+1]/1000.;	// [mol/cm3]
			double C_H2   = c_[hmom_index_H2_+1]/1000.;	// [mol/cm3]
			double C_C2H2 = c_[hmom_index_C2H2_+1]/1000.;	// [mol/cm3]
			double C_O2   = c_[hmom_index_O2_+1]/1000.;	// [mol/cm3]

			double C_PAH = 0.;
			for (unsigned int i=0;i<hmom_index_PAHs_.size();i++)
				C_PAH += c_[hmom_index_PAHs_[i]+1]/1000.;	// [mol/m3]

			double P_atm  = P0_/101325.;							// [atm]
			double mu_lam  = 1.716e-5*Foam::pow(T_/273.15,1.5)*(273.15+110.4)/(T_+110.4);		// [kg/m/s]

			double HMOM_Norm_M00_Source = 0.;	// [mol/m3/s]
			double HMOM_Norm_M10_Source = 0.;	// [mol/m3/s]
			double HMOM_Norm_M01_Source = 0.;	// [mol/m3/s]
			double HMOM_Norm_N0_Source  = 0.;	// [mol/m3/s]

			source_terms_soot_hmom_ (	&Y_H, &Y_OH, 
							&C_OH, &C_H, &C_H2O, &C_H2, &C_C2H2, &C_O2,
							&C_PAH, &T_, &P_atm, &mu_lam, 
							&HMOM_Norm_M00_, &HMOM_Norm_M10_, &HMOM_Norm_M01_, &HMOM_Norm_N0_, 
							&HMOM_Norm_M00_Source, &HMOM_Norm_M10_Source, &HMOM_Norm_M01_Source, &HMOM_Norm_N0_Source);

			// Moments equations
			dy[NC_+2] = HMOM_Norm_N0_Source;
			dy[NC_+3] = HMOM_Norm_M00_Source;
			dy[NC_+4] = HMOM_Norm_M10_Source;
			dy[NC_+5] = HMOM_Norm_M01_Source;
		}
		else
		{
			// Moments equations
			dy[NC_+2] = 0.;
			dy[NC_+3] = 0.;
			dy[NC_+4] = 0.;
			dy[NC_+5] = 0.;
		}

		// ISAT Equation
		if (isat_ == true)
			dy[NC_+2] = 0.;
	}
	else
	{
		// Recover mass fractions
		if (checkMassFractions_ == true)
		{	
			for (unsigned int i=0;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omega_[j] = max(y[i+1], 0.);
			}	
		}
		else
		{
			for (unsigned int i=0;i<drg_->number_important_species();++i)	
			{
				const unsigned int j = drg_->indices_important_species()[i]+1;
				omega_[j] = y[i+1];
			}
		}

		// Recover temperature
		unsigned int index_T = drg_->number_important_species()+1;
		T_ = y[index_T];

		// Calculates the pressure and the concentrations of species
		thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
		cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
	    	rho_ = cTot_*MW_;
		Product(cTot_, x_, &c_);

		// Calculates thermodynamic properties
		thermodynamicsMap_.SetTemperature(T_);
		thermodynamicsMap_.SetPressure(P0_);

		// Calculates kinetics
		kineticsMap_.SetTemperature(T_);
		kineticsMap_.SetPressure(P0_);
		kineticsMap_.KineticConstants();
		kineticsMap_.ReactionRates(c_);

		// Reaction rates
		kineticsMap_.GetReactionRates(&r_);

		// Remove useless reactions
		for (unsigned int i=0;i<drg_->indices_unimportant_reactions().size();++i)
			r_[drg_->indices_unimportant_reactions()[i]+1] = 0.;

		// Formation rates
		kineticsMap_.stoichiometry().FormationRatesFromReactionRates(&R_, r_);

		// Equations
		for (unsigned int i=0;i<drg_->number_important_species();++i)	
		{
			const unsigned int j = drg_->indices_important_species()[i]+1;
			dy[i+1] = thermodynamicsMap_.MW()[j]*R_[j]/rho_;
		}	

	    	// Energy equation
	    	dy[index_T] = 0.;     
	    	if (energyEquation_ == true)
	    	{
			double CpMixMolar; 
			thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
			CpMixMass_ = CpMixMolar / MW_;
			QR_ = kineticsMap_.HeatRelease(R_);
			
			dy[index_T]  = QR_ / (rho_*CpMixMass_);
		}
	}
	
	return 0;
}

void BatchReactorHomogeneousConstantPressure::Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t)
{
	// Recover mass fractions
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = y[i-1];
	
	// Calculates the mole fractions
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);

	// Temperature from enthalpy
	if (energyEquation_ == true)
		T_ = thermodynamicsMap_.GetTemperatureFromEnthalpyAndMoleFractions(enthalpy_*MW_, P0_, x_, T_);

	// Calculates concentrations of species
	cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
    	rho_ = cTot_*MW_;
	Product(cTot_, x_, &c_);

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P0_);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_);
	kineticsMap_.FormationRates(&R_);
	kineticsMap_.ProductionAndDestructionRates(&Rf_, &Rb_);
	
	// Species equations
	for (unsigned int i=1;i<=NC_;++i)	
		q[i-1] = thermodynamicsMap_.MW()[i]*Rf_[i]/rho_;

	for (unsigned int i=1;i<=NC_;++i)	
		d[i-1] = thermodynamicsMap_.MW()[i]*Rb_[i]/rho_;
}

int BatchReactorHomogeneousConstantPressure::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHomogeneousConstantPressure_H

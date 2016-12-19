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

#ifndef BatchReactorHomogeneousConstantVolume_H
#define	BatchReactorHomogeneousConstantVolume_H

class BatchReactorHomogeneousConstantVolume
{
public:

	BatchReactorHomogeneousConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
					OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap);

	void SetReactor( const double V0, const double P0, const double rho0);
	void SetTemperature( const double T, double* y);
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

	void SetHMOMPrimaryParticles(const double HMOM_S0, const double HMOM_V0);

	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() { return NE_; }

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	void Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	double GetTemperature() const;

private:

	double T_;
	double P_;
	double V0_;
	double rho0_;

	double cTot_;
	double MW_;
	double CvMixMass_;
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
	unsigned int NE_;
	unsigned int NHMOM_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN<double>&		kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	
	bool checkMassFractions_;
	bool energyEquation_;

	double enthalpy_;

	void ProcessNormalizedHMOM();
};

BatchReactorHomogeneousConstantVolume::BatchReactorHomogeneousConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NHMOM_ = 4;
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NE_ = NC_+1+NHMOM_;;
		QR_ = 0.;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;

		HMOM_epsilon_M00_ = 1.e-16;
	}

void BatchReactorHomogeneousConstantVolume::SetHMOMPrimaryParticles(const double HMOM_S0, const double HMOM_V0)
{
	AvogadroNumber_ = 6.022e23;	// [#/mol]
	HMOM_S0_ = HMOM_S0;
	HMOM_V0_ = HMOM_V0;
}

void BatchReactorHomogeneousConstantVolume::SetReactor( const double V0, const double P0, const double rho0)
{
	P_    = P0;
	V0_   = V0;
	rho0_ = rho0;
}

void BatchReactorHomogeneousConstantVolume::ProcessNormalizedHMOM()
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


void BatchReactorHomogeneousConstantVolume::SetTemperature( const double T, double* y)
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
	thermodynamicsMap_.SetPressure(P_);
	thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(enthalpy_, x_);
	enthalpy_ /= MW_;
}

double BatchReactorHomogeneousConstantVolume::GetTemperature() const
{
	return T_;
}

int BatchReactorHomogeneousConstantVolume::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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

	// Calculates the pressure and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	cTot_ = rho0_/MW_;
	Product(cTot_, x_, &c_);
	P_ = cTot_ * PhysicalConstants::R_J_kmol * T_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P_);
	
	// Calculates kinetics
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P_);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_);
	kineticsMap_.FormationRates(&R_);
	
	// Species equations
	for (unsigned int i=1;i<=NC_;++i)	
		dy[i] = thermodynamicsMap_.MW()[i]*R_[i]/rho0_;
                
	// Energy equation
	dy[NC_+1] = 0.;   
	if (energyEquation_ == true)
	{
		double CpMixMolar; 
		thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;

		QR_ = kineticsMap_.HeatRelease(R_);
		const double sumMoleFormationRates = R_.SumElements();

		dy[NC_+1]  = (QR_ + PhysicalConstants::R_J_kmol*T_*sumMoleFormationRates) / (rho0_*CvMixMass_);
	}

	// Moments equations
	dy[NC_+2] = 0.;
	dy[NC_+3] = 0.;
	dy[NC_+4] = 0.;
	dy[NC_+5] = 0.;

	return 0;
}

void BatchReactorHomogeneousConstantVolume::Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t)
{
	// TODO
	OpenSMOKE::ErrorMessage("BatchReactorHomogeneousConstantVolume::Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t)", "It is not available (yet)");
}

int BatchReactorHomogeneousConstantVolume::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHomogeneousConstantVolume_H

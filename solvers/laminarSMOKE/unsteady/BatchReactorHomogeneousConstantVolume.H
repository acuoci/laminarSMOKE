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

	BatchReactorHomogeneousConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
					OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	void SetReactor( const double V0, const double P0, const double rho0);
	void SetTemperature( const double T, double* y);
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

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

	unsigned int NC_;
	unsigned int NE_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	
	bool checkMassFractions_;
	bool energyEquation_;

	double enthalpy_;
};

BatchReactorHomogeneousConstantVolume::BatchReactorHomogeneousConstantVolume(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
								OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NE_ = NC_+1;
		QR_ = 0.;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		
		checkMassFractions_ = false;
		energyEquation_ = true;
	}

void BatchReactorHomogeneousConstantVolume::SetReactor( const double V0, const double P0, const double rho0)
{
	P_    = P0;
	V0_   = V0;
	rho0_ = rho0;
}

void BatchReactorHomogeneousConstantVolume::SetTemperature( const double T, double* y)
{
	// Temperature
	T_    = T;

	// Recover mass fractions
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = y[i-1];

	// Mole fractions
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());

	// Enthalpy
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P_);
	enthalpy_ = thermodynamicsMap_.hMolar_Mixture_From_MoleFractions(x_.GetHandle());
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

	// Calculates the pressure and the concentrations of species
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
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
	kineticsMap_.ReactionRates(c_.GetHandle());
	kineticsMap_.FormationRates(R_.GetHandle());
	
	// Species equations
	for (unsigned int i=1;i<=NC_;++i)	
		dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho0_;
                
    // Energy equation
    dy[NC_+1] = 0.;   
    if (energyEquation_ == true)
    {
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		CvMixMass_ = (CpMixMolar - PhysicalConstants::R_J_kmol) / MW_;
	
		QR_ = kineticsMap_.HeatRelease(R_.GetHandle());
		const double sumMoleFormationRates = R_.SumElements();
		
		dy[NC_+1]  = (QR_ + PhysicalConstants::R_J_kmol*T_*sumMoleFormationRates) / (rho0_*CvMixMass_);
	}
	
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

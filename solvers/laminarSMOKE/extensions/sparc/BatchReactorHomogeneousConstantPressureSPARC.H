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

#ifndef BatchReactorHomogeneousConstantPressureSPARC_H
#define	BatchReactorHomogeneousConstantPressureSPARC_H

class BatchReactorHomogeneousConstantPressureSPARC
{
public:

	BatchReactorHomogeneousConstantPressureSPARC(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
							OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap);

	void SetReactor( const double P0, const double CpMixMass, const double MW );
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }

	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() const;

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	void SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega );

private:

	double T_;
	double P0_;

	double MW_;
	double CpMixMass_;
	double QR_;

	unsigned int NC_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::KineticsMap_CHEMKIN&		kineticsMap_;			//!< kinetic map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;

	bool energyEquation_;
	bool debug_;
};

unsigned int BatchReactorHomogeneousConstantPressureSPARC::NumberOfEquations() const
{
	return NC_+1;
}

BatchReactorHomogeneousConstantPressureSPARC::BatchReactorHomogeneousConstantPressureSPARC
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		QR_ = 0.;
		
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &R_, true);
		
		energyEquation_ = true;
		debug_ = false;
	}

void BatchReactorHomogeneousConstantPressureSPARC::SetReactor( const double P0, const double CpMixMass, const double MW)
{
	P0_        = P0;
	MW_        = MW;
	CpMixMass_ = CpMixMass;
}

void BatchReactorHomogeneousConstantPressureSPARC::SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega )
{
	omega_ = omega;
}

int BatchReactorHomogeneousConstantPressureSPARC::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
{
	// Recover mass fractions
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = y[i];

	// Recover temperature
	T_ = y[NC_+1];

	// Calculates the pressure and the concentrations of species
	for (unsigned int i=1;i<=NC_;++i)
		x_[i] = omega_[i]*MW_/thermodynamicsMap_.MW(i-1);

	// Concentrations and density
	const double cTot = P0_/PhysicalConstants::R_J_kmol/T_;
	Product(cTot, x_, &c_);
    	const double rho = cTot*MW_;

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T_);
	thermodynamicsMap_.SetPressure(P0_);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T_);
	kineticsMap_.SetPressure(P0_);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_.GetHandle(), cTot);
	kineticsMap_.FormationRates(R_.GetHandle());

	// Species equations
	for (unsigned int i=1;i<=NC_;++i)	
		dy[i] = thermodynamicsMap_.MW(i-1)*R_[i]/rho;
	   
    	// Energy equation
    	dy[NC_+1] = 0.;     
    	if (energyEquation_ == true)
    	{
		QR_ = kineticsMap_.HeatRelease(R_.GetHandle());
		dy[NC_+1]  = QR_ / (rho*CpMixMass_);
	}

	if (debug_ == true)
	{
		Info << " * T: " << T_ << endl;
		Info << " * P: " << P0_ << endl;
		for(unsigned int i=0;i<NC_;i++)
			Info << " * " << thermodynamicsMap_.NamesOfSpecies()[i] << ": " << omega_[i+1] << endl;

		Info << " * cTot: " << cTot << endl;
		Info << " * xsum: " << x_.SumElements() << endl;
		Info << " * MW: " << MW_ << endl;

		double sum = 0.;
		for(unsigned int i=0;i<NC_;i++)
		{
			Info << " * " << thermodynamicsMap_.NamesOfSpecies()[i] << ": " << thermodynamicsMap_.MW(i)*R_[i+1] << endl;
			sum += thermodynamicsMap_.MW(i)*R_[i+1];
		}
		Info << " * sum: " << sum << endl;

		Info << " * Cp: " << CpMixMass_ << endl;
		Info << " * QR: " << QR_ << endl;
	}
	
	return 0;
}

int BatchReactorHomogeneousConstantPressureSPARC::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHomogeneousConstantPressureSPARC_H

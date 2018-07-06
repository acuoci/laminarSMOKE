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

#ifndef BatchReactorHomogeneousConstantPressureVirtualChemistry_H
#define	BatchReactorHomogeneousConstantPressureVirtualChemistry_H

class BatchReactorHomogeneousConstantPressureVirtualChemistry
{
public:

	BatchReactorHomogeneousConstantPressureVirtualChemistry(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, 
									OpenSMOKE::VirtualChemistry& vc);

	void SetReactor( const double P0 );
	void SetTemperature( const double T, double* y );
	
	void SetEnergyEquation(const bool flag) { energyEquation_ = flag; }
	
	void SetCheckMassFractions( const bool flag ) { checkMassFractions_ = flag; };

	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() const;

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	void Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	double GetTemperature() const;

	void SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega );

private:

	double T_;
	double P0_;
	double rho_;

	double cTot_;
	double MW_;
	double CpMixMass_;
	double QR_;

	unsigned int NC_;

	OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map
	OpenSMOKE::VirtualChemistry& 		vc_;				//!< virtual chemistry map

	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;

	bool checkMassFractions_;
	bool energyEquation_;
	bool debug_;
};

unsigned int BatchReactorHomogeneousConstantPressureVirtualChemistry::NumberOfEquations() const
{
	return NC_+1;
}

BatchReactorHomogeneousConstantPressureVirtualChemistry::BatchReactorHomogeneousConstantPressureVirtualChemistry
(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::VirtualChemistry& vc ):
thermodynamicsMap_(thermodynamicsMap), 
vc_(vc)
{
	NC_ = thermodynamicsMap_.NumberOfSpecies();
	QR_ = 0.;
	
	ChangeDimensions(NC_, &omega_, true);
	ChangeDimensions(NC_, &R_, true);
	
	checkMassFractions_ = false;
	energyEquation_ = true;
	debug_ = false;
}

void BatchReactorHomogeneousConstantPressureVirtualChemistry::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressureVirtualChemistry::SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega )
{
	omega_ = omega;
}

void BatchReactorHomogeneousConstantPressureVirtualChemistry::SetTemperature( const double T, double* y)
{
	Info << "Error: BatchReactorHomogeneousConstantPressureVirtualChemistry::SetTemperature not available" << endl;
	abort();
}

double BatchReactorHomogeneousConstantPressureVirtualChemistry::GetTemperature() const
{
	return T_;
}

int BatchReactorHomogeneousConstantPressureVirtualChemistry::Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy)
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

		// Calculates the density [kg/m3]
		cTot_ = P0_/PhysicalConstants::R_J_kmol/T_;
		MW_   = vc_.MWMix(omega_.GetHandle());
		rho_  = cTot_*MW_;

		// Calculates the formation rates of species [kg/m3/s]
		vc_.FormationRates(cTot_, MW_, T_, omega_.GetHandle(), R_.GetHandle());
		
		// Species equations
		for (unsigned int i=1;i<=NC_;++i)	
			dy[i] = R_[i]/rho_;
		   
	    	// Energy equation
	    	dy[NC_+1] = 0.;     
	    	if (energyEquation_ == true)
	    	{
			// Calculates the specific heat [J/kg/K]
	    		CpMixMass_ = vc_.CpMix(T_,P0_,omega_.GetHandle());

			// Calculates the reaction heat [W/m3]
			QR_ = vc_.Qdot(T_,P0_,R_.GetHandle());

			// Energy equation		
			dy[NC_+1]  = QR_ / (rho_*CpMixMass_);
		}

		if (debug_ == true && T_ > 1500.)
		{
			Info << " * T: " << T_ << endl;
			Info << " * P: " << P0_ << endl;
			for(unsigned int i=0;i<NC_;i++)
				Info << " * " << thermodynamicsMap_.NamesOfSpecies()[i] << ": " << omega_[i+1] << endl;

			Info << " * cTot: " << cTot_ << endl;
			Info << " * MW: " << MW_ << endl;
			Info << " * rho: " << rho_ << endl;

			double sum = 0.;
			for(unsigned int i=0;i<NC_;i++)
			{
				Info << " * " << thermodynamicsMap_.NamesOfSpecies()[i] << ": " << R_[i+1] << endl;
				sum += R_[i+1];
			}
			Info << " * sum: " << sum << endl;

			Info << " * Cp: " << CpMixMass_ << endl;
			Info << " * QR: " << QR_ << endl;

			getchar();
		}
	
	return 0;
}

void BatchReactorHomogeneousConstantPressureVirtualChemistry::Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t)
{
	Info << "BatchReactorHomogeneousConstantPressureVirtualChemistry::Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t) not available" << endl;
	abort();
}

int BatchReactorHomogeneousConstantPressureVirtualChemistry::Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y)
{
}

#endif // BatchReactorHomogeneousConstantPressureVirtualChemistry_H

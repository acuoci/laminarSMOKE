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

	const OpenSMOKE::OpenSMOKEVectorDouble& R() const { return R_; }
	double QR() const { return QR_; }

	unsigned int NumberOfEquations() const;

	virtual int Equations(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y, OpenSMOKE::OpenSMOKEVectorDouble& dy);
	void Equations(double* y, std::vector<double>& q, std::vector<double>& d, const double t);

	virtual int Print(const double t, const OpenSMOKE::OpenSMOKEVectorDouble& y);

	double GetTemperature() const;

	void SetDRG(OpenSMOKE::DRG* drg) { drg_ = drg; drgAnalysis_ = true; }
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
	unsigned int NR_;

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
};

unsigned int BatchReactorHomogeneousConstantPressure::NumberOfEquations() const
{
	if (isat_ == false)	return NC_+1;
	else			return NC_+2;
}

BatchReactorHomogeneousConstantPressure::BatchReactorHomogeneousConstantPressure
(	OpenSMOKE::ThermodynamicsMap_CHEMKIN<double>& thermodynamicsMap, OpenSMOKE::KineticsMap_CHEMKIN<double>& kineticsMap ):
	thermodynamicsMap_(thermodynamicsMap), 
	kineticsMap_(kineticsMap)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NR_ = kineticsMap_.NumberOfReactions();
		QR_ = 0.;
		
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
	}

void BatchReactorHomogeneousConstantPressure::SetReactor( const double P0 )
{
	P0_    = P0;
}

void BatchReactorHomogeneousConstantPressure::SetMassFractions( const OpenSMOKE::OpenSMOKEVectorDouble& omega )
{
	omega_ = omega;
}

void BatchReactorHomogeneousConstantPressure::SetTemperature( const double T, double* y)
{
	// Temperature
	T_    = T;

	// Recover mass fractions
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = y[i-1];

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

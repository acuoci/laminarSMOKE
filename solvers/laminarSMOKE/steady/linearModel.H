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

#include "utilities/soot/hmom/HMOM.h"

class linearModel
{

public:

	linearModel(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_)
	{
		NHMOM_ = 4;
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NE_ = NC_+1+NHMOM_;
		ChangeDimensions(NC_, &omega_, true);
		ChangeDimensions(NC_, &x_, true);
		ChangeDimensions(NC_, &R_, true);
		ChangeDimensions(NC_, &c_, true);
		ChangeDimensions(NC_, &U_, true);
		ChangeDimensions(NC_, &u_, true);
		ChangeDimensions(NE_, &y_plus_, true);
		ChangeDimensions(NE_, &dy_plus_, true);
		ChangeDimensions(NE_, &dy_original_, true);
		Jdiagonal_.resize(NE_);

		pah_gas_consumption_ = true;
		hmom_analysis_ = false;
		hmom_test_ = false;
	}

	void reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
					const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S);

	void reactionJacobian( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       		const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       		Eigen::VectorXd &J);

	void reactionJacobianSparse( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       			const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       			Eigen::VectorXd &J, const bool energyEquation, const bool speciesEquations);

	
	void SetHMOMAnalysis(const bool flag);

	void SetHMOMTest(const bool flag);

	void SetHMOM(OpenSMOKE::HMOM* hmom);

	void SetHMOMOptions(const int index_H, const int index_OH, const int index_H2, const int index_H2O, const int index_C2H2, const int index_O2, const std::vector<int> index_PAH);

	void SetPAHGasConsumption(const bool pah_gas_consumption);

private:

	unsigned int NC_;
	unsigned int NE_;
	unsigned int NHMOM_; 
	OpenSMOKE::OpenSMOKEVectorDouble omega_;
	OpenSMOKE::OpenSMOKEVectorDouble x_;
	OpenSMOKE::OpenSMOKEVectorDouble R_;
	OpenSMOKE::OpenSMOKEVectorDouble c_;
	OpenSMOKE::OpenSMOKEVectorDouble U_;
	OpenSMOKE::OpenSMOKEVectorDouble u_;

	OpenSMOKE::OpenSMOKEVectorDouble y_plus_;
     	OpenSMOKE::OpenSMOKEVectorDouble dy_plus_;
     	OpenSMOKE::OpenSMOKEVectorDouble dy_original_;

	Eigen::VectorXd Jdiagonal_;

	bool hmom_analysis_;
	bool hmom_test_;
	int hmom_index_H_;
	int hmom_index_OH_;
	int hmom_index_H2_;
	int hmom_index_H2O_;
	int hmom_index_C2H2_;
	int hmom_index_O2_;
	std::vector<int> hmom_index_PAH_;

	bool pah_gas_consumption_;

	OpenSMOKE::HMOM* hmom_;
};

void linearModel::SetHMOMAnalysis(const bool flag)
{
	hmom_analysis_ = flag;
}

void linearModel::SetHMOMTest(const bool flag)
{
	hmom_test_ = flag;
}

void linearModel::SetHMOM(OpenSMOKE::HMOM* hmom)
{
	hmom_ = hmom;
}

void linearModel::SetHMOMOptions(const int index_H, const int index_OH, const int index_H2, const int index_H2O, const int index_C2H2, const int index_O2, const std::vector<int> index_PAH)
{
	hmom_index_H_    = index_H;
	hmom_index_OH_   = index_OH;
	hmom_index_H2_   = index_H2;
	hmom_index_H2O_  = index_H2O;
	hmom_index_C2H2_ = index_C2H2;
	hmom_index_O2_   = index_O2;
	hmom_index_PAH_  = index_PAH;
}

void linearModel::SetPAHGasConsumption(const bool pah_gas_consumption)
{
	pah_gas_consumption_ = pah_gas_consumption;
}

// 
void linearModel::reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
				const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S)
{
	
	
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = max(y[i], 0.);
	const double T = y[NC_+1];

	// Recover normalized moments	
	const double HMOM_Norm_M00_ = y[NC_+2];
	const double HMOM_Norm_M10_ = y[NC_+3];
	const double HMOM_Norm_M01_ = y[NC_+4];
	const double HMOM_Norm_N0_  = y[NC_+5];

	// Calculates the pressure and the concentrations of species
	double MW_ = 0.;
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_, MW_, omega_);
	const double cTot_ = P0/PhysicalConstants::R_J_kmol/T;
    	const double rho_ = cTot_*MW_;
	Product(cTot_, x_, &c_);

	// Calculates thermodynamic properties
	thermodynamicsMap_.SetTemperature(T);
	thermodynamicsMap_.SetPressure(P0);

	// Calculates kinetics
	kineticsMap_.SetTemperature(T);
	kineticsMap_.SetPressure(P0);
	kineticsMap_.KineticConstants();
	kineticsMap_.ReactionRates(c_);
	kineticsMap_.FormationRates(&R_);

	for (unsigned int i=1;i<=NC_;++i)	
		S[i] = R_[i]*thermodynamicsMap_.MW()[i];

	// Temperature
    	{
		double CpMixMolar; 
		thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(CpMixMolar,x_);
		const double CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_);
		
		S[NC_+1]  = QR_;
	}

	// HMOM Equations
	if (hmom_analysis_ == true)
	{
		const double viscosity = 1.716e-5*Foam::pow(T/273.15,1.5)*(273.15+110.4)/(T+110.4);		// [kg/m/s]

		double conc_PAH = 0.;
		for (unsigned int j = 0; j<hmom_->pah_species().size(); j++)
			conc_PAH += cTot_*x_[hmom_index_PAH_[j] + 1];

		// Setup
		hmom_->SetNormalizedMoments(HMOM_Norm_M00_, HMOM_Norm_M10_, HMOM_Norm_M01_, HMOM_Norm_N0_);
		hmom_->SetTemperatureAndPressure(T, P0);
		hmom_->SetMassFractions(omega_[hmom_index_OH_+1], omega_[hmom_index_H_+1]);
		hmom_->SetConcentrations("kmol/m3", c_[hmom_index_OH_+1], c_[hmom_index_H_+1], c_[hmom_index_H2O_+1], c_[hmom_index_H2_+1], c_[hmom_index_C2H2_+1], c_[hmom_index_O2_+1], conc_PAH);
		hmom_->SetViscosity(viscosity);
	
		// Calculates
		hmom_->CalculateSourceMoments();

		// Moments equations [mol/m3/s * kg/m3]
		S[NC_+2] = rho_*hmom_->SourceM00();
		S[NC_+3] = rho_*hmom_->SourceM10();
		S[NC_+4] = rho_*hmom_->SourceM01();
		S[NC_+5] = rho_*hmom_->SourceN0();

		// PAH consumption
		if (pah_gas_consumption_ == true)
		{
			if (conc_PAH > 1.e-64)
			{
				const double R_PAH = hmom_->PAHConsumptionRate() / 1000.;	// [kmol/m3/s]
				
				double Omega_PAH = 0.;	// [kg/m3]
				for (unsigned int j = 0; j < hmom_->pah_species().size(); j++)
				{
					const double fraction = cTot_*x_[hmom_index_PAH_[j] + 1] / conc_PAH;
					const double omega = thermodynamicsMap_.MW()[hmom_index_PAH_[j] + 1] * R_PAH * fraction;
					S[hmom_index_PAH_[j] + 1] -= omega;
					Omega_PAH += omega;
				}

				int n2Index  = thermodynamicsMap_.IndexOfSpecies("N2");
				S[n2Index]  += Omega_PAH;
			}
		}
	}
	else
	{
		// Moments equations
		S[NC_+2] = 0.;
		S[NC_+3] = 0.;
		S[NC_+4] = 0.;
		S[NC_+5] = 0.;
	}
}

void linearModel::reactionJacobian( OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       Eigen::VectorXd &J) 
{
     // Calculated as suggested by Buzzi (private communication)
     const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
     const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
     const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
     const double TOLR = 1.e-7;
     const double TOLA = 1.e-12;

     for(unsigned int i=1;i<=NE_;i++)
		y_plus_[i] = y[i];

     // Call equations
     reactionSourceTerms(thermodynamicsMap_, kineticsMap_, y, P0, dy_original_);     

     // Derivatives with respect to y[kd]
     for(int kd=1;kd<=NE_;kd++)
     {
         double hf = 1.e0;
         double error_weight = 1./(TOLA+TOLR*fabs(y[kd]));
         double hJ = ETA2 * fabs(std::max(y[kd], 1./error_weight));
         double hJf = hf/error_weight;
         hJ = std::max(hJ, hJf);
         hJ = std::max(hJ, ZERO_DER);

         // This is what is done by Buzzi
         double dy = std::min(hJ, 1.e-3 + 1e-3*fabs(y[kd]));
         double udy = 1. / dy;
         y_plus_[kd] += dy;

	 reactionSourceTerms(thermodynamicsMap_, kineticsMap_, y_plus_, P0, dy_plus_);

         J(kd-1) = (dy_plus_[kd]-dy_original_[kd]) * udy;

         y_plus_[kd] = y[kd];
     }
 }

void linearModel::reactionJacobianSparse( OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       Eigen::VectorXd &J, const bool energyEquation, const bool speciesEquations) 
{
	J.setConstant(0.);

	// Sparse Jacobian
	if (speciesEquations == true)
	{
		// Recover mass fractions
		for(unsigned int i=1;i<=NC_;++i)
			omega_[i] = max(y[i], 0.);
		const double T = y[NC_+1];

		kineticsMap_.jacobian_sparsity_pattern_map()->SetEpsilon(1e-13); 
		kineticsMap_.jacobian_sparsity_pattern_map()->Jacobian(omega_, T, P0, Jdiagonal_);
		for(unsigned int i=0;i<NC_;++i)
			J(i) = Jdiagonal_(i) * thermodynamicsMap_.MW()[i+1];
	}

	// Only with respect to the temperature
	if (energyEquation == true)
	{
		// Calculated as suggested by Buzzi (private communication)
		const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
		const double ETA2 = std::sqrt(OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE);
		const double BETA = 1.e3 * OpenSMOKE::OPENSMOKE_MACH_EPS_DOUBLE;
		const double TOLR = 1.e-7;
		const double TOLA = 1.e-12;

		for(unsigned int i=1;i<=NE_;i++)
		y_plus_[i] = y[i];

		// Call equations
		reactionSourceTerms(thermodynamicsMap_, kineticsMap_, y, P0, dy_original_);     

		// Derivatives with respect to the temperature
		const int kd=NE_;
		{
			double hf = 1.e0;
			double error_weight = 1./(TOLA+TOLR*fabs(y[kd]));
			double hJ = ETA2 * fabs(std::max(y[kd], 1./error_weight));
			double hJf = hf/error_weight;
			hJ = std::max(hJ, hJf);
			hJ = std::max(hJ, ZERO_DER);

			// This is what is done by Buzzi
			double dy = std::min(hJ, 1.e-3 + 1e-3*fabs(y[kd]));
			double udy = 1. / dy;
			y_plus_[kd] += dy;

			reactionSourceTerms(thermodynamicsMap_, kineticsMap_, y_plus_, P0, dy_plus_);

			J(kd-1) = (dy_plus_[kd]-dy_original_[kd]) * udy;

			y_plus_[kd] = y[kd];
		}
	}
 }

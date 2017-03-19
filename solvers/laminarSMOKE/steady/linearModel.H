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

class linearModel
{

public:

	linearModel(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_)
	{
		NC_ = thermodynamicsMap_.NumberOfSpecies();
		NE_ = NC_+1;
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
	}

	void reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
					const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S);

	void reactionJacobian( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       		const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       		Eigen::VectorXd &J);

	void reactionJacobianSparse( 	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
		       			const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0,
		       			Eigen::VectorXd &J, const bool energyEquation, const bool speciesEquations);

private:

	unsigned int NC_;
	unsigned int NE_; 
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
};

// 
void linearModel::reactionSourceTerms(	OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap_, OpenSMOKE::KineticsMap_CHEMKIN& kineticsMap_,
				const OpenSMOKE::OpenSMOKEVectorDouble& y, const double P0, OpenSMOKE::OpenSMOKEVectorDouble& S)
{
	
	
	for(unsigned int i=1;i<=NC_;++i)
		omega_[i] = max(y[i], 0.);
	const double T = y[NC_+1];

	// Calculates the pressure and the concentrations of species
	double MW_ = 0.;
	thermodynamicsMap_.MoleFractions_From_MassFractions(x_.GetHandle(), MW_, omega_.GetHandle());
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
	kineticsMap_.ReactionRates(c_.GetHandle());
	kineticsMap_.FormationRates(R_.GetHandle());

	for (unsigned int i=1;i<=NC_;++i)	
		S[i] = R_[i]*thermodynamicsMap_.MW(i-1);

	
    	{
		const double CpMixMolar = thermodynamicsMap_.cpMolar_Mixture_From_MoleFractions(x_.GetHandle());
		const double CpMixMass_ = CpMixMolar / MW_;
		const double QR_ = kineticsMap_.HeatRelease(R_.GetHandle());
		
		S[NC_+1]  = QR_;
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
		kineticsMap_.jacobian_sparsity_pattern_map()->Jacobian(omega_.GetHandle(), T, P0, Jdiagonal_);
		for(unsigned int i=0;i<NC_;++i)
			J(i) = Jdiagonal_(i) * thermodynamicsMap_.MW(i);
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

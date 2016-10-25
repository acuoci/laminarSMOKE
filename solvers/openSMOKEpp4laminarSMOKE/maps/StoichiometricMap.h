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

#ifndef OpenSMOKE_StoichiometricMap_H
#define OpenSMOKE_StoichiometricMap_H

#include "math/OpenSMOKEClass.hpp"
#include "math/OpenSMOKEVector.h"
#include <Eigen/Sparse>

namespace OpenSMOKE
{
	class ROPA_Data
	{
	public:
		std::vector<std::vector<unsigned int> >	production_reaction_indices;
		std::vector<std::vector<unsigned int> >	destruction_reaction_indices;
		std::vector<std::vector<double> >	production_coefficients;
		std::vector<std::vector<double> >	destruction_coefficients;

		std::vector<double> production_rates;
		std::vector<double> destruction_rates;
	};

	//!  A class containing the data about the stoichiometry and the reaction orders
	/*!
		 This class provides the tools to manage the stoichiometry and the reaction orders of all the 
		 reactions included in the kinetic scheme. This class is specifically conceived in order to perform 
		 the calculations as fast as possible. To reach this goal the kinetic data are stored in a
		 format which is not so easy to manage and to understand.
	*/

	class StoichiometricMap
	{

	public:

		/**
		*@brief Default constructor
		*@param nspecies number of species 
		*@param nreactions number of reactions
		*/
		StoichiometricMap(const unsigned int nspecies, const unsigned int nreactions);
                
                /**
		*@brief Constructor with verbose option
		*@param nspecies number of species 
		*@param nreactions number of reactions
                *@param verbose show output
		*/
		StoichiometricMap(const unsigned int nspecies, const unsigned int nreactions, bool verbose);

		/**
		*@brief Initialize the stoichiometry from ASCII file (obsolete, TOREMOVE)
		*/
		void ReadFromASCIIFile(std::istream& fInput);

		/**
		*@brief Returns the change of mole for each reaction
		*/
		const OpenSMOKEVectorDouble& changeOfMoles() const { return changeOfMoles_; };

		/**
		*@brief Prints a short summary about useful information
		*/
		void Summary(std::ostream &fOut) const;

		/**
		*@brief Evaluates the equilibrium constants
		*/
		void EquilibriumConstants(OpenSMOKEVectorDouble& Kp, const OpenSMOKEVectorDouble& exp_g_over_RT, const double Patm_over_RT);

		/**
		*@brief Evaluates the rwaction enthalpies and entropies
		*/
		void ReactionEnthalpyAndEntropy(OpenSMOKEVectorDouble& reaction_dh_over_RT, OpenSMOKEVectorDouble& reaction_ds_over_R, const OpenSMOKEVectorDouble& species_h_over_RT, const OpenSMOKEVectorDouble& species_s_over_R);
		
		/**
		*@brief Evaluates the product of concentrations for each reaction (according to the stoichiometry and the reaction orders)
		*/
		void ProductOfConcentrations(OpenSMOKEVectorDouble& productDirect, OpenSMOKEVectorDouble& productReverse, const OpenSMOKEVectorDouble& c);
		
		/**
		*@brief Evaluates the formation rates from the reaction rates
		*/	
		void FormationRatesFromReactionRates(OpenSMOKEVectorDouble* R, const OpenSMOKEVectorDouble& r);

		/**
		*@brief Evaluates the Production and the destruction rates from the reaction rates
		*/
		void ProductionAndDestructionRatesFromReactionRates(OpenSMOKEVectorDouble* F, OpenSMOKEVectorDouble* D, const OpenSMOKEVectorDouble& r);

		/**
		*@brief Evaluates the Production and the destruction rates from the reaction rates (gross definition)
		*/
		void ProductionAndDestructionRatesFromReactionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D, const OpenSMOKEVectorDouble& rF, const OpenSMOKEVectorDouble& rB);

		/**
		*@brief Builds the stoichiometric matrix (sparse matrix of course!)
		*/
		void BuildStoichiometricMatrix();

		/**
		*@brief Builds the reaction order matrix (sparse matrix of course!)
		*/
		void BuildReactionOrdersMatrix();

		/**
		*@brief Internal function (TODO)
		*/
		void CompleteChangeOfMoles(OpenSMOKEVectorBool& isThermodynamicReversible);

		void RateOfProductionAnalysis(const OpenSMOKEVectorDouble& r, const bool iNormalize);
		void RateOfProductionAnalysis(const OpenSMOKEVectorDouble& rf, const OpenSMOKEVectorDouble& rb);


		void WriteRateOfProductionAnalysis(std::ostream& fout);
		void WriteRateOfProductionAnalysis(ROPA_Data& ropa);

		const Eigen::SparseMatrix<double>& stoichiometric_matrix() const { return stoichiometric_matrix_; };

		const Eigen::SparseMatrix<double>& stoichiometric_matrix_reactants() const { return stoichiometric_matrix_reactants_;};
		const Eigen::SparseMatrix<double>& stoichiometric_matrix_products() const { return stoichiometric_matrix_products_; };
		
		const Eigen::SparseMatrix<double>& reactionorders_matrix_reactants() const { return reactionorders_matrix_reactants_; };
		const Eigen::SparseMatrix<double>& reactionorders_matrix_products() const { return reactionorders_matrix_products_; };

		
		
	private:

		unsigned int number_of_species_;
		unsigned int number_of_reactions_;

		bool global_reactions_;
                
                bool verbose_output_;

		OpenSMOKEVectorUnsignedInt	numDir1,	numDir2,	numDir3,	numDir4,	numDir5;
		OpenSMOKEVectorUnsignedInt	numRevTot1, numRevTot2,	numRevTot3,	numRevTot4, numRevTot5;
		OpenSMOKEVectorUnsignedInt	numRevEq1,	numRevEq2,	numRevEq3,	numRevEq4,	numRevEq5;
		OpenSMOKEVectorUnsignedInt	jDir1,		jDir2,		jDir3,		jDir4,		jDir5;
		OpenSMOKEVectorUnsignedInt	jRevTot1,	jRevTot2,	jRevTot3,	jRevTot4,	jRevTot5;
		OpenSMOKEVectorUnsignedInt	jRevEq1,	jRevEq2,	jRevEq3,	jRevEq4,	jRevEq5;

		OpenSMOKEVectorDouble valueDir5;
		OpenSMOKEVectorDouble valueRevTot5;
		OpenSMOKEVectorDouble valueRevEq5;

		OpenSMOKEVectorUnsignedInt	lambda_jDir1, lambda_jDir2, lambda_jDir3, lambda_jDir4, lambda_jDir5;
		OpenSMOKEVectorUnsignedInt	lambda_jRevEq1, lambda_jRevEq2, lambda_jRevEq3, lambda_jRevEq4, lambda_jRevEq5;
		OpenSMOKEVectorUnsignedInt	lambda_numDir1, lambda_numDir2, lambda_numDir3, lambda_numDir4, lambda_numDir5;
		OpenSMOKEVectorUnsignedInt	lambda_numRevEq1, lambda_numRevEq2, lambda_numRevEq3, lambda_numRevEq4, lambda_numRevEq5;
		OpenSMOKEVectorDouble		lambda_valueDir5;
		OpenSMOKEVectorDouble		lambda_valueRevEq5;

		OpenSMOKEVectorDouble		changeOfMoles_;

		OpenSMOKEVectorUnsignedInt	indices_of_reactions_without_change_of_moles_;
		OpenSMOKEVectorUnsignedInt	indices_of_reactions_with_change_of_moles_plus_one_;
		OpenSMOKEVectorUnsignedInt	indices_of_reactions_with_change_of_moles_minus_one_;
		OpenSMOKEVectorUnsignedInt	indices_of_reactions_with_change_of_moles_;

		bool isTheStoichiometricMatrixAvailable_;
		bool isTheReactionOrderMatrixAvailable_;
		bool areTheContributionOfRateOfFormationMatricesAvailable_;

		Eigen::SparseMatrix<double> stoichiometric_matrix_reactants_;
		Eigen::SparseMatrix<double> stoichiometric_matrix_products_;

		Eigen::SparseMatrix<double> reactionorders_matrix_reactants_;
		Eigen::SparseMatrix<double> reactionorders_matrix_products_;

		Eigen::SparseMatrix<double> Cp;
		Eigen::SparseMatrix<double> Cd;
		Eigen::SparseMatrix<double> stoichiometric_matrix_;

	};
}

#include "StoichiometricMap.hpp"

#endif /* StoichiometricMap_H */
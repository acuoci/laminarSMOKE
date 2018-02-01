/*-----------------------------------------------------------------------*\
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
|	License                                                           |
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
		*@brief Buid vectors for non elementary reactions (if any)
		*/
		void BuildNonElementaryReactions();

		/**
		*@brief Returns the change of mole for each reaction
		*/
		const std::vector<double>& ChangeOfMoles() const { return changeOfMoles_; };

		/**
		*@brief Prints a short summary about useful information
		*/
		void Summary(std::ostream &fOut) const;

		/**
		*@brief Evaluates the equilibrium constants
		*/
		void EquilibriumConstants(double* Kp, const double* exp_g_over_RT, const double Patm_over_RT);

		/**
		*@brief Evaluates the rwaction enthalpies and entropies
		*/
		void ReactionEnthalpyAndEntropy(std::vector<double>& reaction_dh_over_RT, std::vector<double>& reaction_ds_over_R, const std::vector<double>& species_h_over_RT, const std::vector<double>& species_s_over_R);

		/**
		*@brief Evaluates the product of concentrations for each reaction (according to the stoichiometry and the reaction orders)
		*/
		void ProductOfConcentrations(std::vector<double>& productDirect, std::vector<double>& productReverse, const double* c);

		/**
		*@brief Evaluates the product of concentrations for each reaction (according to the stoichiometry and the reaction orders)
		*/
		void ProductOfConcentrationsForNonElementaryReactions(std::vector<double>& productDirect, std::vector<double>& productReverse, const double* c);

		/**
		*@brief Evaluates the formation rates from the reaction rates
		*/
		void FormationRatesFromReactionRates(double* R, const double* r);

		/**
		*@brief Evaluates the Production and the destruction rates from the reaction rates
		*/
		void ProductionAndDestructionRatesFromReactionRates(double* F, double* D, const double* r);

		/**
		*@brief Evaluates the Production and the destruction rates from the reaction rates (gross definition)
		*/
		void ProductionAndDestructionRatesFromReactionRatesGross(double* P, double* D, const double* rF, const double* rB);

		/**
		*@brief Builds the stoichiometric matrix (sparse matrix of course!)
		*/
		void BuildStoichiometricMatrix();

		/**
		*@brief Builds the reaction order matrix (sparse matrix of course!)
		*/
		void BuildReactionOrdersMatrix();

		/**
		*@brief Returns the stoichiometric matrix as a sparse matrix
		*/
		const Eigen::SparseMatrix<double>& stoichiometric_matrix() const { return stoichiometric_matrix_; };

		/**
		*@brief Returns the stoichiometric matrix of reactant species only as a sparse matrix
		*/
		const Eigen::SparseMatrix<double>& stoichiometric_matrix_reactants() const { return stoichiometric_matrix_reactants_; };

		/**
		*@brief Returns the stoichiometric matrix of product species only as a sparse matrix
		*/
		const Eigen::SparseMatrix<double>& stoichiometric_matrix_products() const { return stoichiometric_matrix_products_; };

		/**
		*@brief Returns the reaction order matrix of reactant species only as a sparse matrix
		*/
		const Eigen::SparseMatrix<double>& reactionorders_matrix_reactants() const { return reactionorders_matrix_reactants_; };

		/**
		*@brief Returns the reaction order matrix of product species only as a sparse matrix
		*/
		const Eigen::SparseMatrix<double>& reactionorders_matrix_products() const { return reactionorders_matrix_products_; };

		/**
		*@brief Return the sum of stoichiometric coefficients of reactant species (i.e. left side)
		*/
		void GetSumOfStoichiometricCoefficientsOfReactants(Eigen::VectorXd& sum_nu) const;

		/**
		*@brief Return the sum of stoichiometric coefficients of product species (i.e. right side)
		*/
		void GetSumOfStoichiometricCoefficientsOfProducts(Eigen::VectorXd& sum_nu) const;

		/**
		*@brief Internal function (TODO)
		*/
		void CompleteChangeOfMoles(const bool* isThermodynamicReversible);


	public:	// Rate of Production Analysis (ROPA) utilities

		void RateOfProductionAnalysis(const double* r, const bool iNormalize);
		void RateOfProductionAnalysis(const double* rf, const double* rb);
		void WriteRateOfProductionAnalysis(std::ostream& fout);
		void WriteRateOfProductionAnalysis(ROPA_Data& ropa);

	private:

		unsigned int number_of_species_;
		unsigned int number_of_reactions_;

		bool verbose_output_;

		std::vector<unsigned int>	numDir1_, numDir2_, numDir3_, numDir4_, numDir5_;
		std::vector<unsigned int>	numRevTot1_, numRevTot2_, numRevTot3_, numRevTot4_, numRevTot5_;
		std::vector<unsigned int>	numRevEq1_, numRevEq2_, numRevEq3_, numRevEq4_, numRevEq5_;

		std::vector<unsigned int>	jDir1_, jDir2_, jDir3_, jDir4_, jDir5_;
		std::vector<unsigned int>	jRevTot1_, jRevTot2_, jRevTot3_, jRevTot4_, jRevTot5_;
		std::vector<unsigned int>	jRevEq1_, jRevEq2_, jRevEq3_, jRevEq4_, jRevEq5_;

		std::vector<double>	valueDir5_;
		std::vector<double> valueRevTot5_;
		std::vector<double> valueRevEq5_;

		std::vector<unsigned int>	lambda_jDir1_, lambda_jDir2_, lambda_jDir3_, lambda_jDir4_, lambda_jDir5_;
		std::vector<unsigned int>	lambda_jRevEq1_, lambda_jRevEq2_, lambda_jRevEq3_, lambda_jRevEq4_, lambda_jRevEq5_;
		std::vector<unsigned int>	lambda_numDir1_, lambda_numDir2_, lambda_numDir3_, lambda_numDir4_, lambda_numDir5_;
		std::vector<unsigned int>	lambda_numRevEq1_, lambda_numRevEq2_, lambda_numRevEq3_, lambda_numRevEq4_, lambda_numRevEq5_;
		std::vector<double>	lambda_valueDir5_;
		std::vector<double>	lambda_valueRevEq5_;

		std::vector<double>			changeOfMoles_;

		std::vector<unsigned int>	indices_of_reactions_without_change_of_moles_;
		std::vector<unsigned int>	indices_of_reactions_with_change_of_moles_plus_one_;
		std::vector<unsigned int>	indices_of_reactions_with_change_of_moles_minus_one_;
		std::vector<unsigned int>	indices_of_reactions_with_change_of_moles_;

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

		unsigned int non_elementary_reactions_direct_;
		unsigned int non_elementary_reactions_reverse_;
		std::vector<bool>							is_non_elementary_reaction_direct_;
		std::vector<bool>							is_non_elementary_reaction_reverse_;
		std::vector< std::vector<unsigned int> >	non_elementary_reactions_species_indices_direct_;
		std::vector< std::vector<double> >			non_elementary_reactions_orders_direct_;
		std::vector< std::vector<unsigned int> >	non_elementary_reactions_species_indices_reverse_;
		std::vector< std::vector<double> >			non_elementary_reactions_orders_reverse_;

	};
}

#include "StoichiometricMap.hpp"

#endif /* StoichiometricMap_H */
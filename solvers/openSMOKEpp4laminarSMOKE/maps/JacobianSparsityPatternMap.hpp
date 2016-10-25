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

#include "math/OpenSMOKEUtilities.h"

namespace OpenSMOKE
{
	template<typename map>
	JacobianSparsityPatternMap<map>::JacobianSparsityPatternMap(map& kinetics_map) :
		kinetics_map_(kinetics_map)
	{
		nr = kinetics_map_.NumberOfReactions();
		nc = kinetics_map_.thermodynamics().NumberOfSpecies();

		epsilon_ = 1e-15;
		sum_mass_fractions_ = 1. + nc*epsilon_;

		ChangeDimensions(nc, &analytical_omegaStar_, true);
		ChangeDimensions(nc, &analytical_cStar_, true);
		ChangeDimensions(nc, &analytical_xStar_, true);
		ChangeDimensions(nc, &analytical_RStar_, true);
		ChangeDimensions(nr, &analytical_rf_, true);
		ChangeDimensions(nr, &analytical_rb_, true);

		// Reactants
		{
			typedef Eigen::Triplet<double> list_of_values;
			std::vector<list_of_values> tripletList;
			tripletList.reserve(nr * 3);

			for (int k = 0; k < kinetics_map_.stoichiometry().reactionorders_matrix_reactants().outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().reactionorders_matrix_reactants(), k); it; ++it)
				{
					tripletList.push_back(list_of_values(it.row(), it.col(), it.value()));
				}
			}

			drf_over_domega_ = new Eigen::SparseMatrix<double>(nr, nc);
			drf_over_domega_->setFromTriplets(tripletList.begin(), tripletList.end());
		}


		// Products
		{
			typedef Eigen::Triplet<double> list_of_values;
			std::vector<list_of_values> tripletList;
			tripletList.reserve(nr * 3);

			for (int k = 0; k < kinetics_map_.stoichiometry().reactionorders_matrix_products().outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().reactionorders_matrix_products(), k); it; ++it)
				{
					tripletList.push_back(list_of_values(it.row(), it.col(), it.value()));
				}
			}

			drb_over_domega_ = new Eigen::SparseMatrix<double>(nr, nc);
			drb_over_domega_->setFromTriplets(tripletList.begin(), tripletList.end());
		}

		// Third-body: weak coupling
		{
			std::vector<unsigned int> list_reaction;
			std::vector<unsigned int> list_species;

			kinetics_map_.WeakThirdBodyConcentrationEfficiencies(list_reaction, list_species);

			typedef Eigen::Triplet<double> list_of_values;
			std::vector<list_of_values> tripletList;
			tripletList.reserve(list_reaction.size() * 4);

			for (unsigned int k = 0; k < list_reaction.size(); k++)
				tripletList.push_back(list_of_values(list_reaction[k] - 1, list_species[k] - 1, 1.));

			dthirdbody_over_domega_ = new Eigen::SparseMatrix<double>(nr, nc);
			dthirdbody_over_domega_->setFromTriplets(tripletList.begin(), tripletList.end());


			analytical_thirdbody_reactions_.resize(dthirdbody_over_domega_->nonZeros());
			analytical_thirdbody_species_.resize(dthirdbody_over_domega_->nonZeros());

			unsigned int count = 0;
			for (int k = 0; k < dthirdbody_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dthirdbody_over_domega_, k); it; ++it)
				{
					for (unsigned int s = 1; s <= kinetics_map_.NumberOfThirdBodyReactions(); s++)
					{
						const unsigned int j = kinetics_map_.indices_of_thirdbody_reactions()[s];
						if (j == it.row() + 1)
						{
							analytical_thirdbody_reactions_(count) = s;
							break;
						}
					}

					for (int k = 1; k <= kinetics_map_.indices_of_thirdbody_species()[analytical_thirdbody_reactions_(count) - 1].Size(); k++)
					{
						if (kinetics_map_.indices_of_thirdbody_species()[analytical_thirdbody_reactions_(count) - 1][k] == it.col() + 1)
						{
							analytical_thirdbody_species_(count) = k;
							break;
						}
					}

					count++;
				}
			}
		}

		// Fall-Off: weak coupling
		{
			std::vector<unsigned int> list_reaction;
			std::vector<unsigned int> list_species;

			kinetics_map_.WeakFallOffConcentrationEfficiencies(list_reaction, list_species);

			typedef Eigen::Triplet<double> list_of_values;
			std::vector<list_of_values> tripletList;
			tripletList.reserve(list_reaction.size() * 4);

			for (unsigned int k = 0; k < list_reaction.size(); k++)
				tripletList.push_back(list_of_values(list_reaction[k] - 1, list_species[k] - 1, 1.));

			dfalloff_over_domega_ = new Eigen::SparseMatrix<double>(nr, nc);
			dfalloff_over_domega_->setFromTriplets(tripletList.begin(), tripletList.end());

			analytical_falloff_reactions_.resize(dfalloff_over_domega_->nonZeros());

			unsigned int count = 0;
			for (int k = 0; k < dfalloff_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dfalloff_over_domega_, k); it; ++it)
				{
					for (unsigned int k = 1; k <= kinetics_map_.NumberOfFallOffReactions(); k++)
					{
						const unsigned int j = kinetics_map_.indices_of_falloff_reactions()[k];
						if (j == it.row() + 1)
						{
							analytical_falloff_reactions_(count) = k;
							break;
						}
					}

					count++;
				}
			}
		}

		// CABR: weak coupling (TODO)
		{
			std::vector<unsigned int> list_reaction;
			std::vector<unsigned int> list_species;

			kinetics_map_.WeakCABRConcentrationEfficiencies(list_reaction, list_species);

			typedef Eigen::Triplet<double> list_of_values;
			std::vector<list_of_values> tripletList;
			tripletList.reserve(list_reaction.size() * 4);

			for (unsigned int k = 0; k < list_reaction.size(); k++)
				tripletList.push_back(list_of_values(list_reaction[k] - 1, list_species[k] - 1, 1.));

			dcabr_over_domega_ = new Eigen::SparseMatrix<double>(nr, nc);
			dcabr_over_domega_->setFromTriplets(tripletList.begin(), tripletList.end());
		}

		// Stoichiometric map
		Eigen::SparseMatrix<double>* stoichiometric_shadow;
		{
			typedef Eigen::Triplet<double> list_of_values;
			std::vector<list_of_values> tripletList;
			tripletList.reserve(nr * 4);

			// Reactants
			for (int k = 0; k < kinetics_map_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
				{
					tripletList.push_back(list_of_values(it.col(), it.row(), -it.value()));
				}
			}

			// Products
			for (int k = 0; k < kinetics_map_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
				{
					tripletList.push_back(list_of_values(it.col(), it.row(), it.value()));
				}
			}


			stoichiometric_shadow = new Eigen::SparseMatrix<double>(nc, nr);
			stoichiometric_shadow->setFromTriplets(tripletList.begin(), tripletList.end());
		}

		// Jacobian matrix (pattern)
		{
			jacobian_matrix_ = new Eigen::SparseMatrix<double>(nc, nc);
			*jacobian_matrix_ = (*stoichiometric_shadow) * (*drf_over_domega_ - *drb_over_domega_ + *dthirdbody_over_domega_ + *dfalloff_over_domega_ + *dcabr_over_domega_);
		}

			delete stoichiometric_shadow;
	}

	template<typename map>
	JacobianSparsityPatternMap<map>::~JacobianSparsityPatternMap()
	{
		delete drf_over_domega_;
		delete drb_over_domega_;
		delete dthirdbody_over_domega_;
		delete dfalloff_over_domega_;
		delete dcabr_over_domega_;
		delete jacobian_matrix_;		
	}

	template<typename map>
	void JacobianSparsityPatternMap<map>::SetEpsilon(const double epsilon)
	{
		epsilon_ = epsilon;
		sum_mass_fractions_ = 1. + nc*epsilon_;
	}

	template<typename map>
	void JacobianSparsityPatternMap<map>::RecognizeJacobianSparsityPattern(std::vector<unsigned int>& row, std::vector<unsigned int>& col)
	{
		row.resize(jacobian_matrix_->nonZeros());
		col.resize(jacobian_matrix_->nonZeros());
		unsigned int j = 0;
		for (int k = 0; k < jacobian_matrix_->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*jacobian_matrix_, k); it; ++it)
			{
				row[j] = it.row();
				col[j] = it.col();
				j++;
			}
		}

		std::cout << "Jacobian sparsity pattern " << std::endl;
		std::cout << " * black elements:  " << jacobian_matrix_->nonZeros() << " (" << jacobian_matrix_->nonZeros() / double(nc*nc)*100. << "%)" << std::endl;
		std::cout << " * white elements:  " << nc*nc - jacobian_matrix_->nonZeros() << " (" << (nc*nc - jacobian_matrix_->nonZeros()) / double(nc*nc)*100. << "%)" << std::endl;
		std::cout << " * blacks per row:  " << jacobian_matrix_->nonZeros() / double(nc) << std::endl;
		std::cout << std::endl;
	}

	template<typename map>
	void JacobianSparsityPatternMap<map>::Jacobian(const OpenSMOKE::OpenSMOKEVectorDouble& omega, const double T, const double P_Pa, Eigen::SparseMatrix<double> &J)
	{
		for (unsigned int i = 1; i <= nc; ++i)
			analytical_omegaStar_[i] = (omega[i] > epsilon_) ? omega[i] : epsilon_; // (omega[i] + epsilon_) / sum_mass_fractions_;

		double MWStar;
		kinetics_map_.thermodynamics().MoleFractions_From_MassFractions(analytical_xStar_, MWStar, analytical_omegaStar_);
		double cTotStar = P_Pa / PhysicalConstants::R_J_kmol / T;
		double rhoStar = cTotStar*MWStar;

		for (unsigned int i = 1; i <= nc; ++i)
			analytical_cStar_[i] = analytical_omegaStar_[i] * rhoStar / kinetics_map_.thermodynamics().MW()[i];

		// Calculates thermodynamic properties
		kinetics_map_.thermodynamics().SetTemperature(T);
		kinetics_map_.thermodynamics().SetPressure(P_Pa);

		// Calculates kinetics
		kinetics_map_.SetTemperature(T);
		kinetics_map_.SetPressure(P_Pa);
		kinetics_map_.ReactionRates(analytical_cStar_, cTotStar);
		kinetics_map_.FormationRates(&analytical_RStar_);

		// Reaction rates
		kinetics_map_.GetForwardReactionRates(&analytical_rf_);
		kinetics_map_.GetBackwardReactionRates(&analytical_rb_);

		// Reactants
		for (int k = 0; k < kinetics_map_.stoichiometry().reactionorders_matrix_reactants().outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator it_rf(*drf_over_domega_, k);
			for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().reactionorders_matrix_reactants(), k); it; ++it)
			{
				it_rf.valueRef() = it.value()*analytical_rf_[it.row() + 1] / analytical_omegaStar_[it.col() + 1];
				++it_rf;
			}
		}

		// Products
		for (int k = 0; k < kinetics_map_.stoichiometry().reactionorders_matrix_products().outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator it_rb(*drb_over_domega_, k);
			for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().reactionorders_matrix_products(), k); it; ++it)
			{
				it_rb.valueRef() = it.value()*analytical_rb_[it.row() + 1] / analytical_omegaStar_[it.col() + 1];
				++it_rb;
			}
		}

		// Third-body (pure)
		if (kinetics_map_.NumberOfThirdBodyReactions() != 0)
		{
			unsigned count = 0;
			for (int k = 0; k < dthirdbody_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dthirdbody_over_domega_, k); it; ++it)
				{
					it.valueRef() = (analytical_rf_[it.row() + 1] - analytical_rb_[it.row() + 1]) * analytical_cStar_[it.col() + 1] / analytical_omegaStar_[it.col() + 1] *
						kinetics_map_.indices_of_thirdbody_efficiencies()[analytical_thirdbody_reactions_(count) - 1][analytical_thirdbody_species_(count)] /
						kinetics_map_.Meff()[analytical_thirdbody_reactions_(count)];

					count++;
				}
			}
		}

		// Fall-off
		if (kinetics_map_.NumberOfFallOffReactions() != 0)
		{
			const double eps = 1.e-4;
			OpenSMOKE::OpenSMOKEVectorDouble cStarPlus(nc);
			cStarPlus = analytical_cStar_;

			unsigned count = 0;
			for (int k = 0; k < dfalloff_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dfalloff_over_domega_, k); it; ++it)
				{
					cStarPlus[it.col() + 1] = (analytical_omegaStar_[it.col() + 1] + eps) * rhoStar / kinetics_map_.thermodynamics().MW()[it.col() + 1];

					const double GammaFallOff = kinetics_map_.FallOffReactionsCorrection(analytical_falloff_reactions_(count), cTotStar, analytical_cStar_);
					const double GammaFallOffPlus = kinetics_map_.FallOffReactionsCorrection(analytical_falloff_reactions_(count), cTotStar, cStarPlus);
					const double dGammaFallOff = (GammaFallOffPlus - GammaFallOff) / eps;
					cStarPlus[it.col() + 1] = analytical_omegaStar_[it.col() + 1] * rhoStar / kinetics_map_.thermodynamics().MW()[it.col() + 1];

					it.valueRef() = (analytical_rf_[it.row() + 1] - analytical_rb_[it.row() + 1]) * dGammaFallOff / GammaFallOff;

					count++;
				}
			}
		}

		// CABR (TODO)
		if (kinetics_map_.NumberOfCABRReactions() != 0)
		{
			unsigned count = 0;
			for (int k = 0; k < dcabr_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dcabr_over_domega_, k); it; ++it)
				{
					it.valueRef() = (analytical_rf_[it.row() + 1] - analytical_rb_[it.row() + 1]);
					count++;
				}
			}
		}

		// Jacobian matrix (version 2)
		{
			for (int k = 0; k < J.outerSize(); ++k)
			{
				kinetics_map_.netReactionRates() = 0.;
				for (Eigen::SparseMatrix<double>::InnerIterator it(*drf_over_domega_, k); it; ++it)
					kinetics_map_.netReactionRates()[it.row() + 1] = it.value();
				for (Eigen::SparseMatrix<double>::InnerIterator it(*drb_over_domega_, k); it; ++it)
					kinetics_map_.netReactionRates()[it.row() + 1] -= it.value();

				for (Eigen::SparseMatrix<double>::InnerIterator it(*dthirdbody_over_domega_, k); it; ++it)
					kinetics_map_.netReactionRates()[it.row() + 1] += it.value();
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dfalloff_over_domega_, k); it; ++it)
					kinetics_map_.netReactionRates()[it.row() + 1] += it.value();
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dcabr_over_domega_, k); it; ++it)
					kinetics_map_.netReactionRates()[it.row() + 1] += it.value();

				kinetics_map_.FormationRates(&analytical_RStar_);

				for (Eigen::SparseMatrix<double>::InnerIterator it(J, k); it; ++it)
				{
					const int index = it.row() + 1;
					it.valueRef() = analytical_RStar_[index];// *thermodynamics_.MW()[index] / rhoStar;
				}
			}
		}

		/*
		bool output_analytical = false;
		bool output_numerical = false;

		// Output Analytical
		if (output_analytical == true)
		{
			std::ofstream fOut("Jacobian.analytical.out", std::ios::out);
			fOut.setf(std::ios::scientific);

			std::cout << "drf ana" << std::endl;
			for (int k = 0; k < drf_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*drf_over_domega_, k); it; ++it)
				{
					std::cout << " forw " << it.row() << " " << it.col() << " " << it.value() << " " << analytical_rf_[it.row() + 1] << " " << it.value()*analytical_rf_[it.row() + 1] / analytical_omegaStar_[it.col() + 1] << std::endl;;
				}
			}
			std::cout << "drb ana" << std::endl;
			for (int k = 0; k < drb_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*drb_over_domega_, k); it; ++it)
				{
					std::cout << " back " << it.row() << " " << it.col() << " " << it.value() << " " << analytical_rf_[it.row() + 1] << " " << it.value()*analytical_rf_[it.row() + 1] / analytical_omegaStar_[it.col() + 1] << std::endl;;
				}
			}

			Eigen::SparseMatrix<double> diff = *drf_over_domega_ - *drb_over_domega_;
			std::cout << "diff" << std::endl;
			for (int k = 0; k < diff.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(diff, k); it; ++it)
				{
					std::cout << " diff " << it.row() << " " << it.col() << " " << it.value() << std::endl;;
				}
			}

			for (int k = 0; k < J.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(J, k); it; ++it)
				{
					fOut << it.row() + 1 << "\t" << it.col() + 1 << "\t" << it.value() << std::endl;;
				}
			}

			fOut.close();
		}



		// Numerical calculations
		if (output_numerical == true)
		{
			const double ZERO_DER = std::sqrt(OPENSMOKE_TINY_FLOAT);
			const double ETA2 = std::sqrt(OPENSMOKE_MACH_EPS_DOUBLE);
			const double TOLR = 100. * OPENSMOKE_MACH_EPS_FLOAT;
			const double TOLA = 1.e-15;

			// Calculates thermodynamic properties
			thermodynamics_.SetTemperature(T);
			thermodynamics_.SetPressure(P_Pa);

			// Calculates kinetics
			SetTemperature(T);
			SetPressure(P_Pa);
			ReactionRates(analytical_cStar_, cTotStar);
			FormationRates(&analytical_RStar_);

			// Reaction rates
			OpenSMOKE::OpenSMOKEVectorDouble r(this->number_of_reactions_);
			GetReactionRates(&r);
			GetForwardReactionRates(&analytical_rf_);
			GetBackwardReactionRates(&analytical_rb_);

			OpenSMOKE::OpenSMOKEVectorDouble r0(this->number_of_reactions_);
			r0 = r;
			OpenSMOKE::OpenSMOKEVectorDouble RStar0(this->number_of_reactions_);
			RStar0 = analytical_RStar_;
			OpenSMOKE::OpenSMOKEVectorDouble omegaStar0 = analytical_omegaStar_;

			// calculate dr_over_domega
			Eigen::MatrixXd dr_over_domega(this->number_of_reactions_, this->number_of_species_);

			dr_over_domega.setZero();
			for (unsigned int j = 1; j <= this->number_of_species_; ++j)
			{
				double hf = 1.e0;
				double error_weight = 1. / (TOLA + TOLR*std::fabs(omegaStar0[j]));
				double hJ = ETA2 * std::fabs(std::max(omegaStar0[j], 1. / error_weight));
				double hJf = hf / error_weight;
				hJ = std::max(hJ, hJf);
				hJ = std::max(hJ, ZERO_DER);
				double dom = std::min(hJ, 1.e-3 + 1e-3*std::fabs(omegaStar0[j]));


				analytical_omegaStar_[j] = omegaStar0[j] + dom;
				for (unsigned int i = 1; i <= this->number_of_species_; ++i)
					analytical_cStar_[i] = analytical_omegaStar_[i] * rhoStar / thermodynamics_.MW()[i];

				ReactionRates(analytical_cStar_, cTotStar);
				GetReactionRates(&r);
				for (unsigned int i = 1; i <= this->number_of_reactions_; ++i)
				{
					dr_over_domega(i - 1, j - 1) = (r[i] - r0[i]) / dom;
				}
				analytical_omegaStar_[j] = omegaStar0[j];
			}

			// write dr_over_domega
			std::cout << "Matrix drdomega num" << std::endl;
			for (unsigned int i = 0; i < this->number_of_reactions_; ++i)
			{
				for (unsigned int j = 0; j < this->number_of_species_; ++j)
					std::cout << dr_over_domega(i, j) << "\t";
				std::cout << std::endl;
			}

			std::ofstream fOut("Jacobian.numerical.out", std::ios::out);
			fOut.setf(std::ios::scientific);

			// Jacobian
			Eigen::MatrixXd Jac(this->number_of_species_, this->number_of_species_);

			Jac.setZero();
			for (unsigned int j = 1; j <= this->number_of_species_; ++j)
			{
				double hf = 1.e0;
				double error_weight = 1. / (TOLA + TOLR*std::fabs(omegaStar0[j]));
				double hJ = ETA2 * std::fabs(std::max(omegaStar0[j], 1. / error_weight));
				double hJf = hf / error_weight;
				hJ = std::max(hJ, hJf);
				hJ = std::max(hJ, ZERO_DER);
				double dom = std::min(hJ, 1.e-3 + 1e-3*std::fabs(omegaStar0[j]));

				analytical_omegaStar_[j] = omegaStar0[j] + dom;
				for (unsigned int i = 1; i <= this->number_of_species_; ++i)
					analytical_cStar_[i] = analytical_omegaStar_[i] * rhoStar / thermodynamics_.MW()[i];

				ReactionRates(analytical_cStar_, cTotStar);
				FormationRates(&analytical_RStar_);

				for (unsigned int i = 1; i <= this->number_of_species_; ++i)
				{
					Jac(i - 1, j - 1) = (analytical_RStar_[i] - RStar0[i]) / dom;
					//Jac(i - 1, j - 1) *= thermodynamics_.MW()[i] / rhoStar;
				}
				analytical_omegaStar_[j] = omegaStar0[j];
			}

			for (unsigned int i = 0; i < this->number_of_species_; ++i)
			{
				for (unsigned int j = 0; j < this->number_of_species_; ++j)
				if (std::fabs(Jac(j, i)) > 0.)
					fOut << j + 1 << "\t" << i + 1 << "\t" << Jac(j, i) << std::endl;
			}

			fOut.close();
		}

		if (output_analytical == true || output_numerical == true)
			getchar();
		*/

		/*
		{

		J.setZero();
		for (unsigned int i = 0; i < this->number_of_species_; ++i)
		{
		for (unsigned int j = 0; j < this->number_of_reactions_; ++j)
		netReactionRates_[j + 1] = dr_over_domega(j, i);
		FormationRates(&analytical_RStar_);
		for (unsigned int j = 0; j < this->number_of_species_; ++j)
		J(j, i) = analytical_RStar_(j + 1)*thermodynamics_.MW()[j + 1] / rhoStar;
		}


		std::cout << "J an" << std::endl;
		for (unsigned int i = 0; i < this->number_of_species_; ++i)
		{
		for (unsigned int j = 0; j < this->number_of_species_; ++j)
		std::cout << J(i, j) << "\t";
		std::cout << std::endl;
		}
		getchar();

		}
		*/
	}

	template<typename map>
	void JacobianSparsityPatternMap<map>::Jacobian(const OpenSMOKE::OpenSMOKEVectorDouble& omega, const double T, const double P_Pa, Eigen::VectorXd &Jdiagonal)
	{
		for (unsigned int i = 1; i <= nc; ++i)
			analytical_omegaStar_[i] = (omega[i] > epsilon_) ? omega[i] : epsilon_; // (omega[i] + epsilon_) / sum_mass_fractions_;

		double MWStar;
		kinetics_map_.thermodynamics().MoleFractions_From_MassFractions(analytical_xStar_, MWStar, analytical_omegaStar_);
		double cTotStar = P_Pa / PhysicalConstants::R_J_kmol / T;
		double rhoStar = cTotStar*MWStar;

		for (unsigned int i = 1; i <= nc; ++i)
			analytical_cStar_[i] = analytical_omegaStar_[i] * rhoStar / kinetics_map_.thermodynamics().MW()[i];

		// Calculates thermodynamic properties
		kinetics_map_.thermodynamics().SetTemperature(T);
		kinetics_map_.thermodynamics().SetPressure(P_Pa);

		// Calculates kinetics
		kinetics_map_.SetTemperature(T);
		kinetics_map_.SetPressure(P_Pa);
		kinetics_map_.ReactionRates(analytical_cStar_, cTotStar);
		kinetics_map_.FormationRates(&analytical_RStar_);

		// Reaction rates
		kinetics_map_.GetForwardReactionRates(&analytical_rf_);
		kinetics_map_.GetBackwardReactionRates(&analytical_rb_);

		// Reactants
		for (int k = 0; k < kinetics_map_.stoichiometry().reactionorders_matrix_reactants().outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator it_rf(*drf_over_domega_, k);
			for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().reactionorders_matrix_reactants(), k); it; ++it)
			{
				it_rf.valueRef() = it.value()*analytical_rf_[it.row() + 1] / analytical_omegaStar_[it.col() + 1];
				++it_rf;
			}
		}

		// Products
		for (int k = 0; k < kinetics_map_.stoichiometry().reactionorders_matrix_products().outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator it_rb(*drb_over_domega_, k);
			for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().reactionorders_matrix_products(), k); it; ++it)
			{
				it_rb.valueRef() = it.value()*analytical_rb_[it.row() + 1] / analytical_omegaStar_[it.col() + 1];
				++it_rb;
			}
		}

		// Third-body (pure)
		if (kinetics_map_.NumberOfThirdBodyReactions() != 0)
		{
			unsigned count = 0;
			for (int k = 0; k < dthirdbody_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dthirdbody_over_domega_, k); it; ++it)
				{
					it.valueRef() = (analytical_rf_[it.row() + 1] - analytical_rb_[it.row() + 1]) * analytical_cStar_[it.col() + 1] / analytical_omegaStar_[it.col() + 1] *
						kinetics_map_.indices_of_thirdbody_efficiencies()[analytical_thirdbody_reactions_(count) - 1][analytical_thirdbody_species_(count)] /
						kinetics_map_.Meff()[analytical_thirdbody_reactions_(count)];

					count++;
				}
			}
		}

		// Fall-off
		if (kinetics_map_.NumberOfFallOffReactions() != 0)
		{
			const double eps = 1.e-4;
			OpenSMOKE::OpenSMOKEVectorDouble cStarPlus(nc);
			cStarPlus = analytical_cStar_;

			unsigned count = 0;
			for (int k = 0; k < dfalloff_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dfalloff_over_domega_, k); it; ++it)
				{
					cStarPlus[it.col() + 1] = (analytical_omegaStar_[it.col() + 1] + eps) * rhoStar / kinetics_map_.thermodynamics().MW()[it.col() + 1];

					const double GammaFallOff = kinetics_map_.FallOffReactionsCorrection(analytical_falloff_reactions_(count), cTotStar, analytical_cStar_);
					const double GammaFallOffPlus = kinetics_map_.FallOffReactionsCorrection(analytical_falloff_reactions_(count), cTotStar, cStarPlus);
					const double dGammaFallOff = (GammaFallOffPlus - GammaFallOff) / eps;
					cStarPlus[it.col() + 1] = analytical_omegaStar_[it.col() + 1] * rhoStar / kinetics_map_.thermodynamics().MW()[it.col() + 1];

					it.valueRef() = (analytical_rf_[it.row() + 1] - analytical_rb_[it.row() + 1]) * dGammaFallOff / GammaFallOff;

					count++;
				}
			}
		}

		// CABR (TODO)
		if (kinetics_map_.NumberOfCABRReactions() != 0)
		{
			unsigned count = 0;
			for (int k = 0; k < dcabr_over_domega_->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dcabr_over_domega_, k); it; ++it)
				{
					it.valueRef() = (analytical_rf_[it.row() + 1] - analytical_rb_[it.row() + 1]);
					count++;
				}
			}
		}

		// Jacobian matrix (version diagonal)
		Jdiagonal.setConstant(0.);
		Eigen::VectorXd net_reaction_rates(nr);
		{
			for (int k = 0; k < nc; ++k)
			{
				net_reaction_rates.setConstant(0.);
				for (Eigen::SparseMatrix<double>::InnerIterator it(*drf_over_domega_, k); it; ++it)
					net_reaction_rates(it.row()) = it.value();
				for (Eigen::SparseMatrix<double>::InnerIterator it(*drb_over_domega_, k); it; ++it)
					net_reaction_rates(it.row()) -= it.value();

				for (Eigen::SparseMatrix<double>::InnerIterator it(*dthirdbody_over_domega_, k); it; ++it)
					net_reaction_rates(it.row()) += it.value();
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dfalloff_over_domega_, k); it; ++it)
					net_reaction_rates(it.row()) += it.value();
				for (Eigen::SparseMatrix<double>::InnerIterator it(*dcabr_over_domega_, k); it; ++it)
					net_reaction_rates(it.row()) += it.value();

				for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
					Jdiagonal(k) -= it.value() * net_reaction_rates(it.row());
				for (Eigen::SparseMatrix<double>::InnerIterator it(kinetics_map_.stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
					Jdiagonal(k) += it.value() * net_reaction_rates(it.row());
			}
		}
	}

}


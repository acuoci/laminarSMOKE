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

#if __cplusplus <= 199711L
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace OpenSMOKE
{
	StoichiometricMap::StoichiometricMap(const unsigned int nspecies, const unsigned int nreactions)
	{
		number_of_species_ = nspecies;
		number_of_reactions_ = nreactions;

		verbose_output_ = true;

		isTheStoichiometricMatrixAvailable_ = false;
		isTheReactionOrderMatrixAvailable_ = false;
		areTheContributionOfRateOfFormationMatricesAvailable_ = false;
		non_elementary_reactions_direct_ = 0;
		non_elementary_reactions_reverse_ = 0;
	}

	StoichiometricMap::StoichiometricMap(const unsigned int nspecies, const unsigned int nreactions, bool verbose)
	{
		number_of_species_ = nspecies;
		number_of_reactions_ = nreactions;

		verbose_output_ = verbose;

		isTheStoichiometricMatrixAvailable_ = false;
		isTheReactionOrderMatrixAvailable_ = false;
		areTheContributionOfRateOfFormationMatricesAvailable_ = false;
		non_elementary_reactions_direct_ = 0;
		non_elementary_reactions_reverse_ = 0;
	}

	void StoichiometricMap::ReadFromASCIIFile(std::istream& fInput)
	{
		Load(numDir1_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numDir2_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numDir3_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numDir4_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numDir5_, fInput, OPENSMOKE_FORMATTED_FILE);

		Load(numRevTot1_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevTot2_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevTot3_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevTot4_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevTot5_, fInput, OPENSMOKE_FORMATTED_FILE);

		Load(numRevEq1_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevEq2_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevEq3_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevEq4_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(numRevEq5_, fInput, OPENSMOKE_FORMATTED_FILE);

		Load(jDir1_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jDir2_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jDir3_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jDir4_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jDir5_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(valueDir5_, fInput, OPENSMOKE_FORMATTED_FILE);

		#if __cplusplus > 199711L
		std::for_each(jDir1_.begin(), jDir1_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jDir2_.begin(), jDir2_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jDir3_.begin(), jDir3_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jDir4_.begin(), jDir4_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jDir5_.begin(), jDir5_.end(), [](unsigned int& v) { v -= 1; });
		#else
		for (unsigned int i = 0; i<jDir1_.size(); i++)	jDir1_[i] -= 1;
		for (unsigned int i = 0; i<jDir2_.size(); i++)	jDir2_[i] -= 1;
		for (unsigned int i = 0; i<jDir3_.size(); i++)	jDir3_[i] -= 1;
		for (unsigned int i = 0; i<jDir4_.size(); i++)	jDir4_[i] -= 1;
		for (unsigned int i = 0; i<jDir5_.size(); i++)	jDir5_[i] -= 1;
		#endif

		Load(jRevTot1_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevTot2_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevTot3_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevTot4_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevTot5_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(valueRevTot5_, fInput, OPENSMOKE_FORMATTED_FILE);

		#if __cplusplus > 199711L
		std::for_each(jRevTot1_.begin(), jRevTot1_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevTot2_.begin(), jRevTot2_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevTot3_.begin(), jRevTot3_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevTot4_.begin(), jRevTot4_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevTot5_.begin(), jRevTot5_.end(), [](unsigned int& v) { v -= 1; });
		#else
		for (unsigned int i = 0; i<jRevTot1_.size(); i++)	jRevTot1_[i] -= 1;
		for (unsigned int i = 0; i<jRevTot2_.size(); i++)	jRevTot2_[i] -= 1;
		for (unsigned int i = 0; i<jRevTot3_.size(); i++)	jRevTot3_[i] -= 1;
		for (unsigned int i = 0; i<jRevTot4_.size(); i++)	jRevTot4_[i] -= 1;
		for (unsigned int i = 0; i<jRevTot5_.size(); i++)	jRevTot5_[i] -= 1;
		#endif

		Load(jRevEq1_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevEq2_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevEq3_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevEq4_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(jRevEq5_, fInput, OPENSMOKE_FORMATTED_FILE);
		Load(valueRevEq5_, fInput, OPENSMOKE_FORMATTED_FILE);

		#if __cplusplus > 199711L
		std::for_each(jRevEq1_.begin(), jRevEq1_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevEq2_.begin(), jRevEq2_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevEq3_.begin(), jRevEq3_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevEq4_.begin(), jRevEq4_.end(), [](unsigned int& v) { v -= 1; });
		std::for_each(jRevEq5_.begin(), jRevEq5_.end(), [](unsigned int& v) { v -= 1; });
		#else
		for (unsigned int i = 0; i<jRevEq1_.size(); i++)	jRevEq1_[i] -= 1;
		for (unsigned int i = 0; i<jRevEq2_.size(); i++)	jRevEq2_[i] -= 1;
		for (unsigned int i = 0; i<jRevEq3_.size(); i++)	jRevEq3_[i] -= 1;
		for (unsigned int i = 0; i<jRevEq4_.size(); i++)	jRevEq4_[i] -= 1;
		for (unsigned int i = 0; i<jRevEq5_.size(); i++)	jRevEq5_[i] -= 1;
		#endif

		Load(changeOfMoles_, fInput, OPENSMOKE_FORMATTED_FILE);

		bool explicit_reaction_orders;
		fInput >> explicit_reaction_orders;

		if (explicit_reaction_orders == false)
		{
			lambda_numDir1_ = numDir1_;
			lambda_numDir2_ = numDir2_;
			lambda_numDir3_ = numDir3_;
			lambda_numDir4_ = numDir4_;
			lambda_numDir5_ = numDir5_;

			lambda_numRevEq1_ = numRevEq1_;
			lambda_numRevEq2_ = numRevEq2_;
			lambda_numRevEq3_ = numRevEq3_;
			lambda_numRevEq4_ = numRevEq4_;
			lambda_numRevEq5_ = numRevEq5_;

			lambda_jDir1_ = jDir1_;
			lambda_jDir2_ = jDir2_;
			lambda_jDir3_ = jDir3_;
			lambda_jDir4_ = jDir4_;
			lambda_jDir5_ = jDir5_;
			lambda_valueDir5_ = valueDir5_;

			lambda_jRevEq1_ = jRevEq1_;
			lambda_jRevEq2_ = jRevEq2_;
			lambda_jRevEq3_ = jRevEq3_;
			lambda_jRevEq4_ = jRevEq4_;
			lambda_jRevEq5_ = jRevEq5_;
			lambda_valueRevEq5_ = valueRevEq5_;
		}
		else
		{
			Load(lambda_numDir1_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numDir2_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numDir3_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numDir4_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numDir5_, fInput, OPENSMOKE_FORMATTED_FILE);

			Load(lambda_numRevEq1_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numRevEq2_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numRevEq3_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numRevEq4_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_numRevEq5_, fInput, OPENSMOKE_FORMATTED_FILE);

			Load(lambda_jDir1_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jDir2_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jDir3_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jDir4_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jDir5_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_valueDir5_, fInput, OPENSMOKE_FORMATTED_FILE);

			#if __cplusplus > 199711L
			std::for_each(lambda_jDir1_.begin(), lambda_jDir1_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jDir2_.begin(), lambda_jDir2_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jDir3_.begin(), lambda_jDir3_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jDir4_.begin(), lambda_jDir4_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jDir5_.begin(), lambda_jDir5_.end(), [](unsigned int& v) { v -= 1; });
			#else
			for (unsigned int i = 0; i<lambda_jDir1_.size(); i++)	lambda_jDir1_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jDir2_.size(); i++)	lambda_jDir2_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jDir3_.size(); i++)	lambda_jDir3_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jDir4_.size(); i++)	lambda_jDir4_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jDir5_.size(); i++)	lambda_jDir5_[i] -= 1;
			#endif

			Load(lambda_jRevEq1_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jRevEq2_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jRevEq3_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jRevEq4_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_jRevEq5_, fInput, OPENSMOKE_FORMATTED_FILE);
			Load(lambda_valueRevEq5_, fInput, OPENSMOKE_FORMATTED_FILE);

			#if __cplusplus > 199711L
			std::for_each(lambda_jRevEq1_.begin(), lambda_jRevEq1_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jRevEq2_.begin(), lambda_jRevEq2_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jRevEq3_.begin(), lambda_jRevEq3_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jRevEq4_.begin(), lambda_jRevEq4_.end(), [](unsigned int& v) { v -= 1; });
			std::for_each(lambda_jRevEq5_.begin(), lambda_jRevEq5_.end(), [](unsigned int& v) { v -= 1; });
			#else
			for (unsigned int i = 0; i<lambda_jRevEq1_.size(); i++)	lambda_jRevEq1_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jRevEq2_.size(); i++)	lambda_jRevEq2_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jRevEq3_.size(); i++)	lambda_jRevEq3_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jRevEq4_.size(); i++)	lambda_jRevEq4_[i] -= 1;
			for (unsigned int i = 0; i<lambda_jRevEq5_.size(); i++)	lambda_jRevEq5_[i] -= 1;
			#endif
		}

		BuildStoichiometricMatrix();
		BuildReactionOrdersMatrix();
		BuildNonElementaryReactions();
	}

	void StoichiometricMap::BuildNonElementaryReactions()
	{
		// Recognize non elementary reactions: direct reactions
		{
			is_non_elementary_reaction_direct_.resize(number_of_reactions_);
			std::fill(is_non_elementary_reaction_direct_.begin(), is_non_elementary_reaction_direct_.end(), false);

			double* vD5 = lambda_valueDir5_.data();
			unsigned int *jD4 = lambda_jDir4_.data();
			unsigned int *jD5 = lambda_jDir5_.data();

			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				for (unsigned int k = 0; k < lambda_numDir4_[i]; k++)
				{
					is_non_elementary_reaction_direct_[*jD4] = true;
					jD4++;
				}

				for (unsigned int k = 0; k < lambda_numDir5_[i]; k++)
				{
					if (*vD5 < 1.)
						is_non_elementary_reaction_direct_[*jD5] = true;

					jD5++;
					vD5++;
				}
			}

			// Number of non elemetary direct reactions 
			non_elementary_reactions_direct_ = std::count(is_non_elementary_reaction_direct_.begin(), is_non_elementary_reaction_direct_.end(), true);;
		}

		// Recognize non elementary reactions: reverse reactions
		{
			is_non_elementary_reaction_reverse_.resize(number_of_reactions_);
			std::fill(is_non_elementary_reaction_reverse_.begin(), is_non_elementary_reaction_reverse_.end(), false);

			unsigned int *jIE4 = lambda_jRevEq4_.data();
			unsigned int *jIE5 = lambda_jRevEq5_.data();
			double *vIE5 = lambda_valueRevEq5_.data();

			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				for (unsigned int k = 0; k < lambda_numRevEq4_[i]; k++)
				{
					is_non_elementary_reaction_reverse_[*jIE4] = true;
					jIE4++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq5_[i]; k++)
				{
					if (*vIE5 < 1.)
						is_non_elementary_reaction_reverse_[*jIE5] = true;
					jIE5++;
					vIE5++;
				}
			}

			// Number of non elemetary reverse reactions 
			non_elementary_reactions_reverse_ = std::count(is_non_elementary_reaction_reverse_.begin(), is_non_elementary_reaction_reverse_.end(), true);;
		}

		if (verbose_output_ == true)
		{
			std::cout << " * Non elementary direct reactions:  " << non_elementary_reactions_direct_ << std::endl;
			std::cout << " * Non elementary reverse reactions: " << non_elementary_reactions_reverse_ << std::endl;
		}

		// Fill the vectors 
		{
			// Fill vectors for direct reactions
			if (non_elementary_reactions_direct_ != 0)
			{
				non_elementary_reactions_species_indices_direct_.resize(number_of_reactions_);
				non_elementary_reactions_orders_direct_.resize(number_of_reactions_);

				unsigned int *jD1 = lambda_jDir1_.data();
				unsigned int *jD2 = lambda_jDir2_.data();
				unsigned int *jD3 = lambda_jDir3_.data();
				unsigned int *jD4 = lambda_jDir4_.data();
				unsigned int *jD5 = lambda_jDir5_.data();
				double* vD5 = lambda_valueDir5_.data();

				for (unsigned int i = 0; i < number_of_species_; i++)
				{
					for (unsigned int k = 0; k < lambda_numDir1_[i]; k++)
					{
						if (is_non_elementary_reaction_direct_[*jD1] == true)
						{
							non_elementary_reactions_species_indices_direct_[*jD1].push_back(i);
							non_elementary_reactions_orders_direct_[*jD1].push_back(1.);
						}
						jD1++;
					}
					for (unsigned int k = 0; k < lambda_numDir2_[i]; k++)
					{
						if (is_non_elementary_reaction_direct_[*jD2] == true)
						{
							non_elementary_reactions_species_indices_direct_[*jD2].push_back(i);
							non_elementary_reactions_orders_direct_[*jD2].push_back(2.);
						}
						jD2++;
					}
					for (unsigned int k = 0; k < lambda_numDir3_[i]; k++)
					{
						if (is_non_elementary_reaction_direct_[*jD3] == true)
						{
							non_elementary_reactions_species_indices_direct_[*jD3].push_back(i);
							non_elementary_reactions_orders_direct_[*jD3].push_back(3.);
						}
						jD3++;
					}
					for (unsigned int k = 0; k < lambda_numDir4_[i]; k++)
					{
						if (is_non_elementary_reaction_direct_[*jD4] == true)
						{
							non_elementary_reactions_species_indices_direct_[*jD4].push_back(i);
							non_elementary_reactions_orders_direct_[*jD4].push_back(0.5);
						}
						jD4++;
					}
					for (unsigned int k = 0; k < lambda_numDir5_[i]; k++)
					{
						if (is_non_elementary_reaction_direct_[*jD5] == true)
						{
							non_elementary_reactions_species_indices_direct_[*jD5].push_back(i);
							non_elementary_reactions_orders_direct_[*jD5].push_back(*vD5);
						}
						jD5++;
						vD5++;
					}
				}
			}

			// Fill vectors for reverse reactions
			if (non_elementary_reactions_reverse_ != 0)
			{
				non_elementary_reactions_species_indices_reverse_.resize(number_of_reactions_);
				non_elementary_reactions_orders_reverse_.resize(number_of_reactions_);

				unsigned int *jIE1 = lambda_jRevEq1_.data();
				unsigned int *jIE2 = lambda_jRevEq2_.data();
				unsigned int *jIE3 = lambda_jRevEq3_.data();
				unsigned int *jIE4 = lambda_jRevEq4_.data();
				unsigned int *jIE5 = lambda_jRevEq5_.data();
				double *vIE5 = lambda_valueRevEq5_.data();

				for (unsigned int i = 0; i < number_of_species_; i++)
				{
					for (unsigned int k = 0; k<lambda_numRevEq1_[i]; k++)
					{
						if (is_non_elementary_reaction_reverse_[*jIE1] == true)
						{
							non_elementary_reactions_species_indices_reverse_[*jIE1].push_back(i);
							non_elementary_reactions_orders_reverse_[*jIE1].push_back(1.);
						}
						jIE1++;
					}
					for (unsigned int k = 0; k<lambda_numRevEq2_[i]; k++)
					{
						if (is_non_elementary_reaction_reverse_[*jIE2] == true)
						{
							non_elementary_reactions_species_indices_reverse_[*jIE2].push_back(i);
							non_elementary_reactions_orders_reverse_[*jIE2].push_back(2.);
						}
						jIE2++;
					}
					for (unsigned int k = 0; k<lambda_numRevEq3_[i]; k++)
					{
						if (is_non_elementary_reaction_reverse_[*jIE3] == true)
						{
							non_elementary_reactions_species_indices_reverse_[*jIE3].push_back(i);
							non_elementary_reactions_orders_reverse_[*jIE3].push_back(3.);
						}
						jIE3++;
					}
					for (unsigned int k = 0; k<lambda_numRevEq4_[i]; k++)
					{
						if (is_non_elementary_reaction_reverse_[*jIE4] == true)
						{
							non_elementary_reactions_species_indices_reverse_[*jIE4].push_back(i);
							non_elementary_reactions_orders_reverse_[*jIE4].push_back(0.5);
						}
						jIE4++;
					}
					for (unsigned int k = 0; k<lambda_numRevEq5_[i]; k++)
					{
						if (is_non_elementary_reaction_reverse_[*jIE5] == true)
						{
							non_elementary_reactions_species_indices_reverse_[*jIE5].push_back(i);
							non_elementary_reactions_orders_reverse_[*jIE5].push_back(*vIE5);
						}
						jIE5++;
						vIE5++;
					}
				}
			}
		}
	}

	void StoichiometricMap::CompleteChangeOfMoles(const bool* isThermodynamicReversible)
	{
		unsigned number_of_thermodynamic_reversible_reactions = 0;
		unsigned number_of_reactions_without_mole_change = 0;
		unsigned number_of_reactions_with_mole_change_plus_one = 0;
		unsigned number_of_reactions_with_mole_change_minus_one = 0;
		for (unsigned int j = 1; j <= number_of_reactions_; j++)
		{
			if (isThermodynamicReversible[j - 1] == true)
			{
				number_of_thermodynamic_reversible_reactions++;
				if (changeOfMoles_[j - 1] == 0.)		number_of_reactions_without_mole_change++;
				else if (changeOfMoles_[j - 1] == 1.)	number_of_reactions_with_mole_change_plus_one++;
				else if (changeOfMoles_[j - 1] == -1.)	number_of_reactions_with_mole_change_minus_one++;
			}
		}


		unsigned number_of_reactions_with_mole_change_other = number_of_thermodynamic_reversible_reactions - number_of_reactions_without_mole_change -
			number_of_reactions_with_mole_change_plus_one - number_of_reactions_with_mole_change_minus_one;
		indices_of_reactions_without_change_of_moles_.resize(number_of_reactions_without_mole_change);
		indices_of_reactions_with_change_of_moles_plus_one_.resize(number_of_reactions_with_mole_change_plus_one);
		indices_of_reactions_with_change_of_moles_minus_one_.resize(number_of_reactions_with_mole_change_minus_one);
		indices_of_reactions_with_change_of_moles_.resize(number_of_reactions_with_mole_change_other);

		number_of_reactions_without_mole_change = 0;
		number_of_reactions_with_mole_change_plus_one = 0;
		number_of_reactions_with_mole_change_minus_one = 0;
		number_of_reactions_with_mole_change_other = 0;
		for (unsigned int j = 0; j<number_of_reactions_; j++)
		{
			if (isThermodynamicReversible[j] == true)
			{
				if (changeOfMoles_[j] == 0.)	indices_of_reactions_without_change_of_moles_[number_of_reactions_without_mole_change++] = j;
				else if (changeOfMoles_[j] == 1.)	indices_of_reactions_with_change_of_moles_plus_one_[number_of_reactions_with_mole_change_plus_one++] = j;
				else if (changeOfMoles_[j] == -1.)	indices_of_reactions_with_change_of_moles_minus_one_[number_of_reactions_with_mole_change_minus_one++] = j;
				else								indices_of_reactions_with_change_of_moles_[number_of_reactions_with_mole_change_other++] = j;
			}
		}
	}

	void StoichiometricMap::Summary(std::ostream &fOut) const
	{
		unsigned int reactions_thermodynamically_reversible = indices_of_reactions_without_change_of_moles_.size() + indices_of_reactions_with_change_of_moles_plus_one_.size() +
			indices_of_reactions_with_change_of_moles_minus_one_.size() + indices_of_reactions_with_change_of_moles_.size();

		fOut << "----------------------------------------------------------------------------" << std::endl;
		fOut << " Reversible reactions (by thermodynamics)    " << reactions_thermodynamically_reversible << std::endl;
		fOut << "----------------------------------------------------------------------------" << std::endl;
		fOut << " Reactions without change of moles:        " << indices_of_reactions_without_change_of_moles_.size() << " (" << indices_of_reactions_without_change_of_moles_.size() / std::max(1., double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << " Reactions with change of moles (+1):      " << indices_of_reactions_with_change_of_moles_plus_one_.size() << " (" << indices_of_reactions_with_change_of_moles_minus_one_.size() / std::max(1., double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << " Reactions with change of moles (-1):      " << indices_of_reactions_with_change_of_moles_minus_one_.size() << " (" << indices_of_reactions_with_change_of_moles_plus_one_.size() / std::max(1., double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << " Reactions with change of moles (other):   " << indices_of_reactions_with_change_of_moles_.size() << " (" << indices_of_reactions_with_change_of_moles_.size() / std::max(1., double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << std::endl;
	}

	void StoichiometricMap::ProductOfConcentrations(std::vector<double>& productDirect, std::vector<double>& productReverse, const double* c)
	{
		std::fill(productDirect.begin(), productDirect.end(), 1.);
		std::fill(productReverse.begin(), productReverse.end(), 1.);

		double c1, c2, c3, csq;

		unsigned int *jD1 = lambda_jDir1_.data();
		unsigned int *jD2 = lambda_jDir2_.data();
		unsigned int *jD3 = lambda_jDir3_.data();
		unsigned int *jD4 = lambda_jDir4_.data();
		unsigned int *jD5 = lambda_jDir5_.data();
		double* vD5 = lambda_valueDir5_.data();

		unsigned int *jIE1 = lambda_jRevEq1_.data();
		unsigned int *jIE2 = lambda_jRevEq2_.data();
		unsigned int *jIE3 = lambda_jRevEq3_.data();
		unsigned int *jIE4 = lambda_jRevEq4_.data();
		unsigned int *jIE5 = lambda_jRevEq5_.data();
		double *vIE5 = lambda_valueRevEq5_.data();

		for (unsigned int i = 0; i < number_of_species_; i++)
		{
			c1 = c[i];
			c2 = c1 * c1;
			c3 = c2 * c1;
			if (lambda_numDir4_[i] != 0 || lambda_numRevEq4_[i] != 0)
				csq = std::sqrt(c1);

			for (unsigned int k = 0; k < lambda_numDir1_[i]; k++)
			{
				productDirect[*jD1] *= c1;
				jD1++;
			}
			for (unsigned int k = 0; k < lambda_numDir2_[i]; k++)
			{
				productDirect[*jD2] *= c2;
				jD2++;
			}
			for (unsigned int k = 0; k < lambda_numDir3_[i]; k++)
			{
				productDirect[*jD3] *= c3;
				jD3++;
			}
			for (unsigned int k = 0; k < lambda_numDir4_[i]; k++)
			{
				productDirect[*jD4] *= csq;
				jD4++;
			}
			for (unsigned int k = 0; k < lambda_numDir5_[i]; k++)
			{
				productDirect[*jD5] *= std::pow(c1, *vD5);
				jD5++;
				vD5++;
			}

			for (unsigned int k = 0; k < lambda_numRevEq1_[i]; k++)
			{
				productReverse[*jIE1] *= c1;
				jIE1++;
			}
			for (unsigned int k = 0; k < lambda_numRevEq2_[i]; k++)
			{
				productReverse[*jIE2] *= c2;
				jIE2++;
			}
			for (unsigned int k = 0; k < lambda_numRevEq3_[i]; k++)
			{
				productReverse[*jIE3] *= c3;
				jIE3++;
			}
			for (unsigned int k = 0; k < lambda_numRevEq4_[i]; k++)
			{
				productReverse[*jIE4] *= csq;
				jIE4++;
			}
			for (unsigned int k = 0; k < lambda_numRevEq5_[i]; k++)
			{
				productReverse[*jIE5] *= std::pow(c1, *vIE5);
				jIE5++;
				vIE5++;
			}
		}

		if (non_elementary_reactions_direct_ != 0 || non_elementary_reactions_reverse_ != 0)
			ProductOfConcentrationsForNonElementaryReactions(productDirect, productReverse, c);
	}

	void StoichiometricMap::ProductOfConcentrationsForNonElementaryReactions(std::vector<double>& productDirect, std::vector<double>& productReverse, const double* c)
	{
		const double Cstar = 1.e-8;
		const double ALFA = 1.e-5;
		const double H = 1.50*std::log(ALFA / (1. - ALFA));
		const double K = 2.00*std::log((1. - ALFA) / ALFA) / Cstar;
		const double delta = 1.e9;

		for (unsigned int j = 0; j < number_of_reactions_; j++)
		{
			if (is_non_elementary_reaction_direct_[j] == true)
			{
				productDirect[j] = 1.;

				for (unsigned int i = 0; i < non_elementary_reactions_species_indices_direct_[j].size(); i++)
				{
					const double C = c[non_elementary_reactions_species_indices_direct_[j][i]];
					const double lambda = non_elementary_reactions_orders_direct_[j][i];

					if (lambda >= 1.)
					{
						productDirect[j] *= std::pow(C, lambda);
					}
					else
					{
						const double m = (std::tanh(K*C + H) + 1.) / 2.;	// transition is for 9.e-6 < C < 1.5e-5 [kmol/m3]
						const double gamma = m*pow(C + m / delta, lambda) + (1. - m)*pow(Cstar, lambda - 1.)*C;
						productDirect[j] *= gamma;
					}
				}
			}

			if (is_non_elementary_reaction_reverse_[j] == true)
			{
				productReverse[j] = 1.;

				for (unsigned int i = 0; i < non_elementary_reactions_species_indices_reverse_[j].size(); i++)
				{
					const double C = c[non_elementary_reactions_species_indices_reverse_[j][i]];
					const double lambda = non_elementary_reactions_orders_reverse_[j][i];

					if (lambda >= 1.)
					{
						productReverse[j] *= std::pow(C, lambda);
					}
					else
					{
						const double m = (std::tanh(K*C + H) + 1.) / 2.;	// transition is for 9.e-6 < C < 1.5e-5 [kmol/m3]
						const double gamma = m*pow(C + m / delta, lambda) + (1. - m)*pow(Cstar, lambda - 1.)*C;
						productReverse[j] *= gamma;
					}
				}
			}
		}
	}

	void StoichiometricMap::FormationRatesFromReactionRates(double* R, const double* r)
	{
		unsigned int* jD1 = jDir1_.data();
		unsigned int* jD2 = jDir2_.data();
		unsigned int* jD3 = jDir3_.data();
		unsigned int* jD4 = jDir4_.data();
		unsigned int* jD5 = jDir5_.data();
		double* vD5 = valueDir5_.data();

		unsigned int* jIT1 = jRevTot1_.data();
		unsigned int* jIT2 = jRevTot2_.data();
		unsigned int* jIT3 = jRevTot3_.data();
		unsigned int* jIT4 = jRevTot4_.data();
		unsigned int* jIT5 = jRevTot5_.data();
		double* vIT5 = valueRevTot5_.data();

		for (unsigned int i = 0; i < number_of_species_; i++)
		{
			double rate = 0.;
			for (unsigned int k = 0; k < numDir1_[i]; k++)
			{
				rate -= r[*jD1];
				*jD1++;
			}
			for (unsigned int k = 0; k<numDir2_[i]; k++)
			{
				rate -= (r[*jD2] + r[*jD2]);
				*jD2++;
			}
			for (unsigned int k = 0; k<numDir3_[i]; k++)
			{
				rate -= (r[*jD3] + r[*jD3] + r[*jD3]);
				*jD3++;
			}
			for (unsigned int k = 0; k < numDir4_[i]; k++)
			{
				rate -= 0.5 * r[*jD4];
				*jD4++;
			}
			for (unsigned int k = 0; k < numDir5_[i]; k++)
			{
				rate -= (*vD5++) * r[*jD5];
				*jD5++;
			}

			for (unsigned int k = 0; k < numRevTot1_[i]; k++)
			{
				rate += r[*jIT1];
				*jIT1++;
			}
			for (unsigned int k = 0; k<numRevTot2_[i]; k++)
			{
				rate += (r[*jIT2] + r[*jIT2]);
				*jIT2++;
			}
			for (unsigned int k = 0; k<numRevTot3_[i]; k++)
			{
				rate += (r[*jIT3] + r[*jIT3] + r[*jIT3]);
				*jIT3++;
			}
			for (unsigned int k = 0; k < numRevTot4_[i]; k++)
			{
				rate += 0.5 * r[*jIT4];
				*jIT4++;
			}
			for (unsigned int k = 0; k < numRevTot5_[i]; k++)
			{
				rate += (*vIT5++) * r[*jIT5];
				*jIT5++;
			}

			R[i] = rate;
		}
	}

	// This version calculates direct and reverse reaction rates
	void StoichiometricMap::ProductionAndDestructionRatesFromReactionRatesGross(double* P, double* D, const double* rF, const double* rB)
	{
		{
			unsigned int* jD1 = jDir1_.data();
			unsigned int* jD2 = jDir2_.data();
			unsigned int* jD3 = jDir3_.data();
			unsigned int* jD4 = jDir4_.data();
			unsigned int* jD5 = jDir5_.data();
			double* vD5 = valueDir5_.data();

			unsigned int* jIT1 = jRevTot1_.data();
			unsigned int* jIT2 = jRevTot2_.data();
			unsigned int* jIT3 = jRevTot3_.data();
			unsigned int* jIT4 = jRevTot4_.data();
			unsigned int* jIT5 = jRevTot5_.data();
			double* vIT5 = valueRevTot5_.data();

			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				double destruction = 0.;
				double production = 0.;
				for (unsigned int k = 0; k<numDir1_[i]; k++)
				{
					destruction += rF[*jD1];
					production += rB[*jD1];
					*jD1++;
				}
				for (unsigned int k = 0; k<numDir2_[i]; k++)
				{
					destruction += (rF[*jD2] + rF[*jD2]);
					production += (rB[*jD2] + rB[*jD2]);
					*jD2++;
				}
				for (unsigned int k = 0; k<numDir3_[i]; k++)
				{
					destruction += (rF[*jD3] + rF[*jD3] + rF[*jD3]);
					production += (rB[*jD3] + rB[*jD3] + rB[*jD3]);
					*jD3++;
				}
				for (unsigned int k = 0; k<numDir4_[i]; k++)
				{
					destruction += .5 * rF[*jD4];
					production += .5 * rB[*jD4];
					*jD4++;
				}
				for (unsigned int k = 0; k<numDir5_[i]; k++)
				{
					destruction += (*vD5) * rF[*jD5];
					production += (*vD5) * rB[*jD5];
					*vD5++;
					*jD5++;
				}


				for (unsigned int k = 0; k<numRevTot1_[i]; k++)
				{
					production += rF[*jIT1];
					destruction += rB[*jIT1];
					*jIT1++;
				}
				for (unsigned int k = 0; k<numRevTot2_[i]; k++)
				{
					production += (rF[*jIT2] + rF[*jIT2]);
					destruction += (rB[*jIT2] + rB[*jIT2]);
					*jIT2++;
				}
				for (unsigned int k = 0; k<numRevTot3_[i]; k++)
				{
					production += (rF[*jIT3] + rF[*jIT3] + rF[*jIT3]);
					destruction += (rB[*jIT3] + rB[*jIT3] + rB[*jIT3]);
					*jIT3++;
				}
				for (unsigned int k = 0; k<numRevTot4_[i]; k++)
				{
					production += .5 * rF[*jIT4];
					destruction += .5 * rB[*jIT4];
					*jIT4++;
				}
				for (unsigned int k = 0; k<numRevTot5_[i]; k++)
				{
					production += (*vIT5) * rF[*jIT5];
					destruction += (*vIT5) * rB[*jIT5];
					*vIT5++;
					*jIT5++;
				}

				P[i] = production;
				D[i] = destruction;
			}
		}
	}

	void StoichiometricMap::ProductionAndDestructionRatesFromReactionRates(double* P, double* D, const double* r)
	{
		{
			unsigned int* jD1 = jDir1_.data();
			unsigned int* jD2 = jDir2_.data();
			unsigned int* jD3 = jDir3_.data();
			unsigned int* jD4 = jDir4_.data();
			unsigned int* jD5 = jDir5_.data();
			double* vD5 = valueDir5_.data();

			unsigned int* jIT1 = jRevTot1_.data();
			unsigned int* jIT2 = jRevTot2_.data();
			unsigned int* jIT3 = jRevTot3_.data();
			unsigned int* jIT4 = jRevTot4_.data();
			unsigned int* jIT5 = jRevTot5_.data();
			double* vIT5 = valueRevTot5_.data();

			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				double destruction = 0.;
				double production = 0.;
				for (unsigned int k = 0; k<numDir1_[i]; k++)
				{
					if (r[*jD1] > 0.)	destruction += r[*jD1];
					else                production -= r[*jD1];
					*jD1++;
				}
				for (unsigned int k = 0; k<numDir2_[i]; k++)
				{
					if (r[*jD2] > 0.)	destruction += (r[*jD2] + r[*jD2]);
					else                production -= (r[*jD2] + r[*jD2]);
					*jD2++;
				}
				for (unsigned int k = 0; k<numDir3_[i]; k++)
				{
					if (r[*jD3] > 0.)	destruction += (r[*jD3] + r[*jD3] + r[*jD3]);
					else                production -= (r[*jD3] + r[*jD3] + r[*jD3]);
					*jD3++;
				}
				for (unsigned int k = 0; k<numDir4_[i]; k++)
				{
					if (r[*jD4] > 0.)	destruction += .5 * r[*jD4];
					else                production -= .5 * r[*jD4];
					*jD4++;
				}
				for (unsigned int k = 0; k<numDir5_[i]; k++)
				{
					if (r[*jD5] > 0.)	destruction += (*vD5) * r[*jD5];
					else                production -= (*vD5) * r[*jD5];
					*vD5++;
					*jD5++;
				}


				for (unsigned int k = 0; k<numRevTot1_[i]; k++)
				{
					if (r[*jIT1] > 0.)	production += r[*jIT1];
					else				destruction -= r[*jIT1];
					*jIT1++;
				}
				for (unsigned int k = 0; k<numRevTot2_[i]; k++)
				{
					if (r[*jIT2] > 0.)	production += (r[*jIT2] + r[*jIT2]);
					else				destruction -= (r[*jIT2] + r[*jIT2]);
					*jIT2++;
				}
				for (unsigned int k = 0; k<numRevTot3_[i]; k++)
				{
					if (r[*jIT3] > 0.)	production += (r[*jIT3] + r[*jIT3] + r[*jIT3]);
					else				destruction -= (r[*jIT3] + r[*jIT3] + r[*jIT3]);
					*jIT3++;
				}
				for (unsigned int k = 0; k<numRevTot4_[i]; k++)
				{
					if (r[*jIT4] > 0.)	production += .5*r[*jIT4];
					else				destruction -= .5*r[*jIT4];
					*jIT4++;
				}
				for (unsigned int k = 0; k<numRevTot5_[i]; k++)
				{
					if (r[*jIT5] > 0.)	production += (*vIT5)*r[*jIT5];
					else				destruction -= (*vIT5)*r[*jIT5];
					*vIT5++;
					*jIT5++;
				}

				P[i] = production;
				D[i] = destruction;
			}
		}
	}

	void StoichiometricMap::ReactionEnthalpyAndEntropy(std::vector<double>& reaction_dh_over_RT, std::vector<double>& reaction_ds_over_R, const std::vector<double>& species_h_over_RT, const std::vector<double>& species_s_over_R)
	{
		unsigned int *jD1 = jDir1_.data();
		unsigned int *jD2 = jDir2_.data();
		unsigned int *jD3 = jDir3_.data();
		unsigned int *jD4 = jDir4_.data();
		unsigned int *jD5 = jDir5_.data();
		double *vD5 = valueDir5_.data();

		unsigned int *jIT1 = jRevTot1_.data();
		unsigned int *jIT2 = jRevTot2_.data();
		unsigned int *jIT3 = jRevTot3_.data();
		unsigned int *jIT4 = jRevTot4_.data();
		unsigned int *jIT5 = jRevTot5_.data();
		double *vIT5 = valueRevTot5_.data();

		std::fill(reaction_dh_over_RT.begin(), reaction_dh_over_RT.end(), 0.);
		std::fill(reaction_ds_over_R.begin(), reaction_ds_over_R.end(), 0.);

		for (unsigned int i = 0; i < number_of_species_; i++)
		{
			for (unsigned int k = 0; k<numDir1_[i]; k++)
			{
				reaction_dh_over_RT[*jD1] -= species_h_over_RT[i];
				reaction_ds_over_R[*jD1] -= species_s_over_R[i];
				jD1++;
			}
			for (unsigned int k = 0; k<numDir2_[i]; k++)
			{
				reaction_dh_over_RT[*jD2] -= (species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jD2] -= (species_s_over_R[i] + species_s_over_R[i]);
				jD2++;
			}
			for (unsigned int k = 0; k<numDir3_[i]; k++)
			{
				reaction_dh_over_RT[*jD3] -= (species_h_over_RT[i] + species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jD3] -= (species_s_over_R[i] + species_s_over_R[i] + species_s_over_R[i]);
				jD3++;
			}
			for (unsigned int k = 0; k<numDir4_[i]; k++)
			{
				reaction_dh_over_RT[*jD4] -= (0.5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jD4] -= (0.5 * species_s_over_R[i]);
				jD4++;
			}
			for (unsigned int k = 0; k<numDir5_[i]; k++)
			{
				reaction_dh_over_RT[*jD5] -= (*vD5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jD5] -= (*vD5 * species_s_over_R[i]);
				jD5++;
				vD5++;
			}

			for (unsigned int k = 0; k<numRevTot1_[i]; k++)
			{
				reaction_dh_over_RT[*jIT1] += species_h_over_RT[i];
				reaction_ds_over_R[*jIT1] += species_s_over_R[i];
				jIT1++;
			}
			for (unsigned int k = 0; k<numRevTot2_[i]; k++)
			{
				reaction_dh_over_RT[*jIT2] += (species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jIT2] += (species_s_over_R[i] + species_s_over_R[i]);
				jIT2++;
			}
			for (unsigned int k = 0; k<numRevTot3_[i]; k++)
			{
				reaction_dh_over_RT[*jIT3] += (species_h_over_RT[i] + species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jIT3] += (species_s_over_R[i] + species_s_over_R[i] + species_s_over_R[i]);
				jIT3++;
			}
			for (unsigned int k = 0; k<numRevTot4_[i]; k++)
			{
				reaction_dh_over_RT[*jIT4] += (0.5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jIT4] += (0.5 * species_s_over_R[i]);
				jIT4++;
			}
			for (unsigned int k = 0; k<numRevTot5_[i]; k++)
			{
				reaction_dh_over_RT[*jIT5] += (*vIT5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jIT5] += (*vIT5 * species_s_over_R[i]);
				jIT5++;
				vIT5++;
			}
		}
	}

	void StoichiometricMap::BuildStoichiometricMatrix()
	{
		if (isTheStoichiometricMatrixAvailable_ == false)
		{
			if (verbose_output_ == true)
				std::cout << " * Building stoichiometry..." << std::endl;

			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList_reactants;
			std::vector<T> tripletList_products;

			unsigned int estimation_of_entries_reactants = 0;
			unsigned int estimation_of_entries_products = 0;
			for (unsigned int i = 0; i<number_of_species_; i++)
			{
				estimation_of_entries_reactants += numDir1_[i];
				estimation_of_entries_reactants += numDir2_[i];
				estimation_of_entries_reactants += numDir3_[i];
				estimation_of_entries_reactants += numDir4_[i];
				estimation_of_entries_reactants += numDir5_[i];
				estimation_of_entries_products += numRevTot1_[i];
				estimation_of_entries_products += numRevTot2_[i];
				estimation_of_entries_products += numRevTot3_[i];
				estimation_of_entries_products += numRevTot4_[i];
				estimation_of_entries_products += numRevTot5_[i];
			}

			if (verbose_output_ == true)
				std::cout << "   non-zero stoichiometric coefficients: " << estimation_of_entries_reactants + estimation_of_entries_products << " /"
				<< number_of_reactions_*number_of_species_ << " ("
				<< (estimation_of_entries_reactants + estimation_of_entries_products) / std::max(1.,static_cast<double>(number_of_reactions_*number_of_species_))*100. << "%)" << std::endl;

			tripletList_reactants.reserve(estimation_of_entries_reactants);
			tripletList_products.reserve(estimation_of_entries_products);

			unsigned int *jD1 = jDir1_.data();
			unsigned int *jD2 = jDir2_.data();
			unsigned int *jD3 = jDir3_.data();
			unsigned int *jD4 = jDir4_.data();
			unsigned int *jD5 = jDir5_.data();
			double *vD5 = valueDir5_.data();

			for (unsigned int i = 0; i<number_of_species_; i++)
			{
				for (unsigned int k = 0; k<numDir1_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD1, i, 1.));
					jD1++;
				}
				for (unsigned int k = 0; k<numDir2_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD2, i, 2.));
					jD2++;
				}
				for (unsigned int k = 0; k<numDir3_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD3, i, 3.));
					jD3++;
				}
				for (unsigned int k = 0; k<numDir4_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD4, i, 0.5));
					jD4++;
				}
				for (unsigned int k = 0; k<numDir5_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD5, i, (*vD5)));
					jD5++;
					vD5++;
				}
			}

			unsigned int *jIT1 = jRevTot1_.data();
			unsigned int *jIT2 = jRevTot2_.data();
			unsigned int *jIT3 = jRevTot3_.data();
			unsigned int *jIT4 = jRevTot4_.data();
			unsigned int *jIT5 = jRevTot5_.data();
			double *vIT5 = valueRevTot5_.data();

			for (unsigned int i = 0; i<number_of_species_; i++)
			{
				for (unsigned int k = 0; k<numRevTot1_[i]; k++)
				{
					tripletList_products.push_back(T(*jIT1, i, 1.));
					jIT1++;
				}
				for (unsigned int k = 0; k<numRevTot2_[i]; k++)
				{
					tripletList_products.push_back(T(*jIT2, i, 2.));
					jIT2++;
				}
				for (unsigned int k = 0; k<numRevTot3_[i]; k++)
				{
					tripletList_products.push_back(T(*jIT3, i, 3.));
					jIT3++;
				}
				for (unsigned int k = 0; k<numRevTot4_[i]; k++)
				{
					tripletList_products.push_back(T(*jIT4, i, 0.5));
					jIT4++;
				}
				for (unsigned int k = 0; k<numRevTot5_[i]; k++)
				{
					tripletList_products.push_back(T(*jIT5, i, *vIT5));
					jIT5++;
					vIT5++;
				}
			}

			stoichiometric_matrix_reactants_.resize(number_of_reactions_, number_of_species_);
			stoichiometric_matrix_products_.resize(number_of_reactions_, number_of_species_);
			stoichiometric_matrix_reactants_.setFromTriplets(tripletList_reactants.begin(), tripletList_reactants.end());
			stoichiometric_matrix_products_.setFromTriplets(tripletList_products.begin(), tripletList_products.end());

			isTheStoichiometricMatrixAvailable_ = true;
		}
	}

	void StoichiometricMap::BuildReactionOrdersMatrix()
	{
		if (isTheReactionOrderMatrixAvailable_ == false)
		{
			if (verbose_output_ == true)
				std::cout << " * Building reaction orders..." << std::endl;

			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList_reactants;
			std::vector<T> tripletList_products;

			unsigned int estimation_of_entries_reactants = 0;
			unsigned int estimation_of_entries_products = 0;
			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				estimation_of_entries_reactants += lambda_numDir1_[i];
				estimation_of_entries_reactants += lambda_numDir2_[i];
				estimation_of_entries_reactants += lambda_numDir3_[i];
				estimation_of_entries_reactants += lambda_numDir4_[i];
				estimation_of_entries_reactants += lambda_numDir5_[i];
				estimation_of_entries_products += lambda_numRevEq1_[i];
				estimation_of_entries_products += lambda_numRevEq2_[i];
				estimation_of_entries_products += lambda_numRevEq3_[i];
				estimation_of_entries_products += lambda_numRevEq4_[i];
				estimation_of_entries_products += lambda_numRevEq5_[i];
			}

			if (verbose_output_ == true)
				std::cout << "   non-zero reaction-order coefficients: " << estimation_of_entries_reactants + estimation_of_entries_products << " /"
				<< number_of_reactions_*number_of_species_ << " ("
				<< (estimation_of_entries_reactants + estimation_of_entries_products) / std::max(1., static_cast<double>(number_of_reactions_*number_of_species_))*100. << "%)" << std::endl;

			//std::cout	<< "Mean number of species per reaction: " << (estimation_of_entries_reactants+estimation_of_entries_products)/number_of_reactions_ << std::endl;

			tripletList_reactants.reserve(estimation_of_entries_reactants);
			tripletList_products.reserve(estimation_of_entries_products);

			unsigned int *jD1 = lambda_jDir1_.data();
			unsigned int *jD2 = lambda_jDir2_.data();
			unsigned int *jD3 = lambda_jDir3_.data();
			unsigned int *jD4 = lambda_jDir4_.data();
			unsigned int *jD5 = lambda_jDir5_.data();
			double *vD5 = lambda_valueDir5_.data();

			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				for (unsigned int k = 0; k < lambda_numDir1_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD1, i, 1.));
					jD1++;
				}
				for (unsigned int k = 0; k < lambda_numDir2_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD2, i, 2.));
					jD2++;
				}
				for (unsigned int k = 0; k < lambda_numDir3_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD3, i, 3.));
					jD3++;
				}
				for (unsigned int k = 0; k < lambda_numDir4_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD4, i, 0.5));
					jD4++;
				}
				for (unsigned int k = 0; k < lambda_numDir5_[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD5, i, (*vD5)));
					jD5++;
					vD5++;
				}
			}

			unsigned int *jIE1 = lambda_jRevEq1_.data();
			unsigned int *jIE2 = lambda_jRevEq2_.data();
			unsigned int *jIE3 = lambda_jRevEq3_.data();
			unsigned int *jIE4 = lambda_jRevEq4_.data();
			unsigned int *jIE5 = lambda_jRevEq5_.data();
			double *vIE5 = lambda_valueRevEq5_.data();

			for (unsigned int i = 0; i < number_of_species_; i++)
			{
				for (unsigned int k = 0; k < lambda_numRevEq1_[i]; k++)
				{
					tripletList_products.push_back(T(*jIE1, i, 1.));
					jIE1++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq2_[i]; k++)
				{
					tripletList_products.push_back(T(*jIE2, i, 2.));
					jIE2++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq3_[i]; k++)
				{
					tripletList_products.push_back(T(*jIE3, i, 3.));
					jIE3++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq4_[i]; k++)
				{
					tripletList_products.push_back(T(*jIE4, i, 0.5));
					jIE4++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq5_[i]; k++)
				{
					tripletList_products.push_back(T(*jIE5, i, *vIE5));
					jIE5++;
					vIE5++;
				}
			}

			reactionorders_matrix_reactants_.resize(number_of_reactions_, number_of_species_);
			reactionorders_matrix_products_.resize(number_of_reactions_, number_of_species_);
			reactionorders_matrix_reactants_.setFromTriplets(tripletList_reactants.begin(), tripletList_reactants.end());
			reactionorders_matrix_products_.setFromTriplets(tripletList_products.begin(), tripletList_products.end());

			isTheReactionOrderMatrixAvailable_ = true;
		}
	}

	void StoichiometricMap::EquilibriumConstants(double* Kp, const double* exp_g_over_RT, const double Patm_over_RT)
	{
		for (unsigned int i = 0; i<number_of_reactions_; i++)
			Kp[i] = 1.;

		double c1, c2, c3, csq;

		unsigned int *jD1 = jDir1_.data();
		unsigned int *jD2 = jDir2_.data();
		unsigned int *jD3 = jDir3_.data();
		unsigned int *jD4 = jDir4_.data();
		unsigned int *jD5 = jDir5_.data();
		double *vD5 = valueDir5_.data();

		unsigned int *jIT1 = jRevTot1_.data();
		unsigned int *jIT2 = jRevTot2_.data();
		unsigned int *jIT3 = jRevTot3_.data();
		unsigned int *jIT4 = jRevTot4_.data();
		unsigned int *jIT5 = jRevTot5_.data();
		double *vIT5 = valueRevTot5_.data();

		for (unsigned int i = 0; i<number_of_species_; i++)
		{
			c1 = exp_g_over_RT[i];
			c2 = c1 * c1;
			c3 = c2 * c1;
			if (numDir4_[i] != 0 || numRevTot4_[i] != 0)
				csq = std::sqrt(c1);

			for (unsigned int k = 0; k<numDir1_[i]; k++)
			{
				Kp[*jD1] /= c1;
				jD1++;
			}
			for (unsigned int k = 0; k<numDir2_[i]; k++)
			{
				Kp[*jD2] /= c2;
				jD2++;
			}
			for (unsigned int k = 0; k<numDir3_[i]; k++)
			{
				Kp[*jD3] /= c3;
				jD3++;
			}
			for (unsigned int k = 0; k<numDir4_[i]; k++)
			{
				Kp[*jD4] /= csq;
				jD4++;
			}
			for (unsigned int k = 0; k<numDir5_[i]; k++)
			{
				Kp[*jD5] /= std::pow(c1, *vD5);
				jD5++;
				vD5++;
			}

			for (unsigned int k = 0; k<numRevTot1_[i]; k++)
			{
				Kp[*jIT1] *= c1;
				jIT1++;
			}
			for (unsigned int k = 0; k<numRevTot2_[i]; k++)
			{
				Kp[*jIT2] *= c2;
				jIT2++;
			}
			for (unsigned int k = 0; k<numRevTot3_[i]; k++)
			{
				Kp[*jIT3] *= c3;
				jIT3++;
			}
			for (unsigned int k = 0; k<numRevTot4_[i]; k++)
			{
				Kp[*jIT4] *= csq;
				jIT4++;
			}
			for (unsigned int k = 0; k<numRevTot5_[i]; k++)
			{
				Kp[*jIT5] *= std::pow(c1, *vIT5);
				jIT5++;
				vIT5++;
			}
		}

		for (unsigned int j = 0; j<indices_of_reactions_with_change_of_moles_plus_one_.size(); j++)
		{
			const unsigned int k = indices_of_reactions_with_change_of_moles_plus_one_[j] - 1;
			Kp[k] /= Patm_over_RT;
		}
		for (unsigned int j = 0; j<indices_of_reactions_with_change_of_moles_minus_one_.size(); j++)
		{
			const unsigned int k = indices_of_reactions_with_change_of_moles_minus_one_[j] - 1;
			Kp[k] *= Patm_over_RT;
		}
		for (unsigned int j = 0; j<indices_of_reactions_with_change_of_moles_.size(); j++)
		{
			const unsigned int k = indices_of_reactions_with_change_of_moles_[j] - 1;
			Kp[k] *= std::pow(Patm_over_RT, -changeOfMoles_[k]);
		}
	}

	void StoichiometricMap::RateOfProductionAnalysis(const double* r, const bool iNormalize)
	{
		BuildStoichiometricMatrix();

		if (areTheContributionOfRateOfFormationMatricesAvailable_ == false)
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList_;
			tripletList_.reserve(stoichiometric_matrix_reactants_.nonZeros() + stoichiometric_matrix_products_.nonZeros());

			// Reactants
			for (int k = 0; k<stoichiometric_matrix_reactants_.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_reactants_, k); it; ++it)
					tripletList_.push_back(T(it.row(), it.col(), -it.value()));
			}

			// Products	
			for (int k = 0; k<stoichiometric_matrix_products_.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_products_, k); it; ++it)
					tripletList_.push_back(T(it.row(), it.col(), it.value()));
			}

			Cp.resize(number_of_reactions_, number_of_species_);
			Cd.resize(number_of_reactions_, number_of_species_);
			stoichiometric_matrix_.resize(number_of_reactions_, number_of_species_);
			Cp.setFromTriplets(tripletList_.begin(), tripletList_.end());
			Cd.setFromTriplets(tripletList_.begin(), tripletList_.end());
			stoichiometric_matrix_.setFromTriplets(tripletList_.begin(), tripletList_.end());

			areTheContributionOfRateOfFormationMatricesAvailable_ = true;
		}

		// Reset
		for (int k = 0; k<Cd.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
				it.valueRef() = 0.;
		for (int k = 0; k<Cp.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
				it.valueRef() = 0.;

		// Fill (the values are not normalized)
		for (int k = 0; k<stoichiometric_matrix_.outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator itCd(Cd, k);
			Eigen::SparseMatrix<double>::InnerIterator itCp(Cp, k);
			for (Eigen::SparseMatrix<double>::InnerIterator itStoichiometry(stoichiometric_matrix_, k); itStoichiometry; ++itStoichiometry)
			{
				const double value = itStoichiometry.value() * r[itStoichiometry.row()];
				if (value >= 0.)	itCp.valueRef() = value;
				else				itCd.valueRef() = value;
				++itCp;
				++itCd;
			}
		}

		if (iNormalize == true)
		{
			// Reactants
			for (int k = 0; k<Cd.outerSize(); ++k)
			{
				double sum = 0.;
				for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
					sum += it.value();

				if (sum == 0.)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
						it.valueRef() = 0;
				}
				else
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
						it.valueRef() = it.value() / sum;
				}
			}

			// Products
			for (int k = 0; k<Cp.outerSize(); ++k)
			{
				double sum = 0.;
				for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
					sum += it.value();

				if (sum == 0.)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
						it.valueRef() = 0;
				}
				else
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
						it.valueRef() = it.value() / sum;
				}
			}
		}
	}

	void StoichiometricMap::RateOfProductionAnalysis(const double* rf, const double* rb)
	{
		BuildStoichiometricMatrix();

		if (areTheContributionOfRateOfFormationMatricesAvailable_ == false)
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList_;
			tripletList_.reserve(stoichiometric_matrix_reactants_.nonZeros() + stoichiometric_matrix_products_.nonZeros());

			// Reactants
			for (int k = 0; k<stoichiometric_matrix_reactants_.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_reactants_, k); it; ++it)
					tripletList_.push_back(T(it.row(), it.col(), -it.value()));
			}

			// Products	
			for (int k = 0; k<stoichiometric_matrix_products_.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_products_, k); it; ++it)
					tripletList_.push_back(T(it.row(), it.col(), it.value()));
			}

			Cp.resize(number_of_reactions_, number_of_species_);
			Cd.resize(number_of_reactions_, number_of_species_);
			stoichiometric_matrix_.resize(number_of_reactions_, number_of_species_);
			Cp.setFromTriplets(tripletList_.begin(), tripletList_.end());
			Cd.setFromTriplets(tripletList_.begin(), tripletList_.end());
			stoichiometric_matrix_.setFromTriplets(tripletList_.begin(), tripletList_.end());

			areTheContributionOfRateOfFormationMatricesAvailable_ = true;
		}

		// Reset
		for (int k = 0; k<Cd.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
				it.valueRef() = 0.;
		for (int k = 0; k<Cp.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
				it.valueRef() = 0.;

		// Fill (the values are not normalized)
		for (int k = 0; k<stoichiometric_matrix_.outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator itCd(Cd, k);
			Eigen::SparseMatrix<double>::InnerIterator itCp(Cp, k);
			for (Eigen::SparseMatrix<double>::InnerIterator itStoichiometry(stoichiometric_matrix_, k); itStoichiometry; ++itStoichiometry)
			{
				// Production
				{
					const double value = itStoichiometry.value() * rf[itStoichiometry.row()];
					if (value >= 0.)	itCp.valueRef() = value;
					else				itCd.valueRef() = value;
				}

				// Destruction
				{
					const double value = itStoichiometry.value() * rb[itStoichiometry.row()];
					if (value <= 0.)	itCp.valueRef() += -value;
					else				itCd.valueRef() += -value;
				}

				++itCp;
				++itCd;
			}
		}
	}

	void StoichiometricMap::WriteRateOfProductionAnalysis(std::ostream& fout)
	{
		for (int k = 0; k<Cd.outerSize(); ++k)
		{
			// Calculates the sum of destruction rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			fout << sum << std::endl;

			// Writes the reactions involved
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
				fout << it.row() << " ";
			fout << std::endl;

			// Writes the coefficients (not normalized)
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
				fout << it.value() << " ";
			fout << std::endl;
		}

		for (int k = 0; k<Cp.outerSize(); ++k)
		{
			// Calculates the sum of production rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			fout << sum << std::endl;

			// Writes the reactions involved
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
				fout << it.row() << " ";
			fout << std::endl;

			// Writes the coefficients (not normalized)
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
				fout << it.value() << " ";
			fout << std::endl;
		}
	}

	void StoichiometricMap::WriteRateOfProductionAnalysis(ROPA_Data& ropa)
	{
		ropa.destruction_rates.resize(Cd.outerSize());
		ropa.production_rates.resize(Cp.outerSize());

		ropa.destruction_coefficients.resize(Cd.outerSize());
		ropa.destruction_reaction_indices.resize(Cd.outerSize());
		ropa.production_coefficients.resize(Cp.outerSize());
		ropa.production_reaction_indices.resize(Cp.outerSize());

		for (int k = 0; k<Cd.outerSize(); ++k)
		{
			// Calculates the sum of destruction rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			ropa.destruction_rates[k] = sum;

			// Writes the reactions involved
			ropa.destruction_reaction_indices[k].resize(count);
			unsigned int j1 = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
				ropa.destruction_reaction_indices[k][j1++] = it.row();

			// Writes the coefficients (not normalized)
			ropa.destruction_coefficients[k].resize(count);
			unsigned int j2 = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd, k); it; ++it)
				ropa.destruction_coefficients[k][j2++] = it.value();
		}

		for (int k = 0; k<Cp.outerSize(); ++k)
		{
			// Calculates the sum of production rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			ropa.production_rates[k] = sum;

			// Writes the reactions involved
			ropa.production_reaction_indices[k].resize(count);
			unsigned int j1 = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
				ropa.production_reaction_indices[k][j1++] = it.row();

			// Writes the coefficients (not normalized)
			ropa.production_coefficients[k].resize(count);
			unsigned int j2 = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp, k); it; ++it)
				ropa.production_coefficients[k][j2++] = it.value();
		}
	}

	void StoichiometricMap::GetSumOfStoichiometricCoefficientsOfReactants(Eigen::VectorXd& sum_nu) const
	{
		sum_nu.resize(number_of_reactions_);
		sum_nu.setZero();
		for (int k = 0; k<stoichiometric_matrix_reactants_.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_reactants_, k); it; ++it)
				sum_nu(it.row()) += it.value();
		}
	}

	void StoichiometricMap::GetSumOfStoichiometricCoefficientsOfProducts(Eigen::VectorXd& sum_nu) const
	{
		sum_nu.resize(number_of_reactions_);
		sum_nu.setZero();
		for (int k = 0; k<stoichiometric_matrix_products_.outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_products_, k); it; ++it)
				sum_nu(it.row()) += it.value();
		}
	}
}


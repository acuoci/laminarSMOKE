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

namespace OpenSMOKE
{

	AtomicCompositionTable::AtomicCompositionTable(const AtomicCompositionTable& orig)
	{
		 element_weights_list_ = orig.element_weights_list_;
		 element_coefficients_list_ = orig.element_coefficients_list_;
		 element_names_list_.resize(orig.element_weights_list_.size());
		 for (unsigned int i=0;element_names_list_.size();i++)
			 element_names_list_[i] = orig.element_names_list_[i];
	}

	void AtomicCompositionTable::CopyFrom(const AtomicCompositionTable& orig)
	{
		 element_weights_list_.resize(orig.element_weights_list_.size());
		 element_coefficients_list_.resize(orig.element_coefficients_list_.rows(), orig.element_coefficients_list_.cols());
		 element_names_list_.resize(orig.element_names_list_.size());
		 
		 element_weights_list_ = orig.element_weights_list_;
		 element_coefficients_list_ = orig.element_coefficients_list_;
		 for (unsigned int i=0;i<element_names_list_.size();i++)
			 element_names_list_[i] = orig.element_names_list_[i];
	}


	template<typename SpeciesMap>
	void AtomicCompositionTable::CreateAtomicElementLists(const SpeciesMap& species)
	{
		// Recover atomic elements
		for(unsigned int i=0;i<species.size();i++)
		{
			OpenSMOKE::AtomicComposition atomic;
			species[i].AtomicComposition(&atomic);

			for(unsigned int k=0;k<atomic.element_names().size();k++)
			{
				bool iFound = false;
				for(unsigned int j=0;j<element_names_list_.size();j++)
					if ( atomic.element_names()[k] == element_names_list_[j] )
					{
						iFound = true;
						break;
					}

				if (iFound == false)
					element_names_list_.push_back(atomic.element_names()[k]);
			}
		}

		// Reorder
		{
			unsigned current_index = 0;
			std::vector<std::string> target_elements = { "C", "H", "O", "N" };

			for (unsigned int k = 0; k < target_elements.size(); k++)
			{
				for (unsigned int j = 0; j < element_names_list_.size(); j++)
					if (element_names_list_[j] == target_elements[k])
					{
						std::iter_swap(element_names_list_.begin() + current_index, element_names_list_.begin() + j);
						current_index++;
						break;
					}
			}
		}

		// Recover atomic weights
		element_weights_list_.resize(element_names_list_.size());
		element_weights_list_.setZero();
		for(unsigned int j=0;j<element_names_list_.size();j++)
			element_weights_list_(j) = OpenSMOKE::AtomicWeights[element_names_list_[j]];

		// Recover atomic elements
		element_coefficients_list_.resize(species.size(), element_names_list_.size());
		element_coefficients_list_.setZero();
		for(unsigned int i=0;i<species.size();i++)
		{
			OpenSMOKE::AtomicComposition atomic;
			species[i].AtomicComposition(&atomic);

			for(unsigned int k=0;k<atomic.element_names().size();k++)
			{
				bool iFound = false;
				for(unsigned int j=0;j<element_names_list_.size();j++)
					if ( atomic.element_names()[k] == element_names_list_[j] )
					{
						element_coefficients_list_(i,j) = atomic.element_coefficients()[k];
						break;
					}
			}
		}
	}

	template<typename Reaction>
	unsigned int AtomicCompositionTable::CheckStoichiometry(std::ostream& fOut, const Reaction& reaction, const double epsilon)
	{
		Eigen::VectorXd reactant_side(element_weights_list_.size());
		Eigen::VectorXd product_side(element_weights_list_.size());
		reactant_side.setZero();
		product_side.setZero();
		for(unsigned int i=0;i<reaction.reactant_nu_indices().size();i++)
			for(int j=0;j<reactant_side.size();j++)
				reactant_side(j) += element_coefficients_list_(reaction.reactant_nu_indices()[i],j)*reaction.reactant_nu()[i];
		for(unsigned int i=0;i<reaction.product_nu_indices().size();i++)
			for(int j=0;j<product_side.size();j++)
				product_side(j) += element_coefficients_list_(reaction.product_nu_indices()[i],j)*reaction.product_nu()[i];

		Eigen::VectorXd relative_errors(element_weights_list_.size());
		for(int i=0;i<reactant_side.size();i++)
			relative_errors(i) = std::fabs(reactant_side(i)-product_side(i))/std::max((reactant_side(i)+product_side(i))*0.5, 1.e-32);
		
		for(int i=0;i<reactant_side.size();i++)
			if (relative_errors(i)>epsilon)
			{
				fOut << "Error in reaction stoichiometry" << std::endl;
				fOut << "Atom      Reactants       Products        Rel. error(%)" << std::endl;
				for(int k=0;k<reactant_side.size();k++)
				{
					if (relative_errors(k)*100. > epsilon)
					{
						fOut << std::left << std::setw(10) << element_names_list_[k];
						fOut << std::left << std::setw(16) << reactant_side(k);
						fOut << std::left << std::setw(16) << product_side(k);
						fOut << std::left << std::setw(16) << relative_errors(k)*100.;
						fOut << std::endl;
					}
				}
				return 1;
			}

		double max_relative_errors = relative_errors.maxCoeff();
		if (max_relative_errors>1.e-12)
		{
			fOut << "Warning: the reaction is not perfectly balanced!" << std::endl;
			fOut << "Atom      Reactants       Products        Rel. error(%)" << std::endl;
			for(int k=0;k<reactant_side.size();k++)
			{
				if (relative_errors(k)*100. > 1e-10)
				{
					fOut << std::left << std::setw(10) << element_names_list_[k];
					fOut << std::left << std::setw(16) << reactant_side(k);
					fOut << std::left << std::setw(16) << product_side(k);
					fOut << std::left << std::setw(16) << relative_errors(k)*100.;
					fOut << std::endl;
				}
			}

			return 2;
		}

		return 0;
	}

	template<typename Reaction>
	unsigned int AtomicCompositionTable::CorrectStoichiometry(Reaction& reaction)
	{
		Eigen::VectorXd reactant_side(element_weights_list_.size());
		Eigen::VectorXd product_side(element_weights_list_.size());
		reactant_side.setZero();
		product_side.setZero();
		for (unsigned int i = 0; i<reaction.reactant_nu_indices().size(); i++)
		for (int j = 0; j<reactant_side.size(); j++)
			reactant_side(j) += element_coefficients_list_(reaction.reactant_nu_indices()[i], j)*reaction.reactant_nu()[i];
		for (unsigned int i = 0; i<reaction.product_nu_indices().size(); i++)
		for (int j = 0; j<product_side.size(); j++)
			product_side(j) += element_coefficients_list_(reaction.product_nu_indices()[i], j)*reaction.product_nu()[i];

		// Propose correction
		{
			// List of available elements
			std::vector<unsigned int>	indices_available_elements;
			for (int k = 0; k<reactant_side.size(); k++)
			if (reactant_side(k) > 0.)
				indices_available_elements.push_back(k);

			unsigned int m = indices_available_elements.size();
			unsigned int n = reaction.product_nu_indices().size();

			if (m <= n)
			{
				Eigen::VectorXd	beta(m);
				for (unsigned int k = 0; k < m; k++)
				{
					unsigned int j = indices_available_elements[k];
					beta(k) = reactant_side(j) - product_side(j);
				}

				Eigen::MatrixXd nu(m, n);
				for (unsigned int i = 0; i<reaction.product_nu_indices().size(); i++)
				for (unsigned int k = 0; k < m; k++)
				{
					unsigned int j = indices_available_elements[k];
					nu(k, i) = element_coefficients_list_(reaction.product_nu_indices()[i], j);
				}

				Eigen::VectorXd delta(n);
				if (m == n)
					delta = nu.fullPivLu().solve(beta);
				else
					delta = nu.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(beta);

				for (unsigned int i = 0; i < reaction.product_nu_indices().size(); i++)
				{
					const double new_coefficient = reaction.product_nu()[i] + delta(i);
					reaction.set_product_nu(i, new_coefficient);
				}

				return 1;
			}
			else
			{
				return 0;
			}
		}
	}

}
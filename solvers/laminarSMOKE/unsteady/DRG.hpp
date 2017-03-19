#include <queue>

namespace OpenSMOKE
{
	DRG::DRG(OpenSMOKE::ThermodynamicsMap_CHEMKIN* thermodynamicsMapXML,
			 OpenSMOKE::KineticsMap_CHEMKIN* kineticsMapXML) :

			 thermodynamicsMapXML_(*thermodynamicsMapXML),
			 kineticsMapXML_(*kineticsMapXML)
	{
		epsilon_ = 1.e-2;

		NS_ = thermodynamicsMapXML_.NumberOfSpecies();
		NR_ = kineticsMapXML_.NumberOfReactions();

		ChangeDimensions(NR_, &rNet_, true);
		r_.resize(NS_, NS_);
		

		important_species_.resize(NS_);
                important_reactions_.resize(NR_);

		number_important_species_ = NS_;
		number_unimportant_reactions_ = 0;

		// Build a full matrix of net stoichiometric coefficients nu = nuB - nuF
		Eigen::MatrixXd nu_(NR_, NS_);						
		{
			// Be careful: eigen vectors and matrices are 0-index based
			// Be careful: if the kinetic scheme is large, this matrix, since it is full, can be very memory expensive
			//             Example: 10^3 species, 10^4 reactions = size of the matrix 10^7 elements!
			//             This is the reason why we store stoichiometric matrices in sparse format.
			//             Of course te advantage of having a full matrix, is that you access the elements directly, without
			//             using iterators and pointers, as reported above
			nu_.setZero();

			// Loop over all the reactions (product side)
			for (int k = 0; k < kineticsMapXML_.stoichiometry().stoichiometric_matrix_products().outerSize(); ++k)
			{
				// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
				for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_products(), k); it; ++it)
				{
					nu_(it.row(), it.col()) += it.value();
				}
			}

			// Loop over all the reactions (product side)
			for (int k = 0; k < kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants().outerSize(); ++k)
			{
				// Loop over all the non-zero stoichiometric coefficients (product side) of reaction k
				for (Eigen::SparseMatrix<double>::InnerIterator it(kineticsMapXML_.stoichiometry().stoichiometric_matrix_reactants(), k); it; ++it)
				{
					nu_(it.row(), it.col()) -= it.value();
				}
			}
		}

		// Build the delta matrix (dense matrix) used by the DRG method
		Eigen::MatrixXd delta_(NR_, NS_);					
		{
			for (unsigned int i = 0; i < NR_; i++)
			{
				for (unsigned int j = 0; j < NS_; j++)
					delta_(i, j) = (nu_(i, j) == 0) ? 0 : 1;
			}
		}

		// Sparse Algebra
		{			
			// Sparse net stoichiometric matrix
			{
				nu_sparse_.resize(NS_,NR_);

				typedef Eigen::Triplet<double> T;
				std::vector<T> tripletList;
				tripletList.reserve(NR_*4);
				for (unsigned int i = 0; i < NR_; i++)
					for (unsigned int j = 0; j < NS_; j++)
					{
  						if ( nu_(i,j) != 0.)
  						tripletList.push_back(T(j,i,nu_(i,j)));
					}
		
				nu_sparse_.setFromTriplets(tripletList.begin(), tripletList.end());
			}

			// Sparse delta matrix, 1 means the species is involved in the reaction, (NR x NS)
			delta_sparse_.resize(NS_,NR_);
			{
				typedef Eigen::Triplet<double> T;
				std::vector<T> tripletList;
				tripletList.reserve(NR_*4);
				for (unsigned int i = 0; i < NR_; i++)
					for (unsigned int j = 0; j < NS_; j++)
					{
  						if ( delta_(i,j) != 0.)
  						tripletList.push_back(T(j,i,1.));
					}
		
				delta_sparse_.setFromTriplets(tripletList.begin(), tripletList.end());
			}

			// Nu times delta
			{
				nu_times_delta_.resize(NS_);

				for (unsigned int k = 0; k < NS_; k++)
				{
					nu_times_delta_[k].resize(NS_, NR_);

					typedef Eigen::Triplet<double> T;
					std::vector<T> tripletList;
					tripletList.reserve(NR_*4);
					for (unsigned int i = 0; i < NR_; i++)
						for (unsigned int j = 0; j < NS_; j++)
						{
							const double prod = nu_(i,k) * delta_(i,j);
  							if ( prod != 0.)
  								tripletList.push_back(T(j,i, prod));
						}
		
					nu_times_delta_[k].setFromTriplets(tripletList.begin(), tripletList.end());
				}
			}
		}
	}

	void DRG::SetKeySpecies(const std::vector<std::string> names_key_species)
	{
		index_key_species_.resize(names_key_species.size());
		for (unsigned int i = 0; i < names_key_species.size(); i++)
			index_key_species_[i] = thermodynamicsMapXML_.IndexOfSpecies(names_key_species[i]) - 1;
	}

	void DRG::SetKeySpecies(const std::vector<unsigned int> key_species)
	{
		index_key_species_ = key_species;
	}

	void DRG::SetEpsilon(const double epsilon)
	{
		epsilon_ = epsilon;
	}

	void DRG::Analysis(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c)
	{
		PairWiseErrorMatrix(T, P_Pa, c);
		ParsePairWiseErrorMatrix();
	}

	void DRG::PairWiseErrorMatrix(const double T, const double P_Pa, const OpenSMOKE::OpenSMOKEVectorDouble& c)
	{
		// Now we know T, P and composition. 
		// We have to pass those data to the thermodynamic and kinetic maps
		kineticsMapXML_.SetTemperature(T);
		kineticsMapXML_.SetPressure(P_Pa);
		thermodynamicsMapXML_.SetTemperature(T);
		thermodynamicsMapXML_.SetPressure(P_Pa);

		// Now we can calculate (internally) the reaction rates concentrations are needed
		kineticsMapXML_.ReactionRates(c.GetHandle());
		kineticsMapXML_.GiveMeReactionRates(rNet_.GetHandle());	// [kmol/m3/s]

		// Calculate the pair-wise error matrix
		r_.setZero();
		{
			Eigen::VectorXd numerator_(NS_);
			Eigen::VectorXd denominator_(NS_);

			// Denominator			
			denominator_.setConstant(0.);		
			for (int k=0; k<nu_sparse_.outerSize(); ++k)
  			{
				const double rnet = rNet_[k+1];
				for (Eigen::SparseMatrix<double>::InnerIterator it(nu_sparse_,k); it; ++it)
  				{
    					denominator_[it.row()] += std::fabs(it.value() * rnet);
  				}
			}

			// Numerator
			for (int i = 0; i < NS_; i++)
			{
				numerator_.setConstant(0.);
				for (int k=0; k<nu_times_delta_[i].outerSize(); ++k)
  				{
					const double rnet = rNet_[k+1];
					for (Eigen::SparseMatrix<double>::InnerIterator it(nu_times_delta_[i],k); it; ++it)
  					{
    						numerator_[it.row()] += std::fabs(it.value() * rnet);
  					}
				}

				for (int j = 0; j < NS_; j++)
					r_(i,j) = numerator_[j]/(1.e-64+denominator_[i]);
			}
			
		}
        }
        
	void DRG::ParsePairWiseErrorMatrix()
	{
		// Reset important species and important reactions
		important_species_.assign(NS_,false);
		important_reactions_.assign(NR_,true);
     
		// Initialize the queue with key-species
		std::queue <int> Q;
		for ( int i=0; i<index_key_species_.size(); i++)
		{
			Q.push(index_key_species_[i]);
			important_species_[index_key_species_[i]]=true;
		}
           
		// DFS with rAB
		while (!Q.empty())
		{
			for ( int k=0; k<NS_; k++)
			{
				if (important_species_[k] == false)                    
				{
					if (r_(Q.front(), k) > epsilon_)
					{                            
						important_species_[k] = true;
						Q.push(k);
					}
				}
			}
			Q.pop();
		}	
                
		// Important reactions
		for (int k=0; k<delta_sparse_.outerSize(); ++k)
  		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(delta_sparse_,k); it; ++it)
  			{
				if (important_species_[it.row()] == false)
				{
					important_reactions_[k] = false;
				}
  			}
		}

		// Count important species and reactions
		number_important_species_ = std::count (important_species_.begin(), important_species_.end(), true);
		number_unimportant_reactions_ = std::count (important_reactions_.begin(), important_reactions_.end(), false);

		// Vector containing the indices of important species (zero-based)
		{
			indices_important_species_.resize(number_important_species_);
			unsigned int count = 0;
			for(unsigned int k = 0; k < NS_; k++)
				if (important_species_[k] == true)
				{
					indices_important_species_[count] = k;
					count++;
				}
                }

		// Vector containing the indices of unimportant species (zero-based)
		{
			indices_unimportant_reactions_.resize(number_unimportant_reactions_);
			unsigned int count = 0;
			for(unsigned int k = 0; k < NR_; k++)
				if (important_reactions_[k] == false)
				{
					indices_unimportant_reactions_[count] = k;
					count++;
				}
                }
			
	}
}

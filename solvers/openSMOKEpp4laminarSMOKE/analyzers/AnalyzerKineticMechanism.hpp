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

#include <boost/config.hpp>
#include <vector>
#include <iostream>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>

namespace OpenSMOKE
{
	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	AnalyzerKineticMechanism<Kinetics_PreProcessor, Kinetics_Map>::AnalyzerKineticMechanism(Kinetics_PreProcessor& kinetics_preprocessor, Kinetics_Map& kinetics_map) :
		kinetics_preprocessor_(kinetics_preprocessor), kinetics_map_(kinetics_map)
	{
	}

	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	bool AnalyzerKineticMechanism<Kinetics_PreProcessor, Kinetics_Map>::WriteReactionTablesOnASCIIFile(const std::string& file_name, const std::vector<double> list_of_temperatures) const
	{
		std::cout << " * Writing the reaction tables on ascii file" << std::endl;

		const unsigned int nc = (unsigned int)(kinetics_map_.NamesOfSpecies().size());
		
		OpenSMOKE::OpenSMOKEVectorDouble c_bath(nc);
		c_bath = 1./double(nc);

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		fOutput << std::endl;
		fOutput << " ================================================================================================================================" << std::endl;
		fOutput << "   Kinetic data summary @  298.15K"  << std::endl;
		fOutput << " ================================================================================================================================" << std::endl;
		fOutput << "        Reaction                                                                        DG            DH            DS               " << std::endl;          
		fOutput << "                                                                                        [kcal/mol]    [kcal/mol]    [cal/mol/K]      " << std::endl;        
		fOutput << " ================================================================================================================================" << std::endl;

		for (unsigned int i=0;i<kinetics_preprocessor_.reactions().size();i++)
		{
			std::string reaction_string;
			kinetics_preprocessor_.reactions()[i].GetReactionString(kinetics_map_.NamesOfSpecies(), reaction_string);
			boost::erase_all(reaction_string, " ");
			if (reaction_string.size()<80)
				fOutput << " " << std::setw(7) << std::left << i+1 << std::setw(80) << reaction_string;
			else
				fOutput << " " << std::setw(7) << std::left << i+1 << std::setw(80) << reaction_string.substr(0,76) << "...";
			kinetics_map_.WriteKineticData(fOutput, i+1);
			fOutput << std::endl;
		}

		fOutput << " ================================================================================================================================" << std::endl;
		fOutput << std::endl;
		fOutput << std::endl;

		for (unsigned int i=0;i<kinetics_preprocessor_.reactions().size();i++)
		{
			kinetics_preprocessor_.reactions()[i].WriteSummary(fOutput, kinetics_map_.NamesOfSpecies(), i+1);
			const double conversion_forward = kinetics_preprocessor_.reactions()[i].GetForwardConversionFactor();
			const double conversion_backward = kinetics_preprocessor_.reactions()[i].GetBackwardConversionFactor();
			kinetics_map_.WriteKineticData(fOutput, i+1, c_bath.GetHandle(), list_of_temperatures, conversion_forward, conversion_backward);
		}

		fOutput.close();

		return true;

	}

	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	bool AnalyzerKineticMechanism<Kinetics_PreProcessor, Kinetics_Map>::WriteCollisionRatesOnASCIIFile(TransportPropertiesMap_CHEMKIN& transport, const std::string& file_name, const std::vector<double> list_of_temperatures) const
	{
		// Reconstruct reaction names
		std::vector<std::string> reaction_names(kinetics_preprocessor_.reactions().size());
		for (unsigned int i = 0; i<kinetics_preprocessor_.reactions().size(); i++)
		{
			std::string reaction_string;
			kinetics_preprocessor_.reactions()[i].GetReactionString(kinetics_map_.NamesOfSpecies(), reaction_names[i]);
			boost::erase_all(reaction_names[i], " ");
		}

		// Analyze
		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);
		kinetics_map_.WriteCollisionRateConstantForBimolecularReactions(transport, fOutput, reaction_names, list_of_temperatures);
		fOutput.close();

		return true;
	}

	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	bool AnalyzerKineticMechanism<Kinetics_PreProcessor, Kinetics_Map>::WriteFittedInverseKineticConstantsOnASCIIFile(const std::string& file_name) const
	{
		std::cout << " * Fitting the reverse reactions..." << std::endl;

		const unsigned int nc = (unsigned int)(kinetics_map_.NamesOfSpecies().size());
		
		OpenSMOKE::OpenSMOKEVectorDouble x_bath(nc);
		x_bath = 1./double(nc);

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		// Details (every reaction, even if not reversible)
		{
			fOutput << "---------------------------------------------------------------------------------------" << std::endl;
			fOutput << "                                  CHEMICAL REACTIONS                                   " << std::endl;
			fOutput << std::endl;
			fOutput << "                          Units: [mol, cm3, s] and [cal/mol]                           " << std::endl;
			fOutput << "---------------------------------------------------------------------------------------" << std::endl;
			fOutput << std::endl;
			fOutput << std::endl;

			Eigen::MatrixXd fittedKineticParameters;
			kinetics_map_.FittedReverseKineticConstants(x_bath.GetHandle(), 2, fittedKineticParameters, false);

			for (unsigned int k = 1; k <= kinetics_map_.NumberOfReactions(); k++)
			{
				//unsigned int i = kinetics_map_.indices_of_reversible_reactions()[k + 1];

				fOutput << std::setw(7) << std::right << k;
				fOutput << ". ";
				kinetics_preprocessor_.reactions()[k - 1].WriteShortSummary(fOutput, kinetics_map_.NamesOfSpecies(), fittedKineticParameters.col(k-1));
				fOutput << std::endl;
				fOutput << std::endl;
			}
		}

		// Fitting with 2 parameters
		{
			fOutput << " ================================================================================================================================" << std::endl;
			fOutput << "    Reaction        A           Beta           E         Reaction                                                                   " << std::endl;
			fOutput << "     index     [kmol,m3,s]                 [cal/mol]                                                                                " << std::endl;
			fOutput << " ================================================================================================================================" << std::endl;

			std::cout << "   2 parameters fitting..." << std::endl;

			Eigen::MatrixXd fittedKineticParameters;
			kinetics_map_.FittedReverseKineticConstants(x_bath.GetHandle(), 2, fittedKineticParameters, true);

			std::vector<size_t> indices;
			{
				std::vector<double> E(fittedKineticParameters.cols());
				for (unsigned int i = 0; i < E.size();i++)
					E[i] = fittedKineticParameters(1,i);

				indices = SortAndTrackIndicesIncreasing(E);
			}

			for (unsigned int k = 0; k < indices.size(); k++)
			{
				unsigned int i = kinetics_map_.IndicesOfReversibleReactions()[indices[k]];

				fOutput << " " << std::setw(7) << std::left << i;
				kinetics_map_.FittedReverseKineticConstants(i, fOutput, fittedKineticParameters);

				std::string reaction_string;
				kinetics_preprocessor_.reactions()[i - 1].GetReactionString(kinetics_map_.NamesOfSpecies(), reaction_string);
				boost::erase_all(reaction_string, " ");

				if (reaction_string.size() < 80)
					fOutput << std::left << std::setw(80) << reaction_string;
				else
					fOutput << std::left << std::setw(80) << reaction_string.substr(0, 76) << "...";
				fOutput << std::endl;
			}

			fOutput << std::endl;
		}

		// Fitting with 3 parameters
		{
			fOutput << " ================================================================================================================================" << std::endl;
			fOutput << "    Reaction        A           Beta           E         Reaction                                                                   " << std::endl;
			fOutput << "     index     [kmol,m3,s]                 [cal/mol]                                                                                " << std::endl;
			fOutput << " ================================================================================================================================" << std::endl;

			std::cout << "   3 parameters fitting..." << std::endl;

			Eigen::MatrixXd fittedKineticParameters;
			kinetics_map_.FittedReverseKineticConstants(x_bath.GetHandle(), 3, fittedKineticParameters, true);

			std::vector<size_t> indices;
			{
				std::vector<double> E(fittedKineticParameters.cols());
				for (unsigned int i = 0; i < E.size(); i++)
					E[i] = fittedKineticParameters(1, i);

				indices = SortAndTrackIndicesIncreasing(E);
			}

			for (unsigned int k = 0; k < indices.size(); k++)
			{
				unsigned int i = kinetics_map_.IndicesOfReversibleReactions()[indices[k]];

				fOutput << " " << std::setw(7) << std::left << i;
				kinetics_map_.FittedReverseKineticConstants(i, fOutput, fittedKineticParameters);

				std::string reaction_string;
				kinetics_preprocessor_.reactions()[i - 1].GetReactionString(kinetics_map_.NamesOfSpecies(), reaction_string);
				boost::erase_all(reaction_string, " ");

				if (reaction_string.size() < 80)
					fOutput << std::left << std::setw(80) << reaction_string;
				else
					fOutput << std::left << std::setw(80) << reaction_string.substr(0, 76) << "...";
				fOutput << std::endl;
			}
		}

		fOutput.close();
		return true;
	}

	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	bool AnalyzerKineticMechanism<Kinetics_PreProcessor, Kinetics_Map>::WriteFittedChebyshevOnASCIIFile(const std::string& file_name) const
	{
		std::cout << " * Fitting the Chebishev reactions..." << std::endl;

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		fOutput.setf(std::ios::scientific);

		fOutput << "---------------------------------------------------------------------------------------" << std::endl;
		fOutput << "                             CHEBYSHEV CHEMICAL REACTIONS                              " << std::endl;
		fOutput << std::endl;
		fOutput << "                       Units: [atm], [mol, cm3, s] and [cal/mol]                       " << std::endl;
		fOutput << "---------------------------------------------------------------------------------------" << std::endl;
		fOutput << std::endl;
		fOutput << std::endl;

		for (unsigned int j = 0; j < kinetics_map_.number_of_chebyshev_reactions(); j++)
		{
			unsigned int index_reaction = kinetics_map_.indices_of_chebyshev_reactions()[j];
			std::string reaction_string;
			kinetics_preprocessor_.reactions()[index_reaction-1].GetReactionString(kinetics_map_.NamesOfSpecies(), reaction_string);
			boost::erase_all(reaction_string, " ");

			fOutput << "---------------------------------------------------------------------------------------" << std::endl;
			if (reaction_string.size() < 80)	fOutput << std::left << std::setw(80) << reaction_string;
			else                                fOutput << std::left << std::setw(80) << reaction_string.substr(0, 76) << "...";
			fOutput << std::endl;
			fOutput << "---------------------------------------------------------------------------------------" << std::endl;

			kinetics_map_.chebyshev_reactions(j).FitttingArrheniusLaw(fOutput);
		}
		
		fOutput.close();
		return true;
	}

	template<typename Kinetics_PreProcessor, typename Kinetics_Map>
	bool AnalyzerKineticMechanism<Kinetics_PreProcessor, Kinetics_Map>::SparsityPatternAnalysis(const std::string& file_name) const
	{
		const unsigned int nc = (unsigned int)(kinetics_map_.NamesOfSpecies().size());
		const unsigned int nr = (unsigned int)(kinetics_map_.NumberOfReactions());

		std::ofstream fOutput;
		fOutput.open(file_name.c_str(), std::ios::out);
		
		fOutput << "% Jacobian matrix (weak version)" << std::endl;
		fOutput << nc << " " << nc << std::endl;
		fOutput << kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix()->nonZeros() << std::endl;
		for (int k = 0; k < kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix(), k); it; ++it)
			{
				fOutput << it.row() + 1 << " " << it.col() + 1 << std::endl;
			}
		}
		
		
		fOutput << "% Derivative of reaction rates matrix (forward) (weak version)" << std::endl;
		fOutput << nr << " " << nc << std::endl;
		fOutput << kinetics_map_.jacobian_sparsity_pattern_map()->drf_over_domega()->nonZeros() << std::endl;
		for (int k = 0; k < kinetics_map_.jacobian_sparsity_pattern_map()->drf_over_domega()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*kinetics_map_.jacobian_sparsity_pattern_map()->drf_over_domega(), k); it; ++it)
			{
				fOutput << it.row() + 1 << " " << it.col() + 1 << std::endl;
			}
		}
		
		fOutput << "% Derivative of reaction rates matrix (backward) (weak version)" << std::endl;
		fOutput << nr << " " << nc << std::endl;
		fOutput << kinetics_map_.jacobian_sparsity_pattern_map()->drb_over_domega()->nonZeros() << std::endl;
		for (int k = 0; k < kinetics_map_.jacobian_sparsity_pattern_map()->drb_over_domega()->outerSize(); ++k)
		{
			for (Eigen::SparseMatrix<double>::InnerIterator it(*kinetics_map_.jacobian_sparsity_pattern_map()->drb_over_domega(), k); it; ++it)
			{
				fOutput << it.row() + 1 << " " << it.col() + 1 << std::endl;
			}
		}
		
		
		// Graph analysis
		{
			typedef boost::adjacency_list<
				boost::vecS, boost::vecS, boost::undirectedS,
				boost::property<boost::vertex_color_t, boost::default_color_type,
				boost::property<boost::vertex_degree_t, int> > > Graph;

			typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
			typedef boost::graph_traits<Graph>::vertices_size_type size_type;

			typedef std::pair<std::size_t, std::size_t> Pair;

			unsigned int count = 0;
			for (int k = 0; k < kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix()->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix(), k); it; ++it)
				{
						count++;
				}
			}
			
			Pair* edges = new Pair[count];
			count = 0;
			for (int k = 0; k < kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix()->outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(*kinetics_map_.jacobian_sparsity_pattern_map()->jacobian_matrix(), k); it; ++it)
				{
					{
						edges[count] = Pair(it.row(), it.col());
						count++;
					}
				}
			}
			

			Graph G(nc);
			for (unsigned int i = 0; i < count; ++i)
				boost::add_edge(edges[i].first, edges[i].second, G);

			boost::graph_traits<Graph>::vertex_iterator ui, ui_end;

			boost::property_map<Graph, boost::vertex_degree_t>::type deg = get(boost::vertex_degree, G);
			for (boost::tie(ui, ui_end) = boost::vertices(G); ui != ui_end; ++ui)
				deg[*ui] = boost::degree(*ui, G);

			boost::property_map<Graph, boost::vertex_index_t>::type
				index_map = get(boost::vertex_index, G);

			std::cout << " Original bandwidth: " << boost::bandwidth(G) << std::endl;

			
			// Reverse cuthill_mckee_ordering
			{
				std::vector<Vertex> inv_perm(num_vertices(G));
				std::vector<size_type> perm(num_vertices(G));

				cuthill_mckee_ordering(G, inv_perm.rbegin(), get(boost::vertex_color, G), make_degree_map(G));

				std::cout << " Reverse Cuthill-Mckee ordering:" << std::endl;
				for (size_type c = 0; c != inv_perm.size(); ++c)
					perm[index_map[inv_perm[c]]] = c;
				std::cout << " New bandwidth: " << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0])) << std::endl;

				fOutput << "%Reverse Cuthill-Mckee ordering" << std::endl;
				std::cout << "  ";
				for (std::vector<Vertex>::const_iterator i = inv_perm.begin(); i != inv_perm.end(); ++i)
					fOutput << kinetics_map_.thermodynamics().NamesOfSpecies()[index_map[*i]] << std::endl;
			}
			
		}
	
		fOutput.close();
		
		return true;
	}
}
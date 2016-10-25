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
	StoichiometricMap::StoichiometricMap(const unsigned int nspecies, const unsigned int nreactions)
	{
		number_of_species_ = nspecies;
		number_of_reactions_= nreactions;
                
                verbose_output_ = true;

		isTheStoichiometricMatrixAvailable_ = false;
		isTheReactionOrderMatrixAvailable_ = false;
		areTheContributionOfRateOfFormationMatricesAvailable_ = false;
	}
        
        StoichiometricMap::StoichiometricMap(const unsigned int nspecies, const unsigned int nreactions, bool verbose)
	{
		number_of_species_ = nspecies;
		number_of_reactions_= nreactions;
                
                verbose_output_ = verbose;

		isTheStoichiometricMatrixAvailable_ = false;
		isTheReactionOrderMatrixAvailable_ = false;
		areTheContributionOfRateOfFormationMatricesAvailable_ = false;
	}

	void StoichiometricMap::ReadFromASCIIFile(std::istream& fInput)
	{
		numDir1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numDir2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numDir3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numDir4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numDir5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		numRevTot1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevTot2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevTot3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevTot4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevTot5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		numRevEq1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevEq2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevEq3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevEq4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		numRevEq5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		jDir1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jDir2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jDir3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jDir4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jDir5.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		valueDir5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		jRevTot1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevTot2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevTot3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevTot4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevTot5.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		valueRevTot5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		jRevEq1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevEq2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevEq3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevEq4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		jRevEq5.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		valueRevEq5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		changeOfMoles_.Load(fInput, OPENSMOKE_FORMATTED_FILE);

		fInput >> global_reactions_;

		if (global_reactions_ == false)
		{
			lambda_numDir1	 = numDir1;
			lambda_numDir2	 = numDir2;
			lambda_numDir3	 = numDir3;
			lambda_numDir4	 = numDir4;
			lambda_numDir5	 = numDir5;

			lambda_numRevEq1 = numRevEq1;
			lambda_numRevEq2 = numRevEq2;
			lambda_numRevEq3 = numRevEq3;
			lambda_numRevEq4 = numRevEq4;
			lambda_numRevEq5 = numRevEq5;

			lambda_jDir1	 = jDir1;
			lambda_jDir2	 = jDir2;
			lambda_jDir3	 = jDir3;
			lambda_jDir4	 = jDir4;
			lambda_jDir5	 = jDir5;
			lambda_valueDir5 = valueDir5;

			lambda_jRevEq1   = jRevEq1;
			lambda_jRevEq2   = jRevEq2;
			lambda_jRevEq3   = jRevEq3;
			lambda_jRevEq4   = jRevEq4;
			lambda_jRevEq5   = jRevEq5;
			lambda_valueRevEq5 = valueRevEq5;
		}
		else
		{
			lambda_numDir1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numDir2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numDir3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numDir4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numDir5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

			lambda_numRevEq1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numRevEq2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numRevEq3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numRevEq4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_numRevEq5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

			lambda_jDir1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jDir2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jDir3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jDir4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jDir5.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_valueDir5.Load(fInput, OPENSMOKE_FORMATTED_FILE);

			lambda_jRevEq1.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jRevEq2.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jRevEq3.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jRevEq4.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_jRevEq5.Load(fInput, OPENSMOKE_FORMATTED_FILE);
			lambda_valueRevEq5.Load(fInput, OPENSMOKE_FORMATTED_FILE);
		}

		BuildStoichiometricMatrix();
		BuildReactionOrdersMatrix();
	}

	void StoichiometricMap::CompleteChangeOfMoles(OpenSMOKEVectorBool& isThermodynamicReversible)
	{
		unsigned number_of_thermodynamic_reversible_reactions = 0;
		unsigned number_of_reactions_without_mole_change = 0;
		unsigned number_of_reactions_with_mole_change_plus_one = 0;
		unsigned number_of_reactions_with_mole_change_minus_one = 0;
		for(unsigned int j=1;j<=number_of_reactions_;j++)
		{
			if (isThermodynamicReversible[j] == true)
			{
				number_of_thermodynamic_reversible_reactions++;
				if (changeOfMoles_[j] ==  0)		number_of_reactions_without_mole_change++;
				else if (changeOfMoles_[j] ==  1)	number_of_reactions_with_mole_change_plus_one++;
				else if (changeOfMoles_[j] == -1)	number_of_reactions_with_mole_change_minus_one++;
			}
		}

		
		unsigned number_of_reactions_with_mole_change_other = number_of_thermodynamic_reversible_reactions - number_of_reactions_without_mole_change-
				                                                number_of_reactions_with_mole_change_plus_one - number_of_reactions_with_mole_change_minus_one ;
		ChangeDimensions(number_of_reactions_without_mole_change, &indices_of_reactions_without_change_of_moles_, false);
		ChangeDimensions(number_of_reactions_with_mole_change_plus_one, &indices_of_reactions_with_change_of_moles_plus_one_, false);
		ChangeDimensions(number_of_reactions_with_mole_change_minus_one, &indices_of_reactions_with_change_of_moles_minus_one_, false);
		ChangeDimensions(number_of_reactions_with_mole_change_other, &indices_of_reactions_with_change_of_moles_, false);
		
		number_of_reactions_without_mole_change = 1;
		number_of_reactions_with_mole_change_plus_one = 1;
		number_of_reactions_with_mole_change_minus_one = 1;
		number_of_reactions_with_mole_change_other = 1;
		for(unsigned int j=1;j<=number_of_reactions_;j++)
		{
			if (isThermodynamicReversible[j] == true)
			{
				if (changeOfMoles_[j] == 0)			indices_of_reactions_without_change_of_moles_[number_of_reactions_without_mole_change++] = j;
				else if (changeOfMoles_[j] ==  1)	indices_of_reactions_with_change_of_moles_plus_one_[number_of_reactions_with_mole_change_plus_one++] = j;
				else if (changeOfMoles_[j] == -1)	indices_of_reactions_with_change_of_moles_minus_one_[number_of_reactions_with_mole_change_minus_one++] = j;
				else								indices_of_reactions_with_change_of_moles_[number_of_reactions_with_mole_change_other++] = j;
			}
		}
	}

	void StoichiometricMap::Summary(std::ostream &fOut) const
	{
		unsigned int reactions_thermodynamically_reversible = indices_of_reactions_without_change_of_moles_.Size() + indices_of_reactions_with_change_of_moles_plus_one_.Size() + 
																	  indices_of_reactions_with_change_of_moles_minus_one_.Size()  + indices_of_reactions_with_change_of_moles_.Size();

		fOut << "----------------------------------------------------------------------------" << std::endl;
		fOut << " Reversible reactions (by thermodynamics)    " << reactions_thermodynamically_reversible << std::endl;
		fOut << "----------------------------------------------------------------------------" << std::endl;
		fOut << " Reactions without change of moles:        " << indices_of_reactions_without_change_of_moles_.Size() << " (" << indices_of_reactions_without_change_of_moles_.Size()/std::max(1.,double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << " Reactions with change of moles (+1):      " << indices_of_reactions_with_change_of_moles_plus_one_.Size() <<  " (" << indices_of_reactions_with_change_of_moles_minus_one_.Size()/std::max(1.,double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << " Reactions with change of moles (-1):      " << indices_of_reactions_with_change_of_moles_minus_one_.Size() <<  " (" << indices_of_reactions_with_change_of_moles_plus_one_.Size()/std::max(1.,double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << " Reactions with change of moles (other):   " << indices_of_reactions_with_change_of_moles_.Size() <<  " (" << indices_of_reactions_with_change_of_moles_.Size()/std::max(1.,double(reactions_thermodynamically_reversible))*100. << "%)" << std::endl;
		fOut << std::endl;
	}

	void StoichiometricMap::ProductOfConcentrations(OpenSMOKEVectorDouble& productDirect, OpenSMOKEVectorDouble& productReverse, const OpenSMOKEVectorDouble& c)
	{
		productDirect  =  1.;
		productReverse =  1.;

		double c1, c2, c3, csq;

		unsigned int *jD1, *jD2, *jD3, *jD4, *jD5;
		double *vD5;

		unsigned int *jIE1, *jIE2, *jIE3, *jIE4, *jIE5;
		double *vIE5;

		jD1 = lambda_jDir1.GetHandle();
		jD2 = lambda_jDir2.GetHandle();
		jD3 = lambda_jDir3.GetHandle();
		jD4 = lambda_jDir4.GetHandle();
		jD5 = lambda_jDir5.GetHandle();
		vD5 = lambda_valueDir5.GetHandle();

		jIE1 = lambda_jRevEq1.GetHandle();
		jIE2 = lambda_jRevEq2.GetHandle();
		jIE3 = lambda_jRevEq3.GetHandle();
		jIE4 = lambda_jRevEq4.GetHandle();
		jIE5 = lambda_jRevEq5.GetHandle();
		vIE5 = lambda_valueRevEq5.GetHandle();

		
		for(unsigned int i=1;i<=number_of_species_;i++)
		{
			c1 = c[i];
			c2 = c1 * c1;
			c3 = c2 * c1;
			if(lambda_numDir4[i] != 0 || lambda_numRevEq4[i] != 0)
				csq = std::sqrt(c1);
			
			for(unsigned int k=0;k<lambda_numDir1[i];k++)
			{
				productDirect[*jD1] *= c1;
				jD1++;
			}
			for(unsigned int k=0;k<lambda_numDir2[i];k++)
			{
				productDirect[*jD2] *= c2;
				jD2++;
			}
			for(unsigned int k=0;k<lambda_numDir3[i];k++)
			{
				productDirect[*jD3] *= c3;
				jD3++;
			}
			for(unsigned int k=0;k<lambda_numDir4[i];k++)
			{
				productDirect[*jD4] *= csq;
				jD4++;
			}
			for(unsigned int k=0;k<lambda_numDir5[i];k++)
			{
				productDirect[*jD5] *= std::pow(c1,*vD5);
				jD5++;
				vD5++;
			}

			for(unsigned int k=0;k<lambda_numRevEq1[i];k++)
			{
				productReverse[*jIE1] *= c1;
				jIE1++;
			}
			for(unsigned int k=0;k<lambda_numRevEq2[i];k++)
			{
				productReverse[*jIE2] *= c2;
				jIE2++;
			}
			for(unsigned int k=0;k<lambda_numRevEq3[i];k++)
			{
				productReverse[*jIE3] *= c3;
				jIE3++;
			}
			for(unsigned int k=0;k<lambda_numRevEq4[i];k++)
			{
				productReverse[*jIE4] *= csq;
				jIE4++;
			}
			for(unsigned int k=0;k<lambda_numRevEq5[i];k++)
			{
				productReverse[*jIE5] *= std::pow(c1,*vIE5);
				jIE5++;
				vIE5++;
			}
		}
	}

	void StoichiometricMap::FormationRatesFromReactionRates(OpenSMOKEVectorDouble* R, const OpenSMOKEVectorDouble& r)
	{
		unsigned int* jD1 = jDir1.GetHandle();
		unsigned int* jD2 = jDir2.GetHandle();
		unsigned int* jD3 = jDir3.GetHandle();
		unsigned int* jD4 = jDir4.GetHandle();
		unsigned int* jD5 = jDir5.GetHandle();
		double* vD5 = valueDir5.GetHandle();

		unsigned int* jIT1 = jRevTot1.GetHandle();
		unsigned int* jIT2 = jRevTot2.GetHandle();
		unsigned int* jIT3 = jRevTot3.GetHandle();
		unsigned int* jIT4 = jRevTot4.GetHandle();
		unsigned int* jIT5 = jRevTot5.GetHandle();
		double* vIT5 = valueRevTot5.GetHandle();
	
		for(unsigned int i=1;i<=number_of_species_;i++)
		{
			double rate = 0.;
			for(unsigned int k=0;k<numDir1[i];k++)
				rate -= r[*jD1++];
			for(unsigned int k=0;k<numDir2[i];k++)
			{
				rate -= (r[*jD2]+r[*jD2]);
				*jD2++;
			}
			for(unsigned int k=0;k<numDir3[i];k++)
			{
				rate -= (r[*jD3]+r[*jD3]+r[*jD3]);
				*jD3++;
			}
			for(unsigned int k=0;k<numDir4[i];k++)
				rate -= .5 * r[*jD4++];
			for(unsigned int k=0;k<numDir5[i];k++)
				rate -= (*vD5++) * r[*jD5++];

			for(unsigned int k=0;k<numRevTot1[i];k++)
				rate += r[*jIT1++];
			for(unsigned int k=0;k<numRevTot2[i];k++)
			{
				rate += (r[*jIT2]+r[*jIT2]);
				*jIT2++;
			}
			for(unsigned int k=0;k<numRevTot3[i];k++)
			{
				rate += (r[*jIT3]+r[*jIT3]+r[*jIT3]);
				*jIT3++;
			}
			for(unsigned int k=0;k<numRevTot4[i];k++)
				rate += .5 * r[*jIT4++];
			for(unsigned int k=0;k<numRevTot5[i];k++)
				rate += (*vIT5++) * r[*jIT5++];

			(*R)[i] = rate;
		}
	}

	
	// This version calculates direct and reverse reaction rates
	void StoichiometricMap::ProductionAndDestructionRatesFromReactionRatesGross(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D, const OpenSMOKEVectorDouble& rF, const OpenSMOKEVectorDouble& rB)
	{
		{
			unsigned int* jD1 = jDir1.GetHandle();
			unsigned int* jD2 = jDir2.GetHandle();
			unsigned int* jD3 = jDir3.GetHandle();
			unsigned int* jD4 = jDir4.GetHandle();
			unsigned int* jD5 = jDir5.GetHandle();
			double* vD5 = valueDir5.GetHandle();

			unsigned int* jIT1 = jRevTot1.GetHandle();
			unsigned int* jIT2 = jRevTot2.GetHandle();
			unsigned int* jIT3 = jRevTot3.GetHandle();
			unsigned int* jIT4 = jRevTot4.GetHandle();
			unsigned int* jIT5 = jRevTot5.GetHandle();
			double* vIT5 = valueRevTot5.GetHandle();
	
			for(unsigned int i=1;i<=number_of_species_;i++)
			{
				double destruction = 0.;
				double production = 0.;
				for(unsigned int k=0;k<numDir1[i];k++)
				{
					destruction += rF[*jD1];
					production  += rB[*jD1];
					*jD1++;
				}
				for(unsigned int k=0;k<numDir2[i];k++)
				{
					destruction += (rF[*jD2] + rF[*jD2]);
					production  += (rB[*jD2] + rB[*jD2]);
					*jD2++;
				}
				for(unsigned int k=0;k<numDir3[i];k++)
				{
					destruction += (rF[*jD3] + rF[*jD3] + rF[*jD3]);
					production  += (rB[*jD3] + rB[*jD3] + rB[*jD3]);
					*jD3++;
				}
				for(unsigned int k=0;k<numDir4[i];k++)
				{
					destruction += .5 * rF[*jD4];
					production  += .5 * rB[*jD4];
					*jD4++;
				}
				for(unsigned int k=0;k<numDir5[i];k++)
				{
					destruction += (*vD5) * rF[*jD5];
					production  += (*vD5) * rB[*jD5];
					*vD5++;
					*jD5++;
				}

				
				for(unsigned int k=0;k<numRevTot1[i];k++)
				{
					production  += rF[*jIT1];
					destruction += rB[*jIT1];
					*jIT1++;
				}
				for(unsigned int k=0;k<numRevTot2[i];k++)
				{
					production  += (rF[*jIT2] + rF[*jIT2]);
					destruction += (rB[*jIT2] + rB[*jIT2]);
					*jIT2++;
				}
				for(unsigned int k=0;k<numRevTot3[i];k++)
				{
					production  += (rF[*jIT3] + rF[*jIT3] + rF[*jIT3]);
					destruction += (rB[*jIT3] + rB[*jIT3] + rB[*jIT3]);
					*jIT3++;
				}
				for(unsigned int k=0;k<numRevTot4[i];k++)
				{
					production  += .5 * rF[*jIT4];
					destruction += .5 * rB[*jIT4];
					*jIT4++;
				}
				for(unsigned int k=0;k<numRevTot5[i];k++)
				{
					production  += (*vIT5) * rF[*jIT5];
					destruction += (*vIT5) * rB[*jIT5];
					*vIT5++;
					*jIT5++;
				}

				(*P)[i] = production;
				(*D)[i] = destruction;
			}
		}
	}
	

	void StoichiometricMap::ProductionAndDestructionRatesFromReactionRates(OpenSMOKEVectorDouble* P, OpenSMOKEVectorDouble* D, const OpenSMOKEVectorDouble& r)
	{
		{
			unsigned int* jD1 = jDir1.GetHandle();
			unsigned int* jD2 = jDir2.GetHandle();
			unsigned int* jD3 = jDir3.GetHandle();
			unsigned int* jD4 = jDir4.GetHandle();
			unsigned int* jD5 = jDir5.GetHandle();
			double* vD5 = valueDir5.GetHandle();

			unsigned int* jIT1 = jRevTot1.GetHandle();
			unsigned int* jIT2 = jRevTot2.GetHandle();
			unsigned int* jIT3 = jRevTot3.GetHandle();
			unsigned int* jIT4 = jRevTot4.GetHandle();
			unsigned int* jIT5 = jRevTot5.GetHandle();
			double* vIT5 = valueRevTot5.GetHandle();
	
			for(unsigned int i=1;i<=number_of_species_;i++)
			{
       
				double destruction = 0.;
				double production = 0.;
				for(unsigned int k=0;k<numDir1[i];k++)
				{
					if (r[*jD1] > 0.)	destruction += r[*jD1];
					else                production  -= r[*jD1];
					*jD1++;
				}
				for(unsigned int k=0;k<numDir2[i];k++)
				{
					if (r[*jD2] > 0.)	destruction += (r[*jD2] + r[*jD2]);
					else                production  -= (r[*jD2] + r[*jD2]);
					*jD2++;
				}
				for(unsigned int k=0;k<numDir3[i];k++)
				{
					if (r[*jD3] > 0.)	destruction += (r[*jD3] + r[*jD3] + r[*jD3]);
					else                production  -= (r[*jD3] + r[*jD3] + r[*jD3]);
					*jD3++;
				}
				for(unsigned int k=0;k<numDir4[i];k++)
				{
					if (r[*jD4] > 0.)	destruction += .5 * r[*jD4];
					else                production  -= .5 * r[*jD4];
					*jD4++;
				}
				for(unsigned int k=0;k<numDir5[i];k++)
				{
					if (r[*jD5] > 0.)	destruction += (*vD5) * r[*jD5];
					else                production  -= (*vD5) * r[*jD5];
					*vD5++;
					*jD5++;
				}

				
				for(unsigned int k=0;k<numRevTot1[i];k++)
				{
					if (r[*jIT1] > 0.)	production  += r[*jIT1];
					else                destruction -= r[*jIT1];
					*jIT1++;
				}
				for(unsigned int k=0;k<numRevTot2[i];k++)
				{
					if (r[*jIT2] > 0.)	production  += (r[*jIT2] + r[*jIT2]);
					else                destruction -= (r[*jIT2] + r[*jIT2]);
					*jIT2++;
				}
				for(unsigned int k=0;k<numRevTot3[i];k++)
				{
					if (r[*jIT3] > 0.)	production  += (r[*jIT3] + r[*jIT3] + r[*jIT3]);
					else                destruction -= (r[*jIT3] + r[*jIT3] + r[*jIT3]);
					*jIT3++;
				}
				for(unsigned int k=0;k<numRevTot4[i];k++)
				{
					if (r[*jIT4] > 0.)	production  += .5*r[*jIT4];
					else                destruction -= .5*r[*jIT4];
					*jIT4++;
				}
				for(unsigned int k=0;k<numRevTot5[i];k++)
				{
					if (r[*jIT5] > 0.)	production  += (*vIT5)*r[*jIT5];
					else                destruction -= (*vIT5)*r[*jIT5];
					*vIT5++;
					*jIT5++;
				}

				(*P)[i] = production;
				(*D)[i] = destruction;
			}
		}
	}

	void StoichiometricMap::ReactionEnthalpyAndEntropy(OpenSMOKEVectorDouble& reaction_dh_over_RT, OpenSMOKEVectorDouble& reaction_ds_over_R, const OpenSMOKEVectorDouble& species_h_over_RT, const OpenSMOKEVectorDouble& species_s_over_R)
	{
		unsigned int *jD1 = jDir1.GetHandle();
		unsigned int *jD2 = jDir2.GetHandle();
		unsigned int *jD3 = jDir3.GetHandle();
		unsigned int *jD4 = jDir4.GetHandle();
		unsigned int *jD5 = jDir5.GetHandle();
		double *vD5 = valueDir5.GetHandle();

		unsigned int *jIT1 = jRevTot1.GetHandle();
		unsigned int *jIT2 = jRevTot2.GetHandle();
		unsigned int *jIT3 = jRevTot3.GetHandle();
		unsigned int *jIT4 = jRevTot4.GetHandle();
		unsigned int *jIT5 = jRevTot5.GetHandle();
		double *vIT5 = valueRevTot5.GetHandle();
	
		reaction_dh_over_RT = 0.;
		reaction_ds_over_R = 0.;
			
		for(unsigned int i=1;i<=number_of_species_;i++)
		{
			for(unsigned int k=0;k<numDir1[i];k++)
			{
				reaction_dh_over_RT[*jD1] -= species_h_over_RT[i];
				reaction_ds_over_R[*jD1]  -= species_s_over_R[i];
				jD1++;
			}
			for(unsigned int k=0;k<numDir2[i];k++)
			{
				reaction_dh_over_RT[*jD2] -= (species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jD2]  -= (species_s_over_R[i]  + species_s_over_R[i]);
				jD2++;
			}
			for(unsigned int k=0;k<numDir3[i];k++)
			{
				reaction_dh_over_RT[*jD3] -= (species_h_over_RT[i] + species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jD3]  -= (species_s_over_R[i] + species_s_over_R[i] + species_s_over_R[i]);
				jD3++;
			}
			for(unsigned int k=0;k<numDir4[i];k++)
			{
				reaction_dh_over_RT[*jD4] -= (0.5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jD4]  -= (0.5 * species_s_over_R[i]);
				jD4++;
			}
			for(unsigned int k=0;k<numDir5[i];k++)
			{
				reaction_dh_over_RT[*jD5] -= (*vD5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jD5]  -= (*vD5 * species_s_over_R[i]);
				jD5++;
				vD5++;
			}

			for(unsigned int k=0;k<numRevTot1[i];k++)
			{
				reaction_dh_over_RT[*jIT1] += species_h_over_RT[i];
				reaction_ds_over_R[*jIT1]  += species_s_over_R[i];
				jIT1++;
			}
			for(unsigned int k=0;k<numRevTot2[i];k++)
			{
				reaction_dh_over_RT[*jIT2] += (species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jIT2]  += (species_s_over_R[i] + species_s_over_R[i]);
				jIT2++;
			}
			for(unsigned int k=0;k<numRevTot3[i];k++)
			{
				reaction_dh_over_RT[*jIT3] += (species_h_over_RT[i] + species_h_over_RT[i] + species_h_over_RT[i]);
				reaction_ds_over_R[*jIT3]  += (species_s_over_R[i] + species_s_over_R[i] + species_s_over_R[i]);
				jIT3++;
			}
			for(unsigned int k=0;k<numRevTot4[i];k++)
			{
				reaction_dh_over_RT[*jIT4] += (0.5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jIT4]  += (0.5 * species_s_over_R[i]);
				jIT4++;
			}
			for(unsigned int k=0;k<numRevTot5[i];k++)
			{
				reaction_dh_over_RT[*jIT5] += (*vIT5 * species_h_over_RT[i]);
				reaction_ds_over_R[*jIT5]  += (*vIT5 * species_s_over_R[i]);
				jIT5++;
				vIT5++;
			}
		}
	}



	void StoichiometricMap::BuildStoichiometricMatrix()
	{
		if (isTheStoichiometricMatrixAvailable_ == false)
		{
            if(verbose_output_ == true)
				std::cout	<< " * Building stoichiometry..." << std::endl;

			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList_reactants;
			std::vector<T> tripletList_products;

			unsigned int estimation_of_entries_reactants = 0;
			unsigned int estimation_of_entries_products = 0;
			for(unsigned int i=1;i<=number_of_species_;i++)
			{
				estimation_of_entries_reactants += numDir1[i];
				estimation_of_entries_reactants += numDir2[i];
				estimation_of_entries_reactants += numDir3[i];
				estimation_of_entries_reactants += numDir4[i];
				estimation_of_entries_reactants += numDir5[i];
				estimation_of_entries_products += numRevTot1[i];
				estimation_of_entries_products += numRevTot2[i];
				estimation_of_entries_products += numRevTot3[i];
				estimation_of_entries_products += numRevTot4[i];
				estimation_of_entries_products += numRevTot5[i];
			}
                        
                        if(verbose_output_ == true)
			std::cout	<< "   non-zero stoichiometric coefficients: " << estimation_of_entries_reactants+estimation_of_entries_products << " /" 
						<< number_of_reactions_*number_of_species_ << " (" 
						<< (estimation_of_entries_reactants+estimation_of_entries_products)/double(number_of_reactions_*number_of_species_)*100. << "%)" << std::endl;
		
			//std::cout	<< "Mean number of species per reaction: " << (estimation_of_entries_reactants+estimation_of_entries_products)/number_of_reactions_ << std::endl;

			tripletList_reactants.reserve(estimation_of_entries_reactants);
			tripletList_products.reserve(estimation_of_entries_products);

			unsigned int *jD1 = jDir1.GetHandle();
			unsigned int *jD2 = jDir2.GetHandle();
			unsigned int *jD3 = jDir3.GetHandle();
			unsigned int *jD4 = jDir4.GetHandle();
			unsigned int *jD5 = jDir5.GetHandle();
			double *vD5 = valueDir5.GetHandle();
			
			for(unsigned int i=1;i<=number_of_species_;i++)
			{
				for(unsigned int k=0;k<numDir1[i];k++)
				{
					tripletList_reactants.push_back( T(*jD1-1,i-1,1.) );
					jD1++;
				}
				for(unsigned int k=0;k<numDir2[i];k++)
				{
					tripletList_reactants.push_back( T(*jD2-1,i-1,2.) );
					jD2++;
				}
				for(unsigned int k=0;k<numDir3[i];k++)
				{
					tripletList_reactants.push_back( T(*jD3-1,i-1,3.) );
					jD3++;
				}
				for(unsigned int k=0;k<numDir4[i];k++)
				{
					tripletList_reactants.push_back( T(*jD4-1,i-1,0.5) );
					jD4++;
				}
				for(unsigned int k=0;k<numDir5[i];k++)
				{
					tripletList_reactants.push_back( T(*jD5-1,i-1,(*vD5)) );
					jD5++;
					vD5++;
				}
			}

			unsigned int *jIT1 = jRevTot1.GetHandle();
			unsigned int *jIT2 = jRevTot2.GetHandle();
			unsigned int *jIT3 = jRevTot3.GetHandle();
			unsigned int *jIT4 = jRevTot4.GetHandle();
			unsigned int *jIT5 = jRevTot5.GetHandle();
			double *vIT5 = valueRevTot5.GetHandle();
			
			for(unsigned int i=1;i<=number_of_species_;i++)
			{
				for(unsigned int k=0;k<numRevTot1[i];k++)
				{
					tripletList_products.push_back( T(*jIT1-1,i-1,1.) );
					jIT1++;
				}
				for(unsigned int k=0;k<numRevTot2[i];k++)
				{
					tripletList_products.push_back( T(*jIT2-1,i-1,2.) );
					jIT2++;
				}
				for(unsigned int k=0;k<numRevTot3[i];k++)
				{
					tripletList_products.push_back( T(*jIT3-1,i-1,3.) );
					jIT3++;
				}
				for(unsigned int k=0;k<numRevTot4[i];k++)
				{
					tripletList_products.push_back( T(*jIT4-1,i-1,0.5) );
					jIT4++;
				}
				for(unsigned int k=0;k<numRevTot5[i];k++)
				{
					tripletList_products.push_back( T(*jIT5-1,i-1,*vIT5) );
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
			for (unsigned int i = 1; i <= number_of_species_; i++)
			{
				estimation_of_entries_reactants += lambda_numDir1[i];
				estimation_of_entries_reactants += lambda_numDir2[i];
				estimation_of_entries_reactants += lambda_numDir3[i];
				estimation_of_entries_reactants += lambda_numDir4[i];
				estimation_of_entries_reactants += lambda_numDir5[i];
				estimation_of_entries_products += lambda_numRevEq1[i];
				estimation_of_entries_products += lambda_numRevEq2[i];
				estimation_of_entries_products += lambda_numRevEq3[i];
				estimation_of_entries_products += lambda_numRevEq4[i];
				estimation_of_entries_products += lambda_numRevEq5[i];
			}

			if (verbose_output_ == true)
				std::cout << "   non-zero reaction-order coefficients: " << estimation_of_entries_reactants + estimation_of_entries_products << " /"
				<< number_of_reactions_*number_of_species_ << " ("
				<< (estimation_of_entries_reactants + estimation_of_entries_products) / double(number_of_reactions_*number_of_species_)*100. << "%)" << std::endl;

			//std::cout	<< "Mean number of species per reaction: " << (estimation_of_entries_reactants+estimation_of_entries_products)/number_of_reactions_ << std::endl;

			tripletList_reactants.reserve(estimation_of_entries_reactants);
			tripletList_products.reserve(estimation_of_entries_products);

			unsigned int *jD1 = lambda_jDir1.GetHandle();
			unsigned int *jD2 = lambda_jDir2.GetHandle();
			unsigned int *jD3 = lambda_jDir3.GetHandle();
			unsigned int *jD4 = lambda_jDir4.GetHandle();
			unsigned int *jD5 = lambda_jDir5.GetHandle();
			double *vD5 = lambda_valueDir5.GetHandle();

			for (unsigned int i = 1; i <= number_of_species_; i++)
			{
				for (unsigned int k = 0; k < lambda_numDir1[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD1 - 1, i - 1, 1.));
					jD1++;
				}
				for (unsigned int k = 0; k < lambda_numDir2[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD2 - 1, i - 1, 2.));
					jD2++;
				}
				for (unsigned int k = 0; k < lambda_numDir3[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD3 - 1, i - 1, 3.));
					jD3++;
				}
				for (unsigned int k = 0; k < lambda_numDir4[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD4 - 1, i - 1, 0.5));
					jD4++;
				}
				for (unsigned int k = 0; k < lambda_numDir5[i]; k++)
				{
					tripletList_reactants.push_back(T(*jD5 - 1, i - 1, (*vD5)));
					jD5++;
					vD5++;
				}
			}

			unsigned int *jIE1 = lambda_jRevEq1.GetHandle();
			unsigned int *jIE2 = lambda_jRevEq2.GetHandle();
			unsigned int *jIE3 = lambda_jRevEq3.GetHandle();
			unsigned int *jIE4 = lambda_jRevEq4.GetHandle();
			unsigned int *jIE5 = lambda_jRevEq5.GetHandle();
			double *vIE5 = lambda_valueRevEq5.GetHandle();

			for (unsigned int i = 1; i <= number_of_species_; i++)
			{
				for (unsigned int k = 0; k < lambda_numRevEq1[i]; k++)
				{
					tripletList_products.push_back(T(*jIE1 - 1, i - 1, 1.));
					jIE1++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq2[i]; k++)
				{
					tripletList_products.push_back(T(*jIE2 - 1, i - 1, 2.));
					jIE2++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq3[i]; k++)
				{
					tripletList_products.push_back(T(*jIE3 - 1, i - 1, 3.));
					jIE3++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq4[i]; k++)
				{
					tripletList_products.push_back(T(*jIE4 - 1, i - 1, 0.5));
					jIE4++;
				}
				for (unsigned int k = 0; k < lambda_numRevEq5[i]; k++)
				{
					tripletList_products.push_back(T(*jIE5 - 1, i - 1, *vIE5));
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

	void StoichiometricMap::EquilibriumConstants(OpenSMOKEVectorDouble& Kp, const OpenSMOKEVectorDouble& exp_g_over_RT, const double Patm_over_RT)
	{
		Kp  =  1.;

		double c1, c2, c3, csq;

		unsigned int *jD1 = jDir1.GetHandle();
		unsigned int *jD2 = jDir2.GetHandle();
		unsigned int *jD3 = jDir3.GetHandle();
		unsigned int *jD4 = jDir4.GetHandle();
		unsigned int *jD5 = jDir5.GetHandle();
		double *vD5 = valueDir5.GetHandle();

		unsigned int *jIT1 = jRevTot1.GetHandle();
		unsigned int *jIT2 = jRevTot2.GetHandle();
		unsigned int *jIT3 = jRevTot3.GetHandle();
		unsigned int *jIT4 = jRevTot4.GetHandle();
		unsigned int *jIT5 = jRevTot5.GetHandle();
		double *vIT5 = valueRevTot5.GetHandle();

		
		for(unsigned int i=1;i<=number_of_species_;i++)
		{
			c1 = exp_g_over_RT[i];
			c2 = c1 * c1;
			c3 = c2 * c1;
			if(numDir4[i] != 0 || numRevTot4[i] != 0)
				csq = std::sqrt(c1);
			
			for(unsigned int k=0;k<numDir1[i];k++)
			{
				Kp[*jD1] /= c1;
				jD1++;
			}
			for(unsigned int k=0;k<numDir2[i];k++)
			{
				Kp[*jD2] /= c2;
				jD2++;
			}
			for(unsigned int k=0;k<numDir3[i];k++)
			{
				Kp[*jD3] /= c3;
				jD3++;
			}
			for(unsigned int k=0;k<numDir4[i];k++)
			{
				Kp[*jD4] /= csq;
				jD4++;
			}
			for(unsigned int k=0;k<numDir5[i];k++)
			{
				Kp[*jD5] /= std::pow(c1,*vD5);
				jD5++;
				vD5++;
			}

			for(unsigned int k=0;k<numRevTot1[i];k++)
			{
				Kp[*jIT1] *= c1;
				jIT1++;
			}
			for(unsigned int k=0;k<numRevTot2[i];k++)
			{
				Kp[*jIT2] *= c2;
				jIT2++;
			}
			for(unsigned int k=0;k<numRevTot3[i];k++)
			{
				Kp[*jIT3] *= c3;
				jIT3++;
			}
			for(unsigned int k=0;k<numRevTot4[i];k++)
			{
				Kp[*jIT4] *= csq;
				jIT4++;
			}
			for(unsigned int k=0;k<numRevTot5[i];k++)
			{
				Kp[*jIT5] *= std::pow(c1,*vIT5);
				jIT5++;
				vIT5++;
			}
		}

		for(int j=1;j<=indices_of_reactions_with_change_of_moles_plus_one_.Size();j++)
		{
			const unsigned int k = indices_of_reactions_with_change_of_moles_plus_one_[j];
			Kp[k] /= Patm_over_RT;
		}
		for(int j=1;j<=indices_of_reactions_with_change_of_moles_minus_one_.Size();j++)
		{
			const unsigned int k = indices_of_reactions_with_change_of_moles_minus_one_[j];
			Kp[k] *= Patm_over_RT;
		}
		for(int j=1;j<=indices_of_reactions_with_change_of_moles_.Size();j++)
		{
			const unsigned int k = indices_of_reactions_with_change_of_moles_[j];
			Kp[k] *= std::pow(Patm_over_RT, -changeOfMoles_[k]);
		}
	}

	void StoichiometricMap::RateOfProductionAnalysis(const OpenSMOKEVectorDouble& r, const bool iNormalize)
	{
		BuildStoichiometricMatrix();

		if (areTheContributionOfRateOfFormationMatricesAvailable_ == false)
		{
			typedef Eigen::Triplet<double> T;
			std::vector<T> tripletList_;
			tripletList_.reserve(stoichiometric_matrix_reactants_.nonZeros()+stoichiometric_matrix_products_.nonZeros());

			// Reactants
			for (int k=0; k<stoichiometric_matrix_reactants_.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_reactants_,k); it; ++it)
					tripletList_.push_back( T(it.row(),it.col(),-it.value()) );
			}
		
			// Products	
			for (int k=0; k<stoichiometric_matrix_products_.outerSize(); ++k)
			{
				for (Eigen::SparseMatrix<double>::InnerIterator it(stoichiometric_matrix_products_,k); it; ++it)
					tripletList_.push_back( T(it.row(),it.col(),it.value()) );
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
		for (int k=0; k<Cd.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
				it.valueRef() = 0.;
		for (int k=0; k<Cp.outerSize(); ++k)
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
				it.valueRef() = 0.;
		
		// Fill (the values are not normalized)
		for (int k=0; k<stoichiometric_matrix_.outerSize(); ++k)
		{
			Eigen::SparseMatrix<double>::InnerIterator itCd(Cd,k);
			Eigen::SparseMatrix<double>::InnerIterator itCp(Cp,k);
			for (Eigen::SparseMatrix<double>::InnerIterator itStoichiometry(stoichiometric_matrix_,k); itStoichiometry; ++itStoichiometry)
			{
				const double value = itStoichiometry.value() * r[itStoichiometry.row()+1];
				if (value >= 0.)	itCp.valueRef() = value;
				else				itCd.valueRef() = value;
				++itCp;
				++itCd;
			}
		}

		if (iNormalize == true)
		{
			// Reactants
			for (int k=0; k<Cd.outerSize(); ++k)
			{
				double sum = 0.;
				for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
					sum += it.value();

				if (sum == 0.)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
						it.valueRef() = 0;
				}
				else
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
						it.valueRef() =  it.value() / sum;
				}
			}

			// Products
			for (int k=0; k<Cp.outerSize(); ++k)
			{
				double sum = 0.;
				for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
					sum += it.value();

				if (sum == 0.)
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
						it.valueRef() = 0;
				}
				else
				{
					for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
						it.valueRef() =  it.value() / sum;
				}
			}
		}
	}

	void StoichiometricMap::RateOfProductionAnalysis(const OpenSMOKEVectorDouble& rf, const OpenSMOKEVectorDouble& rb)
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
					const double value = itStoichiometry.value() * rf[itStoichiometry.row() + 1];
					if (value >= 0.)	itCp.valueRef() = value;
					else				itCd.valueRef() = value;
				}

				// Destruction
				{
					const double value = itStoichiometry.value() * rb[itStoichiometry.row() + 1];
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
		for (int k=0; k<Cd.outerSize(); ++k)
		{
			// Calculates the sum of destruction rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			fout << sum << std::endl;

			// Writes the reactions involved
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
				fout << it.row() << " ";
			fout << std::endl;

			// Writes the coefficients (not normalized)
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
				fout << it.value() << " ";
			fout << std::endl;
		}

		for (int k=0; k<Cp.outerSize(); ++k)
		{
			// Calculates the sum of production rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			fout << sum << std::endl;

			// Writes the reactions involved
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
				fout << it.row() << " ";
			fout << std::endl;

			// Writes the coefficients (not normalized)
			fout << count << " ";
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
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

		for (int k=0; k<Cd.outerSize(); ++k)
		{
			// Calculates the sum of destruction rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			ropa.destruction_rates[k] = sum;

			// Writes the reactions involved
			ropa.destruction_reaction_indices[k].resize(count);
			unsigned int j1=0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
				ropa.destruction_reaction_indices[k][j1++] = it.row();
			
			// Writes the coefficients (not normalized)
			ropa.destruction_coefficients[k].resize(count);
			unsigned int j2=0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cd,k); it; ++it)
				ropa.destruction_coefficients[k][j2++] =  it.value();
		}

		for (int k=0; k<Cp.outerSize(); ++k)
		{
			// Calculates the sum of production rates and the number of reactions
			double sum = 0.;
			unsigned int count = 0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
			{
				sum += it.value();
				count++;
			}

			// Writes the sum of destruction rates
			ropa.production_rates[k] = sum;

			// Writes the reactions involved
			ropa.production_reaction_indices[k].resize(count);
			unsigned int j1=0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
				ropa.production_reaction_indices[k][j1++] = it.row();

			// Writes the coefficients (not normalized)
			ropa.production_coefficients[k].resize(count);
			unsigned int j2=0;
			for (Eigen::SparseMatrix<double>::InnerIterator it(Cp,k); it; ++it)
				ropa.production_coefficients[k][j2++] = it.value();
		}
		
	}	
}


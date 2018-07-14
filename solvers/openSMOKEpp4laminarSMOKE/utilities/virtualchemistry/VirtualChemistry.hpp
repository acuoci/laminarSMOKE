/*-----------------------------------------------------------------------*\
|    ___                   ____  __  __  ___  _  _______                  |
|   / _ \ _ __   ___ _ __ / ___||  \/  |/ _ \| |/ / ____| _     _         |
|  | | | | '_ \ / _ \ '_ \\___ \| |\/| | | | | ' /|  _| _| |_ _| |_       |
|  | |_| | |_) |  __/ | | |___) | |  | | |_| | . \| |__|_   _|_   _|      |
|   \___/| .__/ \___|_| |_|____/|_|  |_|\___/|_|\_\_____||_|   |_|        |
|        |_|                                                              |
|                                                                         |
|   Authors: Alberto Cuoci, Giampaolo Maio, Benoit Fiorina                |
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
|   Copyright(C) 2018 Alberto Cuoci                                       |
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

#include <queue>

namespace OpenSMOKE
{
	double POW(const double C_mol_cm3, const double lambda, const double conversion = 1000.);

	VirtualChemistry::VirtualChemistry(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, const bool is_active) :
	thermodynamicsMap_(thermodynamicsMap)
	{
		if (is_active == true)
		{
			// Initializing
			ns_ = thermodynamicsMap_.NumberOfSpecies();
			ns_main_ = 8;
			NSVector_.resize(ns_);
			MW_.resize(ns_);
			Le_.resize(ns_);
			Le_.setConstant(1.);
			iReactions_ = true;
			iSubMechanism_CO_ = false;
			on_the_fly_optimization_ = "none";

			// Indices of products
			I_index_ = thermodynamicsMap_.IndexOfSpecies("I") - 1;
			P1_index_ = thermodynamicsMap_.IndexOfSpecies("P1") - 1;
			P2_index_ = thermodynamicsMap_.IndexOfSpecies("P2") - 1;
			P3_index_ = thermodynamicsMap_.IndexOfSpecies("P3") - 1;
			P4_index_ = thermodynamicsMap_.IndexOfSpecies("P4") - 1;

			A1_ = 1.5384796696859927E+18;
			A2_ = 3.9225000838284247E+18;
			E1_ = 3.5075174514637045E+04;
			E2_ = 8.5680603523860715E+04;
			nuF_1_ = 1.7099829697856470;
			nuOX_1_ = 0.8686294702424591;
			nuI_2_ = 2.3399221348980621;
		}
	}

	VirtualChemistry::VirtualChemistry(OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap, OpenSMOKE::OpenSMOKE_Dictionary& dictionary) :
	VirtualChemistry(thermodynamicsMap)
	{
		SetupFromDictionary(dictionary);
	}

	void VirtualChemistry::SetupFromDictionary(OpenSMOKE::OpenSMOKE_Dictionary& dictionary)
	{
		Grammar_VirtualChemistry grammar;
		dictionary.SetGrammar(grammar);

		// Read table
		{
			boost::filesystem::path path_table;
			if (dictionary.CheckOption("@Table") == true)
				dictionary.ReadPath("@Table", path_table);

			table_main_.Setup(path_table);
		}

		// Read fuel name and molecular weight
		{
			std::string name;
			if (dictionary.CheckOption("@FuelName") == true)
				dictionary.ReadString("@FuelName", name);

			double value;
			if (dictionary.CheckOption("@FuelMW") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@FuelMW", value, units);
				if (units == "kg/kmol")		value *= 1.;
				else if (units == "g/mol")	value *= 1;
				else OpenSMOKE::FatalErrorMessage("Unknown molecular weight units");
			}

			SetFuel(name, value);
		}

		// Read oxidizer name and molecular weight
		{
			std::string name;
			if (dictionary.CheckOption("@OxidizerName") == true)
				dictionary.ReadString("@OxidizerName", name);

			double value;
			if (dictionary.CheckOption("@OxidizerMW") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@OxidizerMW", value, units);
				if (units == "kg/kmol")		value *= 1.;
				else if (units == "g/mol")	value *= 1;
				else OpenSMOKE::FatalErrorMessage("Unknown molecular weight units");
			}

			SetOxidizer(name, value);
		}

		// Read inert name and molecular weight
		{
			std::string name;
			if (dictionary.CheckOption("@InertName") == true)
				dictionary.ReadString("@InertName", name);

			double value;
			if (dictionary.CheckOption("@InertMW") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@InertMW", value, units);
				if (units == "kg/kmol")		value *= 1.;
				else if (units == "g/mol")	value *= 1;
				else OpenSMOKE::FatalErrorMessage("Unknown molecular weight units");
			}

			SetInert(name, value);
		}

		// Read viscosity correlation
		{
			double mu0;
			if (dictionary.CheckOption("@Viscosity_mu0") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@Viscosity_mu0", mu0, units);
				if (units == "kg/m/s")		mu0 *= 1.;
				else if (units == "Pa.s")	mu0 *= 1;
				else OpenSMOKE::FatalErrorMessage("Unknown viscosity units");
			}

			double T0;
			if (dictionary.CheckOption("@Viscosity_T0") == true)
			{
				std::string units;
				dictionary.ReadMeasure("@Viscosity_T0", T0, units);
				if (units == "K")	T0 *= 1.;
				else if (units == "C")	T0 += 273.15;
				else OpenSMOKE::FatalErrorMessage("Unknown temperature units");
			}

			double Beta0;
			if (dictionary.CheckOption("@Viscosity_Beta0") == true)
				dictionary.ReadDouble("@Viscosity_Beta0", Beta0);

			double Pr0;
			if (dictionary.CheckOption("@Viscosity_Pr0") == true)
				dictionary.ReadDouble("@Viscosity_Pr0", Pr0);

			SetTransportProperties(mu0, T0, Beta0, Pr0);		
		}

		// Rections on/off
		{
			if (dictionary.CheckOption("@Reactions") == true)
				dictionary.ReadBool("@Reactions", iReactions_);
		}

		// CO submechanism
		{
			if (dictionary.CheckOption("@SubMechanism_CO") == true)
			{
				dictionary.ReadBool("@SubMechanism_CO", iSubMechanism_CO_);

				if (iSubMechanism_CO_ == true)
				{
					CO_index_ = thermodynamicsMap_.IndexOfSpecies("CO") - 1;
					V1_index_ = thermodynamicsMap_.IndexOfSpecies("V1") - 1;
					V2_index_ = thermodynamicsMap_.IndexOfSpecies("V2") - 1;

					if (CO_index_ < ns_main_ || V1_index_ < ns_main_ || V2_index_ < ns_main_)
					{
						std::cout << "Virtual Chemistry: the first 8 species must be the main species" << std::endl;
						std::cout << "                   please revise the order of CO, V1, and V2 species" << std::endl;
						abort();
					}

					// Read CO table
					{
						if (dictionary.CheckOption("@Table_CO") == true)
						{
							boost::filesystem::path path_table_co;
							dictionary.ReadPath("@Table_CO", path_table_co);
							table_co_.Setup(path_table_co);
						}
						else
						{
							std::cout << "Virtual Chemistry: a look-up table for CO submechanism must be provided!" << std::endl;
							std::cout << "                   please use the @Table_CO option" << std::endl;
							abort();
						}
					}
				}
			}
		}

		// NO submechanism
		{
			if (dictionary.CheckOption("@SubMechanism_NO") == true)
			{
				dictionary.ReadBool("@SubMechanism_NO", iSubMechanism_NO_);

				if (iSubMechanism_NO_ == true)
				{
					NO_index_ = thermodynamicsMap_.IndexOfSpecies("NO") - 1;
					W1_index_ = thermodynamicsMap_.IndexOfSpecies("W1") - 1;
					W2_index_ = thermodynamicsMap_.IndexOfSpecies("W2") - 1;
					W3_index_ = thermodynamicsMap_.IndexOfSpecies("W3") - 1;

					if (NO_index_ < ns_main_ || W1_index_ < ns_main_ || W2_index_ < ns_main_ || W3_index_ < ns_main_)
					{
						std::cout << "Virtual Chemistry: the first 8 species must be the main species" << std::endl;
						std::cout << "                   please revise the order of NO, W1, W2, and W3 species" << std::endl;
						abort();
					}

					// Read NO table
					{
						if (dictionary.CheckOption("@Table_NO") == true)
						{
							boost::filesystem::path path_table_no;
							dictionary.ReadPath("@Table_NO", path_table_no);
							const boost::filesystem::path path_table_no_1 = path_table_no.string() + ".1";
							const boost::filesystem::path path_table_no_2 = path_table_no.string() + ".2";
							const boost::filesystem::path path_table_no_3 = path_table_no.string() + ".3";
							const boost::filesystem::path path_table_no_4 = path_table_no.string() + ".4";
							const boost::filesystem::path path_table_no_5 = path_table_no.string() + ".5";
							table_no_1_.Setup(path_table_no_1);
							table_no_2_.Setup(path_table_no_2);
							table_no_3_.Setup(path_table_no_3);
							table_no_4_.Setup(path_table_no_4);
							table_no_5_.Setup(path_table_no_5);
						}
						else
						{
							std::cout << "Virtual Chemistry: a look-up table for NO submechanism must be provided!" << std::endl;
							std::cout << "                   please use the @Table_NO option" << std::endl;
							abort();
						}
					}
				}
			}
		}

		// Read version
		{
			int version;
			if (dictionary.CheckOption("@Version") == true)
				dictionary.ReadInt("@Version", version);

			SetVersion(version);
		}
	}

	void VirtualChemistry::SetTableMain(const boost::filesystem::path path_table)
	{
		table_main_.Setup(path_table);
	}

	void VirtualChemistry::SetTableCO(const boost::filesystem::path path_table)
	{
		iSubMechanism_CO_ = true;

		CO_index_ = thermodynamicsMap_.IndexOfSpecies("CO") - 1;
		V1_index_ = thermodynamicsMap_.IndexOfSpecies("V1") - 1;
		V2_index_ = thermodynamicsMap_.IndexOfSpecies("V2") - 1;

		if (CO_index_ < ns_main_ || V1_index_ < ns_main_ || V2_index_ < ns_main_)
		{
			std::cout << "Virtual Chemistry: the first 8 species must be the main species" << std::endl;
			std::cout << "                   please revise the order of CO, V1, and V2 species" << std::endl;
			abort();
		}

		table_co_.Setup(path_table);
	}

	void VirtualChemistry::SetTableNO(	const boost::filesystem::path path_table_1, const boost::filesystem::path path_table_2,
										const boost::filesystem::path path_table_3, const boost::filesystem::path path_table_4,
										const boost::filesystem::path path_table_5)
	{
		iSubMechanism_NO_ = true;

		NO_index_ = thermodynamicsMap_.IndexOfSpecies("NO") - 1;
		W1_index_ = thermodynamicsMap_.IndexOfSpecies("W1") - 1;
		W2_index_ = thermodynamicsMap_.IndexOfSpecies("W2") - 1;
		W3_index_ = thermodynamicsMap_.IndexOfSpecies("W3") - 1;

		if (NO_index_ < ns_main_ || W1_index_ < ns_main_ || W2_index_ < ns_main_ || W3_index_ < ns_main_)
		{
			std::cout << "Virtual Chemistry: the first 8 species must be the main species" << std::endl;
			std::cout << "                   please revise the order of NO, W1, W2, and W3 species" << std::endl;
			abort();
		}

		table_no_1_.Setup(path_table_1);
		table_no_2_.Setup(path_table_2);
		table_no_3_.Setup(path_table_3);
		table_no_4_.Setup(path_table_4);
		table_no_5_.Setup(path_table_5);
	}

	void VirtualChemistry::SetReactions(const bool flag)
	{
		iReactions_ = flag;
	}

	void VirtualChemistry::SetOnTheFlyOptimization(const std::string flag)
	{
		on_the_fly_optimization_ = flag;

		const double YN2 = 7.246638000E-01;

		if (on_the_fly_optimization_ == "main")
		{
			// Interpolation (table main)
			table_main_.Interpolation(YN2);
			alpha1_ = table_main_.interpolated()(3);
			alpha2_ = table_main_.interpolated()(4);
			alpha3_ = table_main_.interpolated()(5);
			alpha4_ = table_main_.interpolated()(6);
		}
		else if (on_the_fly_optimization_ == "NO")
		{
			// Interpolation (table 2)
			table_no_2_.Interpolation(YN2);
			Kc5_NO_ = table_no_2_.interpolated()(0);
			beta_NO_ = table_no_2_.interpolated()(1);
			alpha_NO_ = table_no_2_.interpolated()(2);
			gamma_NO_ = table_no_2_.interpolated()(3);

			// Interpolation (table 3)
			table_no_3_.Interpolation(YN2);
			E3_NO_ = table_no_3_.interpolated()(0);
			E4_NO_ = table_no_3_.interpolated()(1);
			E5_NO_ = table_no_3_.interpolated()(2);
			E6_NO_ = table_no_3_.interpolated()(3);
			E7_NO_ = table_no_3_.interpolated()(4);

			// Interpolation (table 4)
			table_no_4_.Interpolation(YN2);
			nuF_3_NO_ = table_no_4_.interpolated()(0);
			nuOX_3_NO_ = table_no_4_.interpolated()(1);
			nuW1_4_NO_ = table_no_4_.interpolated()(2);
			nuNO_5f_NO_ = table_no_4_.interpolated()(3);
			nuW2_5f_NO_ = table_no_4_.interpolated()(4);
			nuNO_5b_NO_ = table_no_4_.interpolated()(5);
			nuW2_5b_NO_ = table_no_4_.interpolated()(6);
			nuW3_6_NO_ = table_no_4_.interpolated()(7);
			nuW3_7_NO_ = table_no_4_.interpolated()(8);

			// Interpolation (table 5)
			table_no_5_.Interpolation(YN2);
			A3_NO_ = table_no_5_.interpolated()(0);
			A4_NO_ = table_no_5_.interpolated()(1);
			A5_NO_ = table_no_5_.interpolated()(2);
			A6_NO_ = table_no_5_.interpolated()(3);
			A7_NO_ = table_no_5_.interpolated()(4);
		}
	}

	void VirtualChemistry::GetOriginalParameters(Eigen::VectorXd& parameters)
	{
		if (on_the_fly_optimization_ == "main")
		{
			parameters.resize(4);

			const double YN2 = 7.246638000E-01;

			// Interpolation (table main)
			table_main_.Interpolation(YN2);
			parameters(0) = table_main_.interpolated()(3);
			parameters(1) = table_main_.interpolated()(4);
			parameters(2) = table_main_.interpolated()(5);
			parameters(3) = table_main_.interpolated()(6);
		}
		else if (on_the_fly_optimization_ == "CO")
		{
			// TODO
		}
		else if (on_the_fly_optimization_ == "NO")
		{
			parameters.resize(23);

			const double YN2 = 7.246638000E-01;

			// Interpolation (table 2)
			table_no_2_.Interpolation(YN2);
			parameters(10) = table_no_2_.interpolated()(0);		// Kc5_NO_
			parameters(11) = table_no_2_.interpolated()(2);		// alpha_NO
			parameters(12) = table_no_2_.interpolated()(1);		// beta_NO
			parameters(13) = table_no_2_.interpolated()(3);		// gamma_NO

			// Interpolation (table 3)
			table_no_3_.Interpolation(YN2);
			parameters(5) = table_no_3_.interpolated()(0);		// E3
			parameters(6) = table_no_3_.interpolated()(1);		// E4
			parameters(7) = table_no_3_.interpolated()(2);		// E5
			parameters(8) = table_no_3_.interpolated()(3);		// E6
			parameters(9) = table_no_3_.interpolated()(4);		// E7

			// Interpolation (table 4)
			table_no_4_.Interpolation(YN2);
			parameters(14) = table_no_4_.interpolated()(0);		// stoichiometric coefficients
			parameters(15) = table_no_4_.interpolated()(1);
			parameters(16) = table_no_4_.interpolated()(2);
			parameters(17) = table_no_4_.interpolated()(3);
			parameters(18) = table_no_4_.interpolated()(4);
			parameters(19) = table_no_4_.interpolated()(5);
			parameters(20) = table_no_4_.interpolated()(6);
			parameters(21) = table_no_4_.interpolated()(7);
			parameters(22) = table_no_4_.interpolated()(8);

			// Interpolation (table 5)
			table_no_5_.Interpolation(YN2);
			parameters(0) = table_no_5_.interpolated()(0);		// A1
			parameters(1) = table_no_5_.interpolated()(1);		// A2
			parameters(2) = table_no_5_.interpolated()(2);		// A3
			parameters(3) = table_no_5_.interpolated()(3);		// A4
			parameters(4) = table_no_5_.interpolated()(4);		// A5
		}

	}

	void VirtualChemistry::GetParameters(Eigen::VectorXd& parameters, const std::vector<bool>& flag) const
	{
		unsigned int k = 0;

		if (on_the_fly_optimization_ == "main")
		{
			if (flag[0] == true) parameters(k++) = alpha1_;			// #0
			if (flag[1] == true) parameters(k++) = alpha2_;			// #1
			if (flag[2] == true) parameters(k++) = alpha3_;			// #2
			if (flag[3] == true) parameters(k++) = alpha4_;			// #3
		}
		else if (on_the_fly_optimization_ == "CO")
		{
			// TODO
		}
		else if (on_the_fly_optimization_ == "NO")
		{
			if (flag[0] == true) parameters(k++) = A3_NO_;			// #0
			if (flag[1] == true) parameters(k++) = A4_NO_;			// #1
			if (flag[2] == true) parameters(k++) = A5_NO_;			// #2
			if (flag[3] == true) parameters(k++) = A6_NO_;			// #3
			if (flag[4] == true) parameters(k++) = A7_NO_;			// #4

			if (flag[5] == true) parameters(k++) = E3_NO_;			// #5
			if (flag[6] == true) parameters(k++) = E4_NO_;			// #6
			if (flag[7] == true) parameters(k++) = E5_NO_;			// #7
			if (flag[8] == true) parameters(k++) = E6_NO_;			// #8
			if (flag[9] == true) parameters(k++) = E7_NO_;			// #9
			if (flag[10] == true) parameters(k++) = Kc5_NO_;		// #10

			if (flag[11] == true) parameters(k++) = alpha_NO_;		// #11
			if (flag[12] == true) parameters(k++) = beta_NO_;		// #12
			if (flag[13] == true) parameters(k++) = gamma_NO_;		// #13

			if (flag[14] == true) parameters(k++) = nuF_3_NO_;		// #14
			if (flag[15] == true) parameters(k++) = nuOX_3_NO_;		// #15
			if (flag[16] == true) parameters(k++) = nuW1_4_NO_;		// #16
			if (flag[17] == true) parameters(k++) = nuNO_5f_NO_;	// #17
			if (flag[18] == true) parameters(k++) = nuW2_5f_NO_;	// #18
			if (flag[19] == true) parameters(k++) = nuNO_5b_NO_;	// #19
			if (flag[20] == true) parameters(k++) = nuW2_5b_NO_;	// #20
			if (flag[21] == true) parameters(k++) = nuW3_6_NO_;		// #21
			if (flag[22] == true) parameters(k++) = nuW3_7_NO_;		// #22
		}
	}

	void VirtualChemistry::SetParameters(const Eigen::VectorXd& parameters, const std::vector<bool>& flag)
	{
		unsigned int k = 0;

		if (on_the_fly_optimization_ == "main")
		{
			if (flag[0] == true)	alpha1_ = parameters(k++);		// #0
			if (flag[1] == true)	alpha2_ = parameters(k++);		// #1
			if (flag[2] == true)	alpha3_ = parameters(k++);		// #2
			if (flag[3] == true)	alpha4_ = parameters(k++);		// #3
		}
		else if (on_the_fly_optimization_ == "CO")
		{
			// TODO
		}
		else if (on_the_fly_optimization_ == "NO")
		{
			if (flag[0] == true)	A3_NO_ = parameters(k++);		// #0
			if (flag[1] == true)	A4_NO_ = parameters(k++);		// #1
			if (flag[2] == true)	A5_NO_ = parameters(k++);		// #2
			if (flag[3] == true)	A6_NO_ = parameters(k++);		// #3
			if (flag[4] == true)	A7_NO_ = parameters(k++);		// #4

			if (flag[5] == true)	E3_NO_ = parameters(k++);		// #5
			if (flag[6] == true)	E4_NO_ = parameters(k++);		// #6
			if (flag[7] == true)	E5_NO_ = parameters(k++);		// #7
			if (flag[8] == true)	E6_NO_ = parameters(k++);		// #8
			if (flag[9] == true)	E7_NO_ = parameters(k++);		// #9
			if (flag[10] == true)	Kc5_NO_ = parameters(k++);		// #10

			if (flag[11] == true)	alpha_NO_ = parameters(k++);	// #11
			if (flag[12] == true)	beta_NO_ = parameters(k++);		// #12
			if (flag[13] == true)	gamma_NO_ = parameters(k++);	// #13

			if (flag[14] == true)	nuF_3_NO_ = parameters(k++);	// #14
			if (flag[15] == true)	nuOX_3_NO_ = parameters(k++);	// #15
			if (flag[16] == true)	nuW1_4_NO_ = parameters(k++);	// #16
			if (flag[17] == true)	nuNO_5f_NO_ = parameters(k++);	// #17
			if (flag[18] == true)	nuW2_5f_NO_ = parameters(k++);	// #18
			if (flag[19] == true)	nuNO_5b_NO_ = parameters(k++);	// #19
			if (flag[20] == true)	nuW2_5b_NO_ = parameters(k++);	// #20
			if (flag[21] == true)	nuW3_6_NO_ = parameters(k++);	// #21
			if (flag[22] == true)	nuW3_7_NO_ = parameters(k++);	// #22
		}
	}

	void VirtualChemistry::SetVersion(const unsigned int value)
	{
		if (value != 170911 && value != 171013)
		{
			std::cout << "Virtual Chemistry: error in choosing the version type: " << value << std::endl;
			std::cout << "                   available: 170911 | 171013" << std::endl;
			abort();
		}

		iVersion_ = value;

		if (iVersion_ == 170911)
		{
			MW_[fuel_index_] = fuel_mw_;
			MW_[oxidizer_index_] = oxidizer_mw_;
			MW_[inert_index_] = inert_mw_;
			MW_[P1_index_] = 15.336;
			MW_[P2_index_] = 16.345;
			MW_[P3_index_] = 93.850;
			MW_[P4_index_] = 37.960;
		}

		if (iVersion_ == 171013)
		{
			MW_[fuel_index_] = fuel_mw_;
			MW_[oxidizer_index_] = oxidizer_mw_;
			MW_[inert_index_] = inert_mw_;
			MW_[P1_index_] = 15.336;
			MW_[P2_index_] = 16.345;
			MW_[P3_index_] = 93.850;
			MW_[P4_index_] = 37.960;
		}

		if (iSubMechanism_CO_ == true)
		{
			MW_[CO_index_] = 28.;
			MW_[V1_index_] = 28.;
			MW_[V2_index_] = 28.;
		}

		if (iSubMechanism_NO_ == true)
		{
			MW_[NO_index_] = 30.;
			MW_[W1_index_] = 30.;
			MW_[W2_index_] = 30.;
			MW_[W3_index_] = 30.;
		}
	}

	void VirtualChemistry::SetFuel(const std::string name, const double mw)
	{
		fuel_name_  = name;
		fuel_mw_    = mw;
		fuel_index_ = thermodynamicsMap_.IndexOfSpecies(fuel_name_)-1; 
	}

	void VirtualChemistry::SetOxidizer(const std::string name, const double mw)
	{
		oxidizer_name_  = name;
		oxidizer_mw_    = mw;
		oxidizer_index_ = thermodynamicsMap_.IndexOfSpecies(oxidizer_name_)-1; 
	}

	void VirtualChemistry::SetInert(const std::string name, const double mw)
	{
		inert_name_  = name;
		inert_mw_    = mw;
		inert_index_ = thermodynamicsMap_.IndexOfSpecies(inert_name_)-1; 
	}

	void VirtualChemistry::SetTransportProperties(const double mu0, const double T0, const double Beta0, const double Pr0)
	{
		mu0_   = mu0;
		T0_    = T0;
		Beta0_ = Beta0;
		Pr0_   = Pr0;
	}

	void VirtualChemistry::SetLewisNumber(const unsigned int i, const double Le)
	{
		Le_(i) = Le;
	}

	double VirtualChemistry::MWMix(double* Y)
	{
		double mw = 0.;

		if (iVersion_ == 170911)
		{
			const double mw_table = table_main_.Interpolation(Y[inert_index_], 0);
			mw = 1. / (Y[fuel_index_] / fuel_mw_ + Y[oxidizer_index_] / oxidizer_mw_ + Y[inert_index_] / inert_mw_ + Y[I_index_] / mw_table +
					Y[P1_index_] / MW_[P1_index_] + Y[P2_index_] / MW_[P2_index_] + Y[P3_index_] / MW_[P3_index_] + Y[P4_index_] / MW_[P4_index_]);
		}
		else if (iVersion_ == 171013)
		{
			const double mw_table = table_main_.Interpolation(Y[inert_index_], 0);
			mw = 1./( Y[fuel_index_]/fuel_mw_ + Y[oxidizer_index_]/oxidizer_mw_ + Y[inert_index_]/inert_mw_ + Y[I_index_]/mw_table +
					Y[P1_index_]/MW_[P1_index_] + Y[P2_index_]/ MW_[P2_index_] + Y[P3_index_]/MW_[P3_index_] + Y[P4_index_]/MW_[P4_index_]);
		}
		
		return mw;     
	}

	void VirtualChemistry::MoleFractions(const double mw, double* Y, double* X)
	{
		if (iVersion_ == 170911)
		{
			X[fuel_index_] = Y[fuel_index_] * mw / fuel_mw_;
			X[oxidizer_index_] = Y[oxidizer_index_] * mw / oxidizer_mw_;
			X[inert_index_] = Y[inert_index_] * mw / inert_mw_;
			X[P1_index_] = Y[P1_index_] * mw / MW_[P1_index_];
			X[P2_index_] = Y[P2_index_] * mw / MW_[P2_index_];
			X[P3_index_] = Y[P3_index_] * mw / MW_[P3_index_];
			X[P4_index_] = Y[P4_index_] * mw / MW_[P4_index_];

			const double mw_table_ = table_main_.Interpolation(Y[inert_index_], 0);
			X[I_index_] = Y[I_index_] * mw / mw_table_;
		}
		else if (iVersion_ == 171013)
		{
			X[fuel_index_] = Y[fuel_index_] * mw / fuel_mw_;
			X[oxidizer_index_] = Y[oxidizer_index_] * mw / oxidizer_mw_;
			X[inert_index_] = Y[inert_index_] * mw / inert_mw_;
			X[P1_index_] = Y[P1_index_] * mw / MW_[P1_index_];
			X[P2_index_] = Y[P2_index_] * mw / MW_[P2_index_];
			X[P3_index_] = Y[P3_index_] * mw / MW_[P3_index_];
			X[P4_index_] = Y[P4_index_] * mw / MW_[P4_index_];

			const double mw_table_ = table_main_.Interpolation(Y[inert_index_], 0);
			X[I_index_] = Y[I_index_] * mw / mw_table_;
		}

		if (iSubMechanism_CO_ == true)
		{
			X[CO_index_] = Y[CO_index_] * mw / MW_[CO_index_];
			X[V1_index_] = Y[V1_index_] * mw / MW_[V1_index_];
			X[V2_index_] = Y[V2_index_] * mw / MW_[V2_index_];
		}

		if (iSubMechanism_NO_ == true)
		{
			X[NO_index_] = Y[NO_index_] * mw / MW_[NO_index_];
			X[W1_index_] = Y[W1_index_] * mw / MW_[W1_index_];
			X[W2_index_] = Y[W2_index_] * mw / MW_[W2_index_];
			X[W3_index_] = Y[W3_index_] * mw / MW_[W3_index_];
		}
	}

	double VirtualChemistry::MW(const unsigned int j, double* Y)
	{
		if (iVersion_ == 170911)
		{
			if (j == I_index_)
			{
				const double mw_table_ = table_main_.Interpolation(Y[inert_index_], 0);
				MW_[I_index_] = mw_table_;
			}
		}
		else if (iVersion_ == 171013)
		{
			if (j == I_index_)
			{
				const double mw_table_ = table_main_.Interpolation(Y[inert_index_], 0);
				MW_[I_index_] = mw_table_;
			}
		}

		return MW_[j];
	}

	double VirtualChemistry::CpMix(const double T, const double P_Pa, double* Y)
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the specific heats of species in mass units [J/kg/K]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.cpMolar_Species(NSVector_.data());

		// Return the specific heat of the mixture [J/kg]
		double sum = 0.;
		for (unsigned int i=0;i<ns_main_;i++)
			sum += Y[i]*NSVector_(i);
		return sum;
	}

	void VirtualChemistry::CpSpecies(const double T, const double P_Pa, double* Cp)
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the specific heats of species in mass units [J/kg/K]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.cpMolar_Species(Cp);

		// Exclude sub-mechanisms species
		for (unsigned int i = ns_main_; i<ns_; i++)
			Cp[i] = 0.;
	}

	double VirtualChemistry::HMix(const double T, const double P_Pa, double* Y)
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the enthalpies of species in mass units [J/kg]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.hMolar_Species(NSVector_.data());

		// Return the enthalpy of the mixture [J/kg]
		double sum = 0.;
		for (unsigned int i=0;i<ns_main_;i++)
			sum += Y[i]*NSVector_(i);
		return sum;
	}

	double VirtualChemistry::Qdot(const double T, const double P_Pa, double* Omega) 
	{
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the enthalpies of species in mass units [J/kg]
		thermodynamicsMap_.SetPressure(P_Pa);
		thermodynamicsMap_.SetTemperature(T);
		thermodynamicsMap_.hMolar_Species(NSVector_.data());

		// Returns the reaction heat [J/m3/s]
		double sum = 0.;
		for (unsigned int i=0;i<ns_main_;i++)
			sum -= Omega[i]*NSVector_(i);
		return sum;
	}

	double VirtualChemistry::DynamicViscosity(const double T)
	{
		return mu0_*std::pow(T/T0_, Beta0_);
	}

	double VirtualChemistry::ThermalConductivity(const double mu, const double cp)
	{
		return mu*cp/Pr0_;
	}

	void VirtualChemistry::MassDiffusionCoefficients(const double lambda, const double rho, const double cp, double* Dmix)
	{
		const double alpha = lambda/rho/cp;
		for (unsigned int i=0;i<ns_;i++)
			Dmix[i] = alpha/Le_(i);
	}

	void VirtualChemistry::FormationRates(const double cTot, const double MW, const double T, double* Y, double* Omega)
	{
		// Mass fraction of N2
		const double YN2 = Y[inert_index_];
		// Lookup-table
		double f1 = 1.;
		double f2 = 1.;
		if (on_the_fly_optimization_ != "main")
		{
			table_main_.Interpolation(YN2);
			f1 = table_main_.interpolated()(1);
			f2 = table_main_.interpolated()(2);

			alpha1_ = table_main_.interpolated()(3);
			alpha2_ = table_main_.interpolated()(4);
			alpha3_ = table_main_.interpolated()(5);
			alpha4_ = table_main_.interpolated()(6);
		}

		// Mass fractions
		const double YF = Y[fuel_index_];
		const double YOX = Y[oxidizer_index_];
		const double YI = Y[I_index_];

		// Concentrations [mol/cm3]
		const double CF  = cTot * YF*MW/1. / 1000.;
		const double COX = cTot * YOX*MW/1. / 1000.;
		const double CI  = cTot * YI*MW/1. / 1000.;

		double r1 = 0.;
		double r2 = 0.;
		
		if (iReactions_ == true)
		{
			if (iVersion_ == 170911)
			{
				// Kinetic constants [mol, cm3, s]
				const double k1 = A1_*f1*std::exp(-E1_/ 1.987 / T);
				const double k2 = A2_*std::exp(-E2_/ 1.987 / T);

				// Reaction rates [kmol/m3/s]
				r1 = k1 * POW(CF, nuF_1_)*POW(COX, nuOX_1_) * 1000.;
				r2 = k2 * POW(CI, nuI_2_*f2) * 1000.;

				// 0.2FUEL + 0.8OX => I 		1.5384796696859927E+18	0.	3.5075174514637045E+04
				// FORD /FUEL 1.7099829697856470/
				// FORD /OX   0.8686294702424591/

				// I => 0.25P1+0.25P2+0.25P3+0.25P4 	3.9225000838284247E+18	0.	8.5680603523860715E+04
				// FORD /I 2.3399221348980621/
			}
			else if (iVersion_ == 171013)
			{
				// Kinetic constants [mol, cm3, s]
				const double k1 = 4.99463851138e+18*f1*std::exp(-50292.9081962 / 1.987 / T);
				const double k2 = 5.62191853964e+18*std::exp(-111552.26047 / 1.987 / T);

				// Reaction rates [kmol/m3/s]
				r1 = k1 * POW(CF, 0.000241938688303)*POW(COX, 2.49626471112)*1000.;
				r2 = k2 * POW(CI, 1.83822628643*f2)*1000.;

				// 0.2FUEL + 0.80OX => I 	 		4.99463851138e+18    0.0    50292.9081962
				// FORD /FUEL 0.000241938688303/
				// FORD /OX 2.49626471112/

				// I => 0.25P1+0.25P2+0.25P3+0.25P4 	 5.62191853964e+18   0.0    111552.26047
				// FORD /I 1.83822628643/
			}
		}
		
		// Formation rates [kg/m3/s]
		// Since the molecular weights of species are assumed equal to 1
		// we are calculating the formation rates of species in mass units [kg/m3/s]
		Omega[fuel_index_] = -0.2*r1;			// FUEL
		Omega[oxidizer_index_] = -0.8*r1;		// OX
		Omega[I_index_] = r1-r2;				// I
		Omega[P1_index_] = alpha1_ *r2;			// P1
		Omega[P2_index_] = alpha2_ *r2;			// P2
		Omega[P3_index_] = alpha3_ *r2;			// P3
		Omega[P4_index_] = alpha4_ *r2;			// P4
		Omega[inert_index_] = 0.;				// N2

		if (iReactions_ == true && iSubMechanism_CO_ == true)
		{
			// Interpolation
			table_co_.Interpolation(YN2);
			const double f3 = table_co_.interpolated()(0);
			const double f4 = table_co_.interpolated()(3);
			const double f5f = table_co_.interpolated()(1);
			const double f5b = table_co_.interpolated()(2);
			const double Kc5 = table_co_.interpolated()(4);
			const double alpha = table_co_.interpolated()(5);

			// Mass fractions
			const double YCO = Y[CO_index_];
			const double YV1 = Y[V1_index_];
			const double YV2 = Y[V2_index_];

			// Concentrations [mol/cm3]
			const double CCO = cTot * YCO*MW / 1. / 1000.;
			const double CV1 = cTot * YV1*MW / 1. / 1000.;
			const double CV2 = cTot * YV2*MW / 1. / 1000.;

			// Kinetic constants
			const double k3  = 1.53847967E+18*f3*std::exp(-3.50751745E+04 / 1.987 / T);
			const double k4  = 7.87700000E+18*f4*std::exp(-7.47828901E+04 / 1.987 / T);
			const double k5f = 2.82300000E+17*f5f*std::exp(-3.79155585E+04 / 1.987 / T);
			const double k5b = 2.82300000E+17*f5b*std::exp(-3.79155585E+04 / 1.987 / T);

			// Reaction rates [kmol/m3/s]
			const double r3 = k3 * POW(CF, 1.70998297)*POW(COX, 0.86862947)*1000.;
			const double r4 = k4 * POW(CF, 1.17030023)*POW(CV1, 1.34422510)*1000.;
			const double r5f = k5f * POW(CCO, 3.70916749)*POW(CV2, -1.08960225)*1000.;
			const double r5b = k5b / Kc5 * POW(CCO, 2.70916749)*POW(CV2, -0.08960225)*1000.;
			const double r5 = r5f - r5b;

			// Formation rates [kg/m3/s]
			// Since the molecular weights of species are assumed equal to 1
			// we are calculating the formation rates of species in mass units [kg/m3/s]
			Omega[CO_index_] = alpha*r3 + r4 - r5;		// CO
			Omega[V1_index_] = (1.-alpha)*r3 - r4;		// V1
			Omega[V2_index_] = r5;						// V2
		}

		if (iReactions_ == true && iSubMechanism_NO_ == true)
		{
			double f3 = 1.;
			double f4 = 1.;
			double f5f = 1.;
			double f5b = 1.;
			double f6 = 1.;
			double f7 = 1.;

			if (on_the_fly_optimization_ != "NO")	// interpolate values from tables
			{
				// Interpolation (table 1)
				table_no_1_.Interpolation(YN2);
				f3 = table_no_1_.interpolated()(0);
				f4 = table_no_1_.interpolated()(3);
				f5f = table_no_1_.interpolated()(1);
				f5b = table_no_1_.interpolated()(2);
				f6 = table_no_1_.interpolated()(4);
				f7 = table_no_1_.interpolated()(5);

				// Interpolation (table 2)
				table_no_2_.Interpolation(YN2);
				Kc5_NO_ = table_no_2_.interpolated()(0);
				beta_NO_ = table_no_2_.interpolated()(1);
				alpha_NO_ = table_no_2_.interpolated()(2);
				gamma_NO_ = table_no_2_.interpolated()(3);

				// Interpolation (table 3)
				table_no_3_.Interpolation(YN2);
				E3_NO_ = table_no_3_.interpolated()(0);
				E4_NO_ = table_no_3_.interpolated()(1);
				E5_NO_ = table_no_3_.interpolated()(2);
				E6_NO_ = table_no_3_.interpolated()(3);
				E7_NO_ = table_no_3_.interpolated()(4);

				// Interpolation (table 4)
				table_no_4_.Interpolation(YN2);
				nuF_3_NO_ = table_no_4_.interpolated()(0);
				nuOX_3_NO_ = table_no_4_.interpolated()(1);
				nuW1_4_NO_ = table_no_4_.interpolated()(2);
				nuNO_5f_NO_ = table_no_4_.interpolated()(3);
				nuW2_5f_NO_ = table_no_4_.interpolated()(4);
				nuNO_5b_NO_ = table_no_4_.interpolated()(5);
				nuW2_5b_NO_ = table_no_4_.interpolated()(6);
				nuW3_6_NO_ = table_no_4_.interpolated()(7);
				nuW3_7_NO_ = table_no_4_.interpolated()(8);

				// Interpolation (table 5)
				table_no_5_.Interpolation(YN2);
				A3_NO_ = table_no_5_.interpolated()(0);
				A4_NO_ = table_no_5_.interpolated()(1);
				A5_NO_ = table_no_5_.interpolated()(2);
				A6_NO_ = table_no_5_.interpolated()(3);
				A7_NO_ = table_no_5_.interpolated()(4);
			}

			{
				// Mass fractions
				const double YNO = Y[NO_index_];
				const double YW1 = Y[W1_index_];
				const double YW2 = Y[W2_index_];
				const double YW3 = Y[W3_index_];

				// Concentrations [mol/cm3]
				const double CNO = cTot * YNO*MW / 1. / 1000.;
				const double CW1 = cTot * YW1*MW / 1. / 1000.;
				const double CW2 = cTot * YW2*MW / 1. / 1000.;
				const double CW3 = cTot * YW3*MW / 1. / 1000.;

				// Kinetic constants
				const double k3 = A3_NO_ * f3*std::exp(-E3_NO_ / 1.987 / T);
				const double k4 = A4_NO_ * f4*std::exp(-E4_NO_ / 1.987 / T);
				const double k5f = A5_NO_ * f5f*std::exp(-E5_NO_ / 1.987 / T);
				const double k5b = A5_NO_ * f5b*std::exp(-E5_NO_ / 1.987 / T);
				const double k6 = A6_NO_ * f6*std::exp(-E6_NO_ / 1.987 / T);
				const double k7 = A7_NO_ * f7*std::exp(-E7_NO_ / 1.987 / T);

				// Reaction rates [kmol/m3/s]
				const double r3 = k3 * POW(CF, nuF_3_NO_)*POW(COX, nuOX_3_NO_)*1000.;
				const double r4 = k4 * POW(CW1, nuW1_4_NO_)*1000.;
				const double r5f = k5f * POW(CNO, nuNO_5f_NO_)*POW(CW2, nuW2_5f_NO_)*1000.;
				const double r5b = k5b / Kc5_NO_ * POW(CNO, nuNO_5b_NO_)*POW(CW2, nuW2_5b_NO_)*1000.;
				const double r5 = r5f - r5b;
				const double r6 = k6 * POW(CW3, nuW3_6_NO_)*1000.;
				const double r7 = k7 * POW(CW3, nuW3_7_NO_)*1000.;

				// Formation rates [kg/m3/s]
				// Since the molecular weights of species are assumed equal to 1
				// we are calculating the formation rates of species in mass units [kg/m3/s]
				Omega[NO_index_] = beta_NO_ * r4 + r5 + r6;
				Omega[W1_index_] = alpha_NO_ * r3 - r4;
				Omega[W2_index_] = (1. - alpha_NO_)*r3 + (1. - beta_NO_)*r4 - r5 + r7;
				Omega[W3_index_] = gamma_NO_ * r3 - r6 - r7;
			}
		}
	}

	double POW(const double C_mol_cm3, const double lambda, const double conversion)
	{
		// conversion = 1.e3   : mol/cm3 -> 1e-3/1e-6=1e3  kmol/m3
		// conversion = 1.e-3  : mol/m3  -> 1e-3/1. = 1e-3 kmol/m3

		if (lambda >= 1.)
		{
			return std::pow(C_mol_cm3, lambda);
		}
		else
		{
			const double Cstar = 1.e-8;
			const double ALFA = 1.e-5;
			const double H = 1.50*std::log(ALFA / (1. - ALFA));
			const double K = 2.00*std::log((1. - ALFA) / ALFA) / Cstar;
			const double delta = 1.e9;

			const double eps = 1.e-16;
			const double C = std::max(C_mol_cm3 * conversion, eps);		// [kmol/m3]

			const double m = (std::tanh(K*C + H) + 1.) / 2.;	// transition is for 9.e-6 < C < 1.5e-5 [kmol/m3]
			const double gamma = m * std::pow(C + m / delta, lambda) + (1. - m)*std::pow(Cstar, lambda - 1.)*C; // [kmol/m3]^lambda
			
			return gamma/std::pow(conversion, lambda);					// [mol,cm3]
		}
	}
}

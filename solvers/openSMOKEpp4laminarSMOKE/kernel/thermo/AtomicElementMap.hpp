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

//#include "thermo/AtomicElementMap.h"

namespace OpenSMOKE
{
	const AtomicElement AtomicElementMap::AtomicWeights[AtomicElementMap::nElements] = 
	{
		AtomicElement("[name]",   "E",    0),      
		AtomicElement("Hydrogen", "H",    1.008000016212463),
		AtomicElement("Helium",   "He",   4.002999782562256),
		AtomicElement("[name]",   "Li",   6.93900),
		AtomicElement("[name]",   "Be",   9.01220),
		AtomicElement("[name]",   "B",   10.81100),
		AtomicElement("Carbon",   "C",   12.010999679565430),
		AtomicElement("Nitrogen", "N",   14.0069999694824),
		AtomicElement("Oxygen",   "O",   15.998999595642090),
		AtomicElement("[name]", "F",   18.99840),
		AtomicElement("[name]", "Ne",  20.18300),
		AtomicElement("[name]", "Na",  22.98980),
		AtomicElement("[name]", "Mg",  24.31200),
		AtomicElement("[name]", "Al",  26.98150),
		AtomicElement("[name]", "Si",  28.08600),
		AtomicElement("[name]", "P",   30.97380),
		AtomicElement("[name]", "S",   32.06),
		AtomicElement("[name]", "Cl",  35.45300),
		AtomicElement("[name]", "Ar",  39.948001861572270),
		AtomicElement("[name]", "K",   39.10200),
		AtomicElement("[name]", "Ca",  40.08000),
		AtomicElement("[name]", "Sc",  44.95600),
		AtomicElement("[name]", "Ti",  47.90000),
		AtomicElement("[name]", "V",   50.94200),
		AtomicElement("[name]", "Cr",  51.99600),
		AtomicElement("[name]", "Mn",  54.93800),
		AtomicElement("[name]", "Fe",  55.84700),
		AtomicElement("[name]", "Co",  58.93320),
		AtomicElement("[name]", "Ni",  58.71000),
		AtomicElement("[name]", "Cu",  63.54000),
		AtomicElement("[name]", "Zn",  65.37000),
		AtomicElement("[name]", "Ga",  69.72000),
		AtomicElement("[name]", "Ge",  72.59000),
		AtomicElement("[name]", "As",  74.92160),
		AtomicElement("[name]", "Se",  78.96000),
		AtomicElement("[name]", "Br",  79.90090),
		AtomicElement("[name]", "Kr",  83.80000),
		AtomicElement("[name]", "Rb",  85.47000),
		AtomicElement("[name]", "Sr",  87.62000),
		AtomicElement("[name]", "Y",   88.90500),
		AtomicElement("[name]", "Zr",  91.22000),
		AtomicElement("[name]", "Nb",  92.90600),
		AtomicElement("[name]", "Mo",  95.94000),
		AtomicElement("[name]", "Tc",  99.00000),
		AtomicElement("[name]", "Ru", 101.07000),
		AtomicElement("[name]", "Rh", 102.90500),
		AtomicElement("[name]", "Pd", 106.40000),
		AtomicElement("[name]", "Ag", 107.87000),
		AtomicElement("[name]", "Cd", 112.40000),
		AtomicElement("[name]", "In", 114.82000),
		AtomicElement("[name]", "Sn", 118.69000),
		AtomicElement("[name]", "Sb", 121.75000),
		AtomicElement("[name]", "Te", 127.60000),
		AtomicElement("[name]", "I",  126.90440),
		AtomicElement("[name]", "Xe", 131.30000),
		AtomicElement("[name]", "Cs", 132.90500),
		AtomicElement("[name]", "Ba", 137.34000),
		AtomicElement("[name]", "La", 138.91000),
		AtomicElement("[name]", "Ce", 140.12000),
		AtomicElement("[name]", "Pr", 140.90700),
		AtomicElement("[name]", "Nd", 144.24000),
		AtomicElement("[name]", "Pm", 145.00000),
		AtomicElement("[name]", "Sm", 150.35000),
		AtomicElement("[name]", "Eu", 151.96000),
		AtomicElement("[name]", "Gd", 157.25000),
		AtomicElement("[name]", "Tb", 158.92400),
		AtomicElement("[name]", "Dy", 162.50000),
		AtomicElement("[name]", "Ho", 164.93000),
		AtomicElement("[name]", "Er", 167.26000),
		AtomicElement("[name]", "Tm", 168.93400),
		AtomicElement("[name]", "Yb", 173.04000),
		AtomicElement("[name]", "Lu", 174.99700),
		AtomicElement("[name]", "Hf", 178.49000),
		AtomicElement("[name]", "Ta", 180.94800),
		AtomicElement("[name]", "W",  183.85000),
		AtomicElement("[name]", "Re", 186.20000),
		AtomicElement("[name]", "Os", 190.20000),
		AtomicElement("[name]", "Ir", 192.20000),
		AtomicElement("[name]", "Pt", 195.09000),
		AtomicElement("[name]", "Au", 196.96700),
		AtomicElement("[name]", "Hg", 200.59000),
		AtomicElement("[name]", "Tl", 204.37000),
		AtomicElement("[name]", "Pb", 207.19000),
		AtomicElement("[name]", "Bi", 208.98000),
		AtomicElement("[name]", "Po", 210.00000),
		AtomicElement("[name]", "At", 210.00000),
		AtomicElement("[name]", "Rn", 222.00000),
		AtomicElement("[name]", "Fr", 223.00000),
		AtomicElement("[name]", "Ra", 226.00000),
		AtomicElement("[name]", "Ac", 227.00000),
		AtomicElement("[name]", "Th", 232.03800),
		AtomicElement("[name]", "Pa", 231.00000),
		AtomicElement("[name]", "U",  238.03000),
		AtomicElement("[name]", "Np", 237.00000),
		AtomicElement("[name]", "Pu", 242.00000),
		AtomicElement("[name]", "Am", 243.00000),
		AtomicElement("[name]", "Cm", 247.00000),
		AtomicElement("[name]", "Bk", 249.00000),
		AtomicElement("[name]", "Cf", 251.00000),
		AtomicElement("[name]", "Es", 254.00000),
		AtomicElement("[name]", "Fm", 253.00000),
		AtomicElement("[name]", "D",    2.01410),
		AtomicElement("[name]", "e",    5.45e-4) 
	};


	AtomicElementMap::AtomicElementMap()
	{
		for (int i=0; i<nElements; i++)
		{
			insert(make_pair(AtomicWeights[i].name(), AtomicWeights[i].mw()));
		}
	}

	// * * * * * * * * * * * * * * * * Global data  * * * * * * * * * * * * * * //
	AtomicElementMap AtomicWeights;

}

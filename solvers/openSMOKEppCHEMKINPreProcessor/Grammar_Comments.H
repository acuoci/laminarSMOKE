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
|   This file is part of OpenSMOKE++ Suite.                               |
|                                                                         |
|   Copyright(C) 2014, 2013  Alberto Cuoci                                |
|   Source-code or binary products cannot be resold or distributed        |
|   Non-commercial use only                                               |
|   Cannot modify source-code for any purpose (cannot create              |
|   derivative works)                                                     |
|                                                                         |
\*-----------------------------------------------------------------------*/

#include "dictionary/OpenSMOKE_DictionaryManager.h"
#include "dictionary/OpenSMOKE_DictionaryGrammar.h"
#include "dictionary/OpenSMOKE_DictionaryKeyWord.h"

class GrammarComments : public OpenSMOKE::OpenSMOKE_DictionaryGrammar
{
protected:

	virtual void DefineRules()
	{
		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Author", 
															OpenSMOKE::SEQUENCE_STRING, 
															"Name of the author of the pre-processing operations", 
															true) );
		
		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Place", 
															OpenSMOKE::SEQUENCE_STRING, 
															"Name of the location where the kinetic scheme was pre-processed", 
															true) );

		AddKeyWord( OpenSMOKE::OpenSMOKE_DictionaryKeyWord("@Comments", 
															OpenSMOKE::SEQUENCE_STRING, 
															"Comments about the kinetic scheme", 
															true) );
	}
};
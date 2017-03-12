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


#define COMPLETE_ODESOLVERINTERFACE_MEBDF(classname)\
bool classname::instanceFlag = false;\
classname* classname::m_this = NULL;

#define DEFINE_ODESOLVERINTERFACE_MEBDF(classname)\
private:\
	\
	static bool instanceFlag;\
	static classname* m_this;\
	classname()\
	{\
	}\
	\
public:\
	   \
	static classname* GetInstance()\
	{\
		if(instanceFlag == false)\
		{\
			m_this = new classname();\
			instanceFlag = true;\
			return m_this;\
		}\
		else\
		{\
			return m_this;\
		}\
	}\
	\
    	~classname()\
	{\
		instanceFlag = false;\
		std::cout << "~classname..." << std::endl; \
	}\
	\
	static void GetSystemFunctionsStatic(int *n, double *x, double *y, double *fy, int *ipar, double *rpar)\
	{\
		m_this->GetSystemFunctionsCallBack(n,x,y,fy, ipar, rpar);\
	}\
	\
	static void GetAnalyticalJacobianStatic(double *x, double *y, double *pd, int *n, int *meband, int *ipar, double *rpar)\
	{\
		m_this->GetAnalyticalJacobianCallBack(x, y, pd, n, meband, ipar, rpar);\
	}\
	\
	static void GetWriteFunctionStatic(double *x, double *y)\
	{\
		m_this->GetWriteFunctionCallBack(x,y);\
	}\
	\
        static void GetMassMatrixMEBDFStatic(int *n, double *am, int *masbnd, int *ipar, double *rpar)\
        {\
        }\
        \
	void GetSystemFunctionsCallBack(int *n, double *x, double *y, double *fy, int *ipar, double *rpar)\
	{\
		GetSystemFunctions(*x, y, fy);\
	}\
	\
	void GetAnalyticalJacobianCallBack(double *x, double *y, double *pd, int *n, int *meband, int *ipar, double *rpar)\
	{\
		GetAnalyticalJacobian(*x, y, pd);\
	}\
	\
	void GetWriteFunctionCallBack(double *x, double *y)\
	{\
		GetWriteFunction(*x, y);\
	}\
	\



     
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


#define COMPLETE_ODESOLVERINTERFACE_DASPK(classname)\
bool classname::instanceFlag = false;\
classname* classname::m_this = NULL; \
unsigned int classname::n_; \
double* classname::fy_;

#define DEFINE_ODESOLVERINTERFACE_DASPK(classname)\
private:\
	\
	static bool instanceFlag;\
	static classname* m_this;\
	static unsigned int n_; \
	static double* fy_; \
	\
	classname()\
	{\
	}\
	\
public:\
	   \
	static classname* GetInstance(const unsigned int n)\
	{\
		if(instanceFlag == false)\
		{\
			m_this = new classname();\
			instanceFlag = true;\
			m_this->n_ = n; \
			m_this->fy_ = new double[m_this->n_]; \
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
	static void GetSystemFunctionsStatic(double *x, double *y, double *dy, double *cj, double *delta, int *ires, double *rpar, int *ipar)\
	{\
		m_this->GetSystemFunctionsCallBack(x, y, dy, cj, delta, ires, rpar, ipar);\
	}\
	\
	static void GetSystemDerivativesStatic(double *x, double *y, double *dy)\
	{\
		m_this->GetSystemDerivativesCallBack(x, y, dy);\
	}\
	\
	static void GetAnalyticalJacobianStatic(double *x, double *y, double *dy, double *pd, double *cj, double *rpar, int *ipar)\
	{\
		m_this->GetAnalyticalJacobianCallBack(x, y, dy, pd, cj, rpar, ipar);\
	}\
	\
	static void GetKrylovSolverStatic(int *n, double *x, double *y, double *dy, double *savr, double *wk, double *cj, double *wght, double *wp, int *iwp, double *b, double *eplin, int *ier, double *rpar, int *ipar)\
	{\
		m_this->GetKrylovSolverCallBack(n, x, y, dy, savr, wk, cj, wght, wp, iwp, b, eplin, ier, rpar, ipar);\
	}\
	\
	static void GetWriteFunctionStatic(double *x, double *y)\
	{\
		m_this->GetWriteFunctionCallBack(x,y);\
	}\
	\
	void GetSystemFunctionsCallBack(double *x, double *y, double *dy, double *cj, double *delta, int *ires, double *rpar, int *ipar)\
	{\
		GetSystemFunctions(*x, y, m_this->fy_);\
		for(unsigned int i=0;i<m_this->n_;++i) \
			delta[i] = dy[i] - m_this->fy_[i]; \
	}\
	\
	void GetSystemDerivativesCallBack(double *x, double *y, double *dy)\
	{\
		GetSystemFunctions(*x, y, dy);\
	}\
	\
	void GetAnalyticalJacobianCallBack(double *x, double *y, double *dy, double *pd, double *cj, double *rpar, int *ipar)\
	{\
		GetAnalyticalJacobian(*x, y, pd);\
	}\
	\
	void GetKrylovSolverCallBack(int *n, double *x, double *y, double *dy, double *savr, double *wk, double *cj, double *wght, double *wp, int *iwp, double *b, double *eplin, int *ier, double *rpar, int *ipar)\
	{\
	}\
	\
	void GetWriteFunctionCallBack(double *x, double *y)\
	{\
		GetWriteFunction(*x, y);\
	}\
	\


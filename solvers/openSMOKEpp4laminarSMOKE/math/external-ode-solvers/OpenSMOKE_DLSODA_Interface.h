/***************************************************************************
 *   Copyright (C) 2012 by Alberto Cuoci								   *
 *   alberto.cuoci@polimi.it                                               *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#define COMPLETE_ODESOLVERINTERFACE_DLSODA(classname)\
bool classname::instanceFlag = false;\
classname* classname::m_this = NULL;

#define DEFINE_ODESOLVERINTERFACE_DLSODA(classname)\
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
	}\
	\
	static void GetSystemFunctionsStatic(int *n, double *x, double *y, double *fy)\
	{\
		m_this->GetSystemFunctionsCallBack(n,x,y,fy);\
	}\
	\
	static void GetAnalyticalJacobianStatic(int *n, double *x, double *y, int *ml, int *mu, double *pd, int *nrowpd)\
	{\
		m_this->GetAnalyticalJacobianCallBack(n,x,y,ml,mu,pd,nrowpd);\
	}\
	\
	static void GetWriteFunctionStatic(double *x, double *y)\
	{\
		m_this->GetWriteFunctionCallBack(x,y);\
	}\
	\
	void GetSystemFunctionsCallBack(int *n, double *x, double *y, double *fy)\
	{\
		GetSystemFunctions(*x, y, fy);\
	}\
	\
	void GetAnalyticalJacobianCallBack(int *n, double *x, double *y, int *ml, int *mu, double *pd, int *nrowpd)\
	{\
		GetAnalyticalJacobian(*x, y, pd);\
	}\
	\
	void GetWriteFunctionCallBack(double *x, double *y)\
	{\
		GetWriteFunction(*x, y);\
	}\
	\



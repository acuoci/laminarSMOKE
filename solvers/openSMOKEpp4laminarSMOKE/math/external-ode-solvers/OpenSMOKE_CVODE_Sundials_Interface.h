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


#define COMPLETE_ODESOLVERINTERFACE_CVODE_Sundials(classname)\
bool classname::instanceFlag = false;\
classname* classname::m_this = NULL;

#define DEFINE_ODESOLVERINTERFACE_CVODE_Sundials(classname)\
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
	static int GetSystemFunctionsStatic(realtype t, N_Vector y, N_Vector ydot, void *user_data)\
	{\
		return m_this->GetSystemFunctionsCallBack(t,y,ydot,user_data);\
	}\
	\
	static int GetAnalyticalJacobianStatic(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\
	{\
		return m_this->GetAnalyticalJacobianCallBack(N,t,y,fy,J,user_data,tmp1, tmp2, tmp3);\
	}\
	\
	static int GetWriteFunctionStatic(realtype t, N_Vector y)\
	{\
		return m_this->GetWriteFunctionCallBack(t,y);\
	}\
	\
	int GetSystemFunctionsCallBack(realtype t, N_Vector y, N_Vector ydot, void *user_data)\
	{\
		double* ydata = N_VGetArrayPointer(y);\
		double* fdata = N_VGetArrayPointer(ydot);\
		int flag = GetSystemFunctions(t, ydata, fdata);\
		return(flag);\
	}\
	\
	int GetAnalyticalJacobianCallBack(int N, realtype t, N_Vector y, N_Vector fy, DlsMat J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3)\
	{\
		return(0);\
	}\
	\
	int GetWriteFunctionCallBack(realtype t, N_Vector y)\
	{\
		double* ydata = N_VGetArrayPointer(y);\
		int flag = GetWriteFunction(t, ydata);\
		return flag; \
	}\

/*
		for(int i=0;i<n_;i++) y_[i] = NV_Ith_S(y,i);\
		int flag = GetSystemFunctions(t, y_, f_);\
		for(int i=0;i<n_;i++) NV_Ith_S(ydot,i)=f_[i];\
		return(flag);\
		*/

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

#ifndef MethodAdamsMoulton_H
#define MethodAdamsMoulton_H

#include <Eigen/Dense>
#include "OdeSolverUtilities.h"

namespace OdeSMOKE
{
	//!  A class implementing the Gear methods to solve stiff ODE systems
	/*!
	The purpose of this class is to solve a stiff ODE system using the Gear methods.
	The class is based on the ODESystemKernel policy, which implements the functions to solve the linear systems
	requested by the Gear methods.
	*/

	template <typename ODESystemKernel>
	class MethodAdamsMoulton : public ODESystemKernel
	{
	public:

		/**
		*@brief Default constructor
		*/
		MethodAdamsMoulton();

		/**
		*@brief Default destructor
		*/
		~MethodAdamsMoulton();

		/**
		*@brief Summary
		*@param out output stream where to write the output
		*/
		void OdeMethodSummary(std::ostream& out);

	protected:

		/**
		*@brief Performs the initialization of the class (memory allocation)
		*/
		void MemoryAllocationMethod();

		/**
		*@brief Reset of counters and preparation of parameters (default values)
		*/
		void ResetMethod();

		/**
		*@brief Core function of the class: calculates the correction b to be applied to the vector v
		*@param t current value of independent variable
		*@param one_over_epsilon error vector: one_over_epsilon = (tolAbs + tolRel*abs(y))^(-1)
		*/
		unsigned int FindCorrection(const double t, const Eigen::VectorXd& one_over_epsilon);

		/**
		*@brief Function internallyused to estimate the convergence rate of the method and apply changes on the
		integration parameters (if needed): order and step size
		*/
		void ConvergenceRate();

		/**
		*@brief In case of failure of FindCorrection function, the class try to find a possible, alternative solution
		*@param tInMeshPoint ???
		*@param tStabilize ???
		*/
		void WhatToDoInCaseOfConvergenceFailure(const double tInMeshPoint, double& tStabilize);

		/**
		*@brief Check the constraints imposed by the user on minimum and maximum values
		*@param y the vector to be constrained
		*/
		void CheckConstraints(Eigen::VectorXd& y);

		/**
		*@brief Check if the current solution does not satisfy the imposed constraints
		*@param y the vector which is supposed to satisfy the contraints
		*/
		unsigned int CheckIllegalConstraints(const Eigen::VectorXd& y);

	protected:

		// Parameters of the Gear methods
		Eigen::VectorXd* r_;			//!< vector r characteristic of the method (there is a vector for each order which was coded, up to MAX_ORDER)
		double* Ep_;					//!< error constant Ep (see Buzzi-Ferraris, eq. 29.163 and 29.204)

		// Internal vectors
		Eigen::VectorXd* z_;			//!< z vectors z
		Eigen::VectorXd* v_;			//!< v vectors v
		Eigen::VectorXd b_;				//!< correction vector b
		Eigen::VectorXd deltab_;		//!< correction vector for the correction vector b
		Eigen::VectorXd y_;
		Eigen::VectorXd f_;

		// Auxiliary vectors
		Eigen::VectorXd va_;				//!< auxiliary vector
		Eigen::VectorXd vb_;				//!< auxiliary vector

		// Parameter to choose the next order
		double deltaAlfa1_;					//!< reduction of safety coefficient for the current order minus 1 (see Buzzi-Ferraris, eq. 29.187)
		double deltaAlfa3_;					//!< reduction of safety coefficient for the current order plus 1 (see Buzzi-Ferraris, eq. 29.189)
		double* alfa2_;						//!< safety coefficient for the current order (see Buzzi-Ferraris, eq. 29.188)

		// Policy about the order and the step size
		OdeHStatus odeHStatus_;				//!< status of the step size (decreased, constant, or decreased)
		OdeOrderStatus odeOrderStatus_;		//!< status of the order (decreased, constant, or decreased)

		// Step size
		double h_;					//!< current step
		double min_step_size_;		//!< minimum step allowed
		double hNextStep_;			//!< next step to be adopted
		double hScale_;				//!< current scaling factor

		// Order
		unsigned int maximum_order_;	//!< maximum order to be used during the integration
		unsigned int maxOrderUsed_;		//!< the maximum order used during the integration
		unsigned int p_;				//!< current order
		unsigned int orderInNextStep_;	//!< order to be adopted in the next step

		// Number of steps
		unsigned int numberOfSteps_;				//!< current number of steps
		unsigned int maxConvergenceIterations_;		//!< maximum number of iterations to reach the convergence in the Newton's method (can be changed by the user)

		// Cumulative counters
		unsigned int numberOfFunctionCalls_;					//!< total number of calls to the system of equations
		unsigned int numberOfDecreasedSteps_;					//!< how many times the step size was decreased
		unsigned int numberOfIncreasedSteps_;					//!< how many times the step size was increased
		unsigned int numberOfConvergenceFailuresForOrderMax_;

		// Local counters
		unsigned int iterErrorFailure_;							//!< how many times, in the current step, the error resulted too large
		unsigned int iterConvergence_;							//!< number of iterations performed during iterative procedure
		unsigned int iterConvergenceFailure_;					//!< how many times, in the current step, the convergence procedure to find the correction failed
		unsigned int iterOrder_;								//!< how many steps were performed using the same order
		unsigned int iterConvergenceRate_;						//!< ???

		// Status
		OdeConvergence convergenceStatus_;						//!< status of the convergence procedure (success or failure)
		OdeStatus status_;										//!< status of the ODE system

		// Default values
		static const unsigned int MIN_ORDER;					//!< minimum order of the method
		static const unsigned int MAX_CONVERGENCE_FAILURE;		//!< maximum number of allowed successive convergence error
		static const unsigned int MAX_ORDER;					//!< maximum order which was coded in the current implementation 

		// Constraints on minimum and maximum values
		bool min_constraints_;									//!< constraints on minimum values are enabled
		bool max_constraints_;									//!< constraints on maximum values are enabled
		Eigen::VectorXd min_values_;							//!< allowed minimum values 
		Eigen::VectorXd max_values_;							//!< allowed maximum values 

	private:

		// Default values
		static const unsigned int DEFAULT_MAX_CONVERGENCE_ITER;		//!< maximum number of iteration for the Newton's method
		static const double DELTA_ALFA1;							//!< reduction of safety coefficient for the current order minus 1 (see Buzzi-Ferraris, eq. 29.187)
		static const double ALFA2;									//!< reduction of safety coefficient for the current order (see Buzzi-Ferraris, eq. 29.188)
		static const double DELTA_ALFA3;							//!< reduction of safety coefficient for the current order plus 1 (see Buzzi-Ferraris, eq. 29.189)
	};
}

#include "MethodAdamsMoulton.hpp"

#endif
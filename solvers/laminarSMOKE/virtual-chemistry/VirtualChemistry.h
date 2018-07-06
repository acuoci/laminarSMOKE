#ifndef OpenSMOKE_VirtualChemistry
#define OpenSMOKE_VirtualChemistry

namespace OpenSMOKE
{
	//!  A class to manage Virtual Chemistry Approach
	/*!
	This class provides the tools to manage Virtual Chemistry
	*/

	class VirtualChemistry
	{
	public:

		/**
		*@brief Default constructor
		*/
		VirtualChemistry( OpenSMOKE::ThermodynamicsMap_CHEMKIN& thermodynamicsMap,
				  boost::filesystem::path file_name );

		void SetOptimization(const unsigned int value);
		void SetFuel(const std::string name, const double mw);
		void SetOxidizer(const std::string name, const double mw);
		void SetInert(const std::string name, const double mw);

		void SetTransportProperties(const double mu0, const double T0, const double Beta0_, const double Pr0_);
		void SetLewisNumber(const unsigned int i, const double Le);

		double MWMix(double* Y);
		double CpMix(const double T, const double P_Pa, double* Y);
		double HMix(const double T, const double P_Pa, double* Y);
		double Qdot(const double T, const double P_Pa, double* Omega);

		double DynamicViscosity(const double T);
		double ThermalConductivity(const double mu, const double cp);
		void MassDiffusionCoefficients(const double lambda, const double rho, const double Cp, double* Dmix);

		void FormationRates(const double cTot, const double MW, const double T, double* Y, double* Omega);


		/**
		*@brief Interpolation
		*/
		unsigned int NumberOfTabulatedVariables() const { return nv_; };

		/**
		*@brief Interpolation
		*/
		void Interpolation(const double x, double* v);

		double InterpolationSingleValue(const double x, const unsigned int k);
		
                
	private:

		unsigned int nv_;
		unsigned int np_;

		double min_x_;
		double max_x_;

		Eigen::VectorXd	x_;
		Eigen::MatrixXd v_;
		Eigen::MatrixXd ratios_;
		Eigen::VectorXd table_;


		OpenSMOKE::ThermodynamicsMap_CHEMKIN&	thermodynamicsMap_;		//!< thermodynamic map

		std::string fuel_name_;
		std::string oxidizer_name_;
		std::string inert_name_;

		double fuel_mw_;
		double oxidizer_mw_;
		double inert_mw_;

		unsigned int fuel_index_;
		unsigned int oxidizer_index_;
		unsigned int inert_index_;
		unsigned int I_index_;
		unsigned int P1_index_;
		unsigned int P2_index_;
		unsigned int P3_index_;
		unsigned int P4_index_;

		unsigned int iOptimization_;

		unsigned int ns_;
		Eigen::VectorXd NSVector_;

		double mu0_;	// [kg/m/s]
		double T0_;	// [K]
		double Beta0_;
		double Pr0_;
		Eigen::VectorXd Le_;
	};
}

#include "VirtualChemistry.hpp"

#endif /* OpenSMOKE_VirtualChemistry */

#ifndef PHYSICAL_CONSTANTS_H
#define	PHYSICAL_CONSTANTS_H

#include <boost/math/constants/constants.hpp> 
#include <boost/math/special_functions/pow.hpp> 

namespace PhysicalConstants
{
	const double pi          		= boost::math::constants::pi<double>();			// pi
	const double pi_times_2    		= 2.*boost::math::constants::pi<double>();		// 2pi
	const double pi_over_2			= 0.5*boost::math::constants::pi<double>();		// pi/2
	const double pi_over_4			= 0.25*boost::math::constants::pi<double>();	// pi/4
	const double ln_10				= std::log(10.);										// ln(10)

	// Boltzmann constant
	const double kBoltzmann  		= 1.3806504e-23;    // Boltzmann constant [J/K]

	// Gas constant
	const double R_J_mol    		= 8.3144621;			// Ideal gas constant [J/mol/K]
	const double R_cal_mol			= 8.3144621/4.18443;	// Ideal gas constant [cal/mol/K]
	const double R_kcal_mol			= 8.3144621/4184.43;	// Ideal gas constant [kcal/mol/K]
	const double R_eV_mol			= 5.189479288e19;		// Ideal gas constant [eV/mol/K]
	const double R_J_kmol    		= 8314.4621;			// Ideal gas constant [J/kmol/K]
	const double R_cal_kmol    		= 8314.4621/4.18443;	// Ideal gas constant [cal/kmol/K]

	// Avogadro number
    const double Nav_mol         	= 6.0221417930e23;		// Avogadro constant [1/mol]
	const double Nav_kmol      		= 6.0221417930e26;  // Avogadro constant [1/kmol]

/*
	// Conversions // Moved to Conversions.h
    const double J_from_cal		= R_J_mol/R_cal_mol;
    const double J_from_eV		= R_J_mol/R_eV_mol;
    const double J_from_kcal	= 1.e3*J_from_cal;
*/
	// Constants
	const double sigmaStefanBoltzmann = 5.670373e-8;		// Stefan-Boltzmann constant [W/m2/K4]
    

	enum TAG_REACTION { REACTION_SIMPLE, REACTION_THIRDBODY, REACTION_LINDEMANN_FALLOFF, REACTION_SRI_FALLOFF, REACTION_TROE_FALLOFF, REACTION_LINDEMANN_CABR, REACTION_SRI_CABR, REACTION_TROE_CABR, REACTION_CHEBYSHEV };
	enum TAG_REACTION_SURFACE { REACTION_SURFACE_SIMPLE, REACTION_SURFACE_STICK, REACTION_SURFACE_LANGMUIR, REACTION_SURFACE_COVERAGE_DEPENDENT, REACTION_SURFACE_STICK_COVERAGE_DEPENDENT, REACTION_SURFACE_LUMPED };
	
	enum UNITS_REACTION { UNITS_CAL_MOLE, UNITS_EVOLTS, UNITS_JOULES_MOLE, UNITS_KCAL_MOLE, UNITS_KELVINS, UNITS_KJOULES_MOLE, UNITS_MOLECULES, UNITS_MOLES };
	enum UNITS_REACTION_COMPOSITION { UNITS_STD, UNITS_ATM, UNITS_BAR, UNITS_DYNES, UNITS_PASCALS, UNITS_TORR };
	
	enum UBIQEP_TYPE { UBIQEP_TYPE_ADSORPTION, UBIQEP_TYPE_DESORPTION, UBIQEP_TYPE_SURFACE, UBIQEP_TYPE_DUMMY };
	enum sensitivity_type { SENSITIVITY_FREQUENCY_FACTOR, SENSITIVITY_ACTIVATION_ENERGY, SENSITIVITY_KINETIC_CONSTANT };

	const double MAX_MW_THERMALDIFFUSION_RATIOS = 20;

	enum OpenSMOKE_GasMixture_Viscosity_Model { OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_WILKE, OPENSMOKE_GASMIXTURE_VISCOSITYMODEL_HERNING};
}

namespace Conversions
{
	const double from_atm_to_Pa = 101325.;
	const double from_bar_to_Pa = 100000.;
}

#endif	/* PHYSICAL_CONSTANTS_H */


/*
   
    const double Nav_mol         	= 6.02214129e23;	// Avogadro constant [1/mol]
    const double Nav_kmol      		= 6.02214129e26;	// Avogadro constant [1/kmol]
	const double hPlanck     		= 6.62606957e-34;	// Planck constant [J.s]
	
	const double g           		= 9.80665;			// gravitational acceleration [m/s2]
	const double m_u         		= 1.660538921e-27;  // atomic mass constant [kg]

*/


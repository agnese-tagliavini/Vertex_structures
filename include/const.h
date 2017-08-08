
/******************************************************************************************//** @file
 *  		
 * 	file: 		const.h
 * 	contents:  	Holds relevant constants for the grids and the numerical calculations
 * 
 ****************************************************************************************************/


#pragma once

#include <complex>
#include <string>

const std::complex<double> I( 0.0, 1.0 );			///< Imaginary unit
const double PI = 3.14159265358979323846; 			///< PI
const double LN_10 = 2.30258509299;				///< Natural log of 10

const int REDU = 0; 
const int COUNT = 50-REDU; 

// ----- SE dimensions

const int POS_FFREQ_COUNT_SIG = 200;				///< Amount of positive frequencies in self-energy grid
const int FFREQ_COUNT_SIG = 2 * POS_FFREQ_COUNT_SIG;		///< Amount of frequencies in self-energy grid

// ----- phi dimensions

const int POS_FFREQ_COUNT_PHI = COUNT;				///< Amount of positive fermionic frequencies in phi grid 
const int FFREQ_COUNT_PHI = 2 * POS_FFREQ_COUNT_PHI;		///< Amount of fermionic frequencies in phi grid

const int POS_BFREQ_COUNT_PHI = 75-REDU;			 	///< Amount of positive bosonic frequencies in phi grid 
const int BFREQ_COUNT_PHI = 2 * POS_BFREQ_COUNT_PHI + 1;	///< Amount of bosonic frequencies in phi grid

// ----- P dimensions

#ifdef SELFCONSISTENCY
const int POS_FFREQ_COUNT_P = 2 * COUNT; 
#elif defined ONESHOT
const int POS_FFREQ_COUNT_P = 25;
#else
const int POS_FFREQ_COUNT_P = 120;
#endif
const int FFREQ_COUNT_P = 2 * POS_FFREQ_COUNT_P;		///< Amount of fermionic frequencies in P grid                                                                                   
#ifdef SELFCONSISTENCY
const int POS_BFREQ_COUNT_P = 3 * POS_FFREQ_COUNT_P / 2;
#elif defined ONESHOT
const int POS_BFREQ_COUNT_P = POS_FFREQ_COUNT_P / 2;
#else
const int POS_BFREQ_COUNT_P = 180;
#endif
const int BFREQ_COUNT_P = 2 * POS_BFREQ_COUNT_P + 1;		///< Amount of bosonic frequencies in P grid

// ----- internal integration range and green function grid
#ifdef SELFCONSISTENCY 
const int POS_INT_RANGE = 2 * FFREQ_COUNT_SIG;
#elif defined ONESHOT
const int POS_INT_RANGE = COUNT;
#else
const int POS_INT_RANGE = 2 * FFREQ_COUNT_SIG;		///< Positive range for internal integrations
#endif
const int TAIL_LENGTH = POS_INT_RANGE/8;
const int FIT_ORDER = 4; 					///< Fit tail function has exponents one lower than this constant

// ----- chi dimensions
#ifdef SELFCONSISTENCY 
const int POS_BFREQ_COUNT_CHI =  POS_INT_RANGE;
#elif defined ONESHOT
const int POS_BFREQ_COUNT_CHI = POS_BFREQ_COUNT_PHI/2;
#else 
const int POS_BFREQ_COUNT_CHI = 1000; 		///< Amount of positive bosonic frequencies in chi grid
#endif
const int BFREQ_COUNT_CHI = 2 * POS_BFREQ_COUNT_CHI + 1;	///< Amount of bosonic frequencies in chi grid

// ----- Green function and single scale propagator dimension

const int POS_1P_RANGE = POS_BFREQ_COUNT_CHI + POS_INT_RANGE; 	///< Positive range for Green functions and Single scale propagator

// ----- Output ranges

const int POS_PLOT_RANGE_PHI = 2* POS_FFREQ_COUNT_PHI; 	///< Amount of positive frequencies in phi output grid
const int POS_PLOT_RANGE_VERT = 2*POS_PLOT_RANGE_PHI; 		///< Amount of positive frequencies in vertex output grid

//------ Inversion ranges Bethe-Salpeter equations

//const int  POS_INV_RANGE = COUNT; 
const int  POS_ASY_RANGE = 3 * POS_INV_RANGE;	 

const int TAIL_LENGTH_ASY = POS_ASY_RANGE/8;
const int FIT_ORDER_ASY = 4; 					///< Fit tail function has exponents one lower than this constant

#ifdef NO_MOMENTA
const int PATCH_COUNT = 1;					///< Amount of k-patches
#else
const int PATCH_COUNT = 8;					///< Amount of k-patches
#endif

const int QN_COUNT = 1;						///< Amount of possible tuples of the discrete quantum numbers


// ODE SOlVER

const double ERR_ABS = 1e-8;		///< Absolute error for ODE routines
const double ERR_REL = 1e-6;		///< Relative error for ODE routines
const double MAX_COUPLING = 1e+3; 	///< If this coupling is exceeded in the vertex function the ODE solver will abort

#define 	ERR_STEPPER 		runge_kutta_cash_karp54
const std::string ERR_STEPPER_STRING( "runge_kutta_cash_karp54" ); 


/*********************  INTERACTION FLOW ********************/

#ifdef INT_FLOW

const std::string FLOW_SCHEME_STRING( "Interaction flow" ); 
const std::string FLOW_SCHEME_ABBREV( "INTFL" ); 

const double LAM_START = 0.0;		///< Starting value of the flow parameter
const double LAM_FIN = 1.0;		///< Final value of the flow parameter
const double INIT_STEP = 0.25;		///< Initial step of the ODE solver
const double Lam = 1.0;
/*********************  EBERLEIN FLOW ********************/

#elif EBERL_FLOW

const std::string FLOW_SCHEME_STRING( "Eberlein flow" ); 
const std::string FLOW_SCHEME_ABBREV( "EBFL" ); 

const double LAM_START = 10.0;		///< Starting value of the flow parameter t -> Flow starts at energy scale 10^LAM_START
const double LAM_FIN = -10.0;		///< Final value of the flow parameter t -> Flow endsd at energy scale 10^LAM_FIN
const double INIT_STEP = -2.0;		///< Initial step of the ODE solver for the flow parameter t

/*********************  OMEGA FLOW ********************/

#elif OMEGA_FLOW

const std::string FLOW_SCHEME_STRING( "Omega flow" ); 
const std::string FLOW_SCHEME_ABBREV( "OMFL" ); 

const double LAM_START = 10.0;		///< Starting value of the flow parameter t -> Flow starts at energy scale 10^LAM_START
const double LAM_FIN = -10.0;		///< Final value of the flow parameter t -> Flow endsd at energy scale 10^LAM_FIN
const double INIT_STEP = -2.0;		///< Initial step of the ODE solver for the flow parameter t

/*********************  RESERVOIR FLOW ********************/

#elif RES_FLOW

const std::string FLOW_SCHEME_STRING( "Reservoir flow" ); 
const std::string FLOW_SCHEME_ABBREV( "RESFL" ); 

const double LAM_START = 10.0;		///< Starting value of the flow parameter t -> Flow starts at energy scale 10^LAM_START
const double LAM_FIN = -10.0;		///< Final value of the flow parameter t -> Flow endsd at energy scale 10^LAM_FIN
const double INIT_STEP = -2.0;		///< Initial step of the ODE solver for the flow parameter t

#endif


/******************************************************************************************//** @file
 *  		
 * 	file: 		params.h
 * 	contents:  	Declare global parameters as extern
 * 
 ****************************************************************************************************/


#pragma once

#include <vector>

/********************* PHYSICAL PARAMETERS  ********************/	

// -------- System independent 

extern double UINT;	///< Bare interaction
extern double BETA; 	///< Inverse Temperature 
extern double B;	///< Zeeman splitting


// --------- Give Bath description
#ifdef 	QMC_BATH
#define	BATH_DOS_STRING	"Semi-Circular bath DOS"
#elif 	ED_BATH
#define	BATH_DOS_STRING	"Discrete bath DOS according to ED parameters"
#else
#define	BATH_DOS_STRING	"Structureless bath DOS ( Wide-band limit )"
#endif

// --------- Specify ED Bath parameters

extern std::vector<double> energies; 		///< ED Bath energies
extern std::vector<double> hybridizations; 	///< ED Bath hybridizations


// -------- SIAM / SQDJJ, compile with -D NO_MOMENTA

	// INDEPENDENT
extern double GAM_L;	///< Constant hybridization of left lead, CAUTION: lead hybridization asymmetry currently not implemented
extern double DEL;	///< Superconducting gap
extern double EPS;	///< Level position, shifted such that 0 corresponds to particle-hole symmetric case
extern double PHI;	///< Phase difference between left and right lead in units of PI

	// DEPENDENT
extern double PHI_L;	///< Phase of left superconducting order parameter
extern double PHI_R;	///< Phase of right superconducting order parameter 
extern double GAM_R;	///< Constant hybridization of right lead, CAUTION: lead hybridization assymetry currently not implemented
extern double DD;	///< Constant appearing naturally in atomic limit

void update_dep_params(); ///< Function to allow the update of all dependent parameters after a change of the independent ones

// -------- HUBBARD MODEL, compile without -D NO_MOMENTA

extern double MU;	///< Chemical potential
extern double T_PRIME;	///< Next nearest neighbour hopping

/********************* Parameter list for read in ********************/

extern std::vector<double*> readIn_lst; 

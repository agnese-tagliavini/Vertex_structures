
/************************************************************************************************//**
 *  		
 * 	file: 		params.o
 * 	contents:   	Initialize system parameters
 * 
 ****************************************************************************************************/


#include <params.h>
#include <math.h>


/********************* PHYSICAL PARAMETERS  ********************/	

// -------- System independent 

double UINT	=	1.0;	///< Bare interaction
double BETA 	=	50.0;	///< Inverse Temperature 
double B	=	0.0;	///< Zeeman splitting


// --------- Specify ED Bath parameters
//std::vector<double> energies = { -0.32, 0.32 } ;               	// BATH Paper ( optimized for Beta=20 ) to fit box DOS
//std::vector<double> energies = { -0.39, 0.39 } ;               	// BATH Paper ( optimized for Beta=10 ) to fit box DOS
//std::vector<double> hybridizations = { 0.56, 0.56 } ;
std::vector<double> energies = { -0.7, -0.15, 0.15, 0.7 } ;  		// BATH Paper ( optimized for Beta=20 ) to fit box DOS
std::vector<double> hybridizations = { 0.45, 0.34, 0.34, 0.45} ;                   

// -------- SIAM / SQDJJ, compile with -D NO_MOMENTA

	// INDEPENDENT
double GAM_L	= 	0.5; 	///< Constant hybridization of left lead, CAUTION: lead hybridization asymmetry currently not implemented
double DEL	=	0.0;	///< Superconducting gap
double EPS	= 	0.0; 	///< Level position, shifted such that 0 corresponds to particle-hole symmetric case
double PHI	=	0.0;	///< Phase difference between left and right lead in units of PI

	// DEPENDENT
double PHI_L	= PHI/2.0;	///< Phase of left superconducting order parameter
double PHI_R	= -PHI_L;	///< Phase of right superconducting order parameter 
double GAM_R	=1.0 - GAM_R;	///< Constant hybridization of right lead, CAUTION: lead hybridization assymetry currently not implemented
double DD	=cos( PHI_L );	///< Constant appearing naturally in atomic limit

void update_dep_params()
{
   PHI_L	= PHI/2.0;	
   PHI_R	= -PHI_L;	
   GAM_R	=1.0 - GAM_R;	
   DD		=cos( PHI_L );	
}

// -------- HUBBARD MODEL, compile without -D NO_MOMENTA

double MU	=	0.0;	///< Chemical potential
double T_PRIME	=	0.0;	///< Next nearest neighbour hopping


/********************* Parameter list for read in ********************/

std::vector<double*> readIn_lst = { &UINT, &BETA }; //, &EPS }; 



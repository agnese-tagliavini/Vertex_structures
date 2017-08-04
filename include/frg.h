
/*******************************************************************************************//** @file
 *  		
 * 	file: 		frg.h
 * 	contents:  	Definition of Green function ( i.e. system ) and fRG specific functions 
 * 
 ****************************************************************************************************/


#pragma once

#include <def.h>
#include <grid.h>
#include <params.h>

// --- Exact data arrays
extern gf_1p_t Sig_exact; 

extern gf_phi_t vert_exact_pp;
extern gf_phi_t vert_exact_ph; 
extern gf_phi_t vert_exact_xph;

extern gf_P_t P_exact_pp;
extern gf_P_t P_exact_ph; 
extern gf_P_t P_exact_xph;

extern gf_chi_t chi_exact_pp;
extern gf_chi_t chi_exact_ph; 
extern gf_chi_t chi_exact_xph;

extern gf_phi_t genchi_exact_pp;
extern gf_phi_t genchi_exact_ph; 
extern gf_phi_t genchi_exact_xph;

extern const int POS_FERM_VERT_COUNT_EXACT_SMALL, POS_BOS_VERT_COUNT_EXACT_SMALL;
extern const int POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT, POS_SIG_COUNT_EXACT; 
extern const int POS_FERM_P_COUNT_EXACT, POS_BOS_P_COUNT_EXACT, POS_BOS_CHI_COUNT_EXACT; 

void read_exact(); 

#ifdef NO_MOMENTA
MatQN Gam( double w );				///< Hybridization function according to the compile flags in makefile
#endif

MatQN G0inv( double w, double kx, double ky );					///< Inverse non-interacting Greens function ( scale-independent! )
MatQN G( double w, double kx, double ky, double Lam, const MatQN& selfEn );	///< Scale-dependent Greens function, introduce regulator here!
MatQN S( double w, double kx, double ky, double Lam, const MatQN& selfEn ); 	///< Single-scale propagator 

dcomplex vert_bare( int s1_in, int s2_in, int s1_out, int s2_out );  	///< Initial vertex values in Nambu basis
dcomplex vert_bare( const idx_2p_t& idx );				///< Return initial vertex value for given index object

dcomplex Sig_init( const idx_1p_t& idx );					///< Return initial self-energy value for given index object
dcomplex phi_init( const idx_phi_t& idx ); 

dcomplex P_pp_init( const idx_P_t& idx ); 
dcomplex P_ph_init( const idx_P_t& idx ); 
dcomplex P_xph_init( const idx_P_t& idx ); 

dcomplex chi_pp_init( const idx_chi_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex chi_ph_init( const idx_chi_t& idx );				///< Return initial vertex value for chiasch functions
dcomplex chi_xph_init( const idx_chi_t& idx );				///< Return initial vertex value for chiasch functions

dcomplex asympt_GG_pp( int W, double Lam = LAM_FIN );	///< Give estimate for the integral 1/2/PI * G(W/2-w-1,Lam) * G(W/2+w,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE )
dcomplex asympt_GG_ph( int W, double Lam = LAM_FIN );	///< Give estimate for the integral 1/2/PI * G(w-W/2,Lam) * G(w+W/2,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE )

dcomplex FUNC_PP( int W, int POS_INV );
dcomplex FUNC_PH( int W, int POS_INV );

inline double w_val( int w_idx )					///< Return value of fermionic matsubara frequency for a given index
{
   return PI / BETA * ( 2 * w_idx + 1 ); 
}

inline double W_val( int W_idx )					///< Return value of bosonic matsubara frequency for a given index
{
   return 2.0 * PI / BETA * W_idx; 
}

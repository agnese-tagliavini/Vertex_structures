
/*******************************************************************************************//** @file
 *  		
 * 	file: 		translate.h
 * 	contents:  	Functions that allow for translations between different notations
 * 
 ****************************************************************************************************/


#pragma once

#include <def.h>

/********************* Index translations P->phi, chi->phi, P->chi  ********************/

idx_phi_t iP_to_iphi( const idx_P_t& idx ); 
idx_P_t iphi_to_iP( const idx_phi_t& idx ); 
idx_phi_t ichi_to_iphi( const idx_chi_t& idx ); 
idx_chi_t iP_to_ichi( const idx_P_t& idx ); 


/********************* Translations between different notations  ********************/

struct FERM {}; 

template< typename notation = FERM >
idx_2p_t to_ferm( const idx_2p_t& idx )
{
   return idx; 
}

template< typename notation = FERM >
idx_2p_t from_ferm( const idx_2p_t& idx )
{
   return idx; 
}

template< typename to_notation = FERM, typename from_notation = FERM >
idx_2p_t translate( const idx_2p_t& idx )
{
   return from_ferm<to_notation>( to_ferm<from_notation>( idx ) ); 
}

#define DECLARE_NOTATION(X) 				\
   struct X {}; 					\
							\
template<>						\
idx_2p_t to_ferm<X>( const idx_2p_t& idx ); 		\
							\
template<>						\
idx_2p_t from_ferm<X>( const idx_2p_t& idx ); 	

//DECLARE_NOTATION(PP)
//DECLARE_NOTATION(PH)
//DECLARE_NOTATION(XPH)
DECLARE_NOTATION(PP_S)
DECLARE_NOTATION(PH_S)
DECLARE_NOTATION(XPH_S)
//DECLARE_NOTATION(BOS)

//idx_2p_t (*const iphi_pp_to_i2p)( const idx_2p_t& idx ) = to_ferm<PP_S>; 
//idx_2p_t (*const iphi_ph_to_i2p)( const idx_2p_t& idx ) = to_ferm<PH_S>; 
//idx_2p_t (*const iphi_xph_to_i2p)( const idx_2p_t& idx ) = to_ferm<XPH_S>; 

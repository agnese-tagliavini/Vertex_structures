
/*******************************************************************************************//** @file
 *  		
 * 	file: 		translate.cpp
 * 	contents:  	Functions that allow for translations between different notations
 * 
 ****************************************************************************************************/


#include <translate.h>
#include <grid.h>
#include <cmath>
#include <mymath.h>


idx_phi_t iP_to_iphi( const idx_P_t& idx )
{
   return  idx_phi_t( { idx( IP::W ), idx( IP::w ), 10*BFREQ_COUNT_CHI, idx( IP::K ), idx( IP::k ), 0, idx( IP::s1_in ), idx( IP::s2_in ), idx( IP::s1_out ), idx( IP::s2_out )} ); 
}

idx_P_t iphi_to_iP( const idx_phi_t& idx )
{
   return  idx_P_t( { idx( IPHI::W ), idx( IPHI::w_in ), idx( IPHI::K ), idx( IPHI::k_in ), idx( IPHI::s1_in ), idx( IPHI::s2_in ), idx( IPHI::s1_out ), idx( IPHI::s2_out )} ); 
}

idx_phi_t ichi_to_iphi( const idx_chi_t& idx )
{
   return  idx_phi_t( { idx( ICHI::W ), 10*BFREQ_COUNT_CHI, 10*BFREQ_COUNT_CHI, idx( ICHI::K ), 0, 0, idx( ICHI::s1_in ), idx( ICHI::s2_in ), idx( ICHI::s1_out ), idx( ICHI::s2_out )} );
}

idx_chi_t iP_to_ichi( const idx_P_t& idx )
{
   return  idx_chi_t( { idx( IP::W ), idx( IP::K ), idx( IP::s1_in ), idx( IP::s2_in ), idx( IP::s1_out ), idx( IP::s2_out )} ); 
}

// -- PP Shifted
template<>
idx_2p_t to_ferm<PP_S>( const idx_2p_t& idx_pp )
{
   int W = idx_pp( IPHI::W ); 
   int w_in = idx_pp( IPHI::w_in ); 
   int w_out = idx_pp( IPHI::w_out ); 

   int K = idx_pp( IPHI::K ); 
   int k_in = idx_pp( IPHI::k_in ); 
   int k_out = idx_pp( IPHI::k_out ); 

   idx_2p_t idx_2p( idx_pp.idx_arr ); 

   idx_2p( I2P::w1_in ) = w_in + div2_ceil( W ); 
   idx_2p( I2P::w2_in ) = div2_floor( W ) - w_in - 1;  // care for shift 
   idx_2p( I2P::w1_out ) = w_out + div2_ceil( W );  

   idx_2p( I2P::k1_in ) = k_in; 
   idx_2p( I2P::k2_in ) = dif_k(K, k_in);  
   idx_2p( I2P::k1_out ) = k_out;  

   return idx_2p; 
}

template<>
idx_2p_t from_ferm<PP_S>( const idx_2p_t& idx_2p )
{
   int w1_in = idx_2p( I2P::w1_in ); 
   int w2_in = idx_2p( I2P::w2_in ); 
   int w1_out = idx_2p( I2P::w1_out ); 

   int k1_in = idx_2p( I2P::k1_in ); 
   int k2_in = idx_2p( I2P::k2_in ); 
   int k1_out = idx_2p( I2P::k1_out ); 

   idx_phi_t idx_pp( idx_2p.idx_arr ); 

   idx_pp( IPHI::W ) = w1_in + w2_in + 1; 
   idx_pp( IPHI::w_in ) = w1_in - div2_ceil( w1_in + w2_in + 1 ); 
   idx_pp( IPHI::w_out ) = w1_out - div2_ceil( w1_in + w2_in + 1 );  

   idx_pp( IPHI::K ) = add_k( k2_in, k1_in ); ; 
   idx_pp( IPHI::k_in ) = k1_in;  
   idx_pp( IPHI::k_out ) = k1_out;  

   return idx_pp; 
}

// -- PH Shifted
template<>
idx_2p_t to_ferm<PH_S>( const idx_2p_t& idx_ph )
{
   int W = idx_ph( IPHI::W ); 
   int w_in = idx_ph( IPHI::w_in ); 
   int w_out = idx_ph( IPHI::w_out ); 

   int K = idx_ph( IPHI::K ); 
   int k_in = idx_ph( IPHI::k_in ); 
   int k_out = idx_ph( IPHI::k_out ); 

   idx_2p_t idx_2p( idx_ph.idx_arr ); 

   idx_2p( I2P::w1_in ) = w_in - div2_floor( W ); 
   idx_2p( I2P::w2_in ) = w_out + div2_ceil( W );  
   idx_2p( I2P::w1_out ) = w_in + div2_ceil( W );  

   idx_2p( I2P::k1_in ) = k_in; 
   idx_2p( I2P::k2_in ) = add_k(K,k_out);  
   idx_2p( I2P::k1_out ) = add_k(K,k_in);  

   return idx_2p; 
}

template<>
idx_2p_t from_ferm<PH_S>( const idx_2p_t& idx_2p )
{
   int w1_in = idx_2p( I2P::w1_in ); 
   int w2_in = idx_2p( I2P::w2_in ); 
   int w1_out = idx_2p( I2P::w1_out ); 

   int k1_in = idx_2p( I2P::k1_in ); 
   int k2_in = idx_2p( I2P::k2_in ); 
   int k1_out = idx_2p( I2P::k1_out ); 

   idx_phi_t idx_ph( idx_2p.idx_arr ); 

   idx_ph( IPHI::W ) = w1_out - w1_in; 
   idx_ph( IPHI::w_in ) = w1_in + div2_floor( w1_out - w1_in ); 
   idx_ph( IPHI::w_out ) = w2_in - div2_ceil( w1_out - w1_in );  

   idx_ph( IPHI::K ) = dif_k(k1_out, k1_in); 
   idx_ph( IPHI::k_in ) = k1_in;  
   idx_ph( IPHI::k_out ) = dif_k(k2_in, dif_k(k1_out, k1_in) );  

   return idx_ph; 
}

// -- XPH Shifted
template<>
idx_2p_t to_ferm<XPH_S>( const idx_2p_t& idx_xph )
{
   int W = idx_xph( IPHI::W ); 
   int w_in = idx_xph( IPHI::w_in ); 
   int w_out = idx_xph( IPHI::w_out ); 

   int K = idx_xph( IPHI::K ); 
   int k_in = idx_xph( IPHI::k_in ); 
   int k_out = idx_xph( IPHI::k_out ); 

   idx_2p_t idx_2p( idx_xph.idx_arr ); 

   idx_2p( I2P::w1_in ) = w_in - div2_floor( W ); 
   idx_2p( I2P::w2_in ) = w_out + div2_ceil( W );  
   idx_2p( I2P::w1_out ) = w_out - div2_floor( W );  

   idx_2p( I2P::k1_in ) = k_in; 
   idx_2p( I2P::k2_in ) = add_k(K,k_out);  
   idx_2p( I2P::k1_out ) = k_out;  

   return idx_2p; 
}

template<>
idx_2p_t from_ferm<XPH_S>( const idx_2p_t& idx_2p )
{
   int w1_in = idx_2p( I2P::w1_in ); 
   int w2_in = idx_2p( I2P::w2_in ); 
   int w1_out = idx_2p( I2P::w1_out ); 

   int k1_in = idx_2p( I2P::k1_in ); 
   int k2_in = idx_2p( I2P::k2_in ); 
   int k1_out = idx_2p( I2P::k1_out ); 

   idx_phi_t idx_xph( idx_2p.idx_arr ); 

   idx_xph( IPHI::W ) = w2_in - w1_out; 
   idx_xph( IPHI::w_in ) = w1_in + div2_floor( w2_in - w1_out ); 
   idx_xph( IPHI::w_out ) = w2_in - div2_ceil( w2_in - w1_out );  

   idx_xph( IPHI::K ) = dif_k(k2_in, k1_out); 
   idx_xph( IPHI::k_in ) = k1_in;  
   idx_xph( IPHI::k_out ) = k1_out;  

   return idx_xph; 
}

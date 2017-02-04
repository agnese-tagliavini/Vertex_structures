
/************************************************************************************************//**
 *  		
 * 	file: 		observables.cpp
 * 	contents:   	See observables.h
 * 
 ****************************************************************************************************/


#include <observables.h>
#include <complex>
#include <params.h>
#include <frg.h>
#include <mymath.h>
#include <grid.h>

using namespace std; 

flow_obs_func_lst_t  flow_obs_func_lst = { { "EFF_MASS", eff_mass }, { "ERR_EFF_MASS", err_eff_mass } }; 	///< List to contain all observable functions to be tracked during the flow

double eff_mass( const state_t& state_vec )
{
   return 1.0 - BETA/2.0/PI * ( imag( state_vec.Sig( 1, 0, 0, 0 ) - state_vec.Sig( 0, 0, 0, 0 ) ) ); 
}

double err_eff_mass( const state_t& state_vec )
{
   return eff_mass( state_vec ) - ( 1.0 - BETA/4.0/PI * ( imag( state_vec.Sig( 2, 0, 0, 0 ) - state_vec.Sig( 0, 0, 0, 0 ) ) ) ) ; 
}

double filling( int s_in, int s_out, const state_t& state_vec )
{
   dcomplex val( 0.0, 0.0 );

   gf_1p_mat_t Gvec( POS_1P_RANGE ); 
   Gvec.init( boost::bind( &state_t::GMat, boost::ref(state_vec), _1, LAM_FIN ) ); // Initialize big Green function vector 

   for( int w = -POS_1P_RANGE; w < POS_1P_RANGE; ++w )
      for( int k = 0; k < PATCH_COUNT; ++k )
	    val += Gvec[w][k]( s_out, s_in );
   
   val *= 1.0 / BETA / PATCH_COUNT; 

   if( s_in == s_out )
      val += 0.5; // Large frequency contribution
      
   if( abs( imag( val )) > 0.00001 && s_in == s_out ) 
      cout << " CAUTION, nonvanishing imaginary part of filling: " << imag( val ) << endl; 
		  
   return real( val ) ; 
}

double jos_curr( const state_t& state_vec )
{
   double val( 0.0 );

   gf_1p_mat_t Gvec( POS_1P_RANGE ); 
   Gvec.init( boost::bind( &state_t::GMat, boost::ref(state_vec), _1, LAM_FIN ) ); // Initialize big Green function vector 

#pragma omp parallel for schedule(dynamic) reduction(+:val)
   for( int w = -POS_1P_RANGE; w < POS_1P_RANGE; ++w ) 
      for( int k = 0; k < PATCH_COUNT; ++k )
      {
	 double w_val = ( 0.5 + w ) * 2*PI / BETA; 
	 val += 4.0 * GAM_L * imag( polar(DEL, PHI_L) / sqrt(w_val*w_val+DEL*DEL) * Gvec[w][k]( 1, 0 ) );
      }

   val *= 1.0 / BETA / PATCH_COUNT; 
		  
   return val; 
}

double get_max_cpl( const state_t& state_vec )
{
   const dcomplex* max_pp_ptr = max_element( state_vec.gf_phi_pp().data(), state_vec.gf_phi_pp().data() + state_vec.gf_phi_pp().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_pp_distance = distance( state_vec.gf_phi_pp().data(), max_pp_ptr ); 
   cout  << " Max Phi pp " << *max_pp_ptr << " at idx " << state_vec.gf_phi_pp().get_idx( max_pp_distance ) << endl; 

   const dcomplex* max_ph_ptr = max_element( state_vec.gf_phi_ph().data(), state_vec.gf_phi_ph().data() + state_vec.gf_phi_ph().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_ph_distance = distance( state_vec.gf_phi_ph().data(), max_ph_ptr ); 
   cout  << " Max Phi ph " << *max_ph_ptr << " at idx " << state_vec.gf_phi_ph().get_idx( max_ph_distance ) << endl; 

   const dcomplex* max_xph_ptr = max_element( state_vec.gf_phi_xph().data(), state_vec.gf_phi_xph().data() + state_vec.gf_phi_xph().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_xph_distance = distance( state_vec.gf_phi_xph().data(), max_xph_ptr ); 
   cout << " Max Phi xph " << *max_xph_ptr << " at idx " << state_vec.gf_phi_xph().get_idx( max_xph_distance ) << endl; 

   return max( max( abs(*max_pp_ptr), abs(*max_ph_ptr) ), abs(*max_xph_ptr) );
}

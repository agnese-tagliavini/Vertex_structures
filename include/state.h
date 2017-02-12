
/******************************************************************************************//** @file
 *  		
 * 	file: 		state.h
 * 	contents:  	Definition of ODE state_t class including its norm
 * 
 ****************************************************************************************************/


#pragma once 

#include <arithmetic_tuple.h>
#include <grid.h>
#include <boost/numeric/odeint.hpp>
#include <def.h>

class state_t: public ReaK::arithmetic_tuple< gf_1p_t, gf_phi_t, gf_phi_t, gf_phi_t, gf_P_t, gf_P_t, gf_P_t, gf_chi_t, gf_chi_t, gf_chi_t > ///< Type of the state vector of the ODE solver
{
   public:
      using base_t = arithmetic_tuple< gf_1p_t, gf_phi_t, gf_phi_t, gf_phi_t, gf_P_t, gf_P_t, gf_P_t, gf_chi_t, gf_chi_t, gf_chi_t >; 
      using base_base_t = base_t::base_t; 
      using Sig_t = gf_1p_t; 
      using phi_t = gf_phi_t; 
      using P_t = gf_P_t; 
      using chi_t = gf_chi_t; 

      inline Sig_t& gf_Sig()			{ return std::get<0>( *this ); }
      inline const Sig_t& gf_Sig() const	{ return std::get<0>( *this ); }

      inline phi_t& gf_phi_pp()			{ return std::get<1>( *this ); }
      inline const phi_t& gf_phi_pp() const	{ return std::get<1>( *this ); }

      inline phi_t& gf_phi_ph()			{ return std::get<2>( *this ); }
      inline const phi_t& gf_phi_ph() const	{ return std::get<2>( *this ); }

      inline phi_t& gf_phi_xph()		{ return std::get<3>( *this ); }
      inline const phi_t& gf_phi_xph() const	{ return std::get<3>( *this ); }

      inline P_t& gf_P_pp()			{ return std::get<4>( *this ); }
      inline const P_t& gf_P_pp() const		{ return std::get<4>( *this ); }

      inline P_t& gf_P_ph()			{ return std::get<5>( *this ); }
      inline const P_t& gf_P_ph() const		{ return std::get<5>( *this ); }

      inline P_t& gf_P_xph()			{ return std::get<6>( *this ); }
      inline const P_t& gf_P_xph() const	{ return std::get<6>( *this ); }

      inline chi_t& gf_chi_pp()			{ return std::get<7>( *this ); }
      inline const chi_t& gf_chi_pp() const	{ return std::get<7>( *this ); }

      inline chi_t& gf_chi_ph()			{ return std::get<8>( *this ); }
      inline const chi_t& gf_chi_ph() const	{ return std::get<8>( *this ); }

      inline chi_t& gf_chi_xph()		{ return std::get<9>( *this ); }
      inline const chi_t& gf_chi_xph() const	{ return std::get<9>( *this ); }

      K_Grid mom_grid;								///< Momentum grid

      state_t():
	 base_t(),
	 mom_grid()
   {}; 

      state_t( const state_t& state_vec ):
	 base_t( state_vec ), 
	 mom_grid()
   {}; 

      state_t( state_t&& state_vec ):
	 base_t( std::move(state_vec) ), 
	 mom_grid()
   {}; 

      state_t& operator=( const state_t& state_vec )
      {
	 base_t::operator=( state_vec );
      }

      state_t& operator=( state_t&& state_vec )
      {
	 base_t::operator=( std::move( state_vec ) );
      }

      /********************* Interfacing gf containers  ********************/

      dcomplex Sig( int w, int k, int s_in, int s_out ) const; 		///< Return self-energy value for a specific set of index

      MatQN SigMat( int w, int k ) const; 				///< Return self-energy quantum number matrix for specific momentum and frequency given the current state vector

      dcomplex genchi_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex genchi_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex genchi_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element

      dcomplex genchi_0_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex genchi_0_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      
      dcomplex genchis_plus_30_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;
      dcomplex genchit_minus_30_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex genchi_d( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex genchi_m( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element

      dcomplex vertx( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element for up, down, up, down
      dcomplex vertx_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element
      dcomplex vertx_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element

      dcomplex vertx_upup( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element for up, up, up, up
      dcomplex vertx_ph_upup( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return vertex-tensor element for up, up, up, up

      dcomplex lambda( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const; /// < Return lambda-tensor element

      dcomplex gam_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return gamma-tensor element
      dcomplex gam_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return gamma-tensor element
      dcomplex gam_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return gamma-tensor element

      dcomplex gam_ph_upup( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return gamma-tensor element

      dcomplex phi_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return Phi function for PP-channel
      dcomplex phi_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return Phi function for PH-channel
      dcomplex phi_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return Phi function for XPH-channel

      dcomplex phi_pp_outside( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Asymptotic part of Phi function for PP-channel
      dcomplex phi_ph_outside( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Asymptotic part of Phi function for PH-channel
      dcomplex phi_xph_outside( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Asymptotic part of Phi function for XPH-channel

      dcomplex chi_pp( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PP-channel
      dcomplex chi_ph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex chi_xph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for XPH-channel

      dcomplex chi_pp_t( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PP-channel
      dcomplex chi_pp_s( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PP-channel
      dcomplex chi_pp_d( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PP-channel
      dcomplex chi_pp_m( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PP-channel
      
      dcomplex chi_ph_m( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex chi_ph_d( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for XPH-channel
      dcomplex chi_ph_s( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex chi_ph_t( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for XPH-channel
      
      dcomplex chi_xph_m( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex chi_xph_d( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for XPH-channel
      dcomplex chi_xph_s( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for PH-channel
      dcomplex chi_xph_t( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return chi function for XPH-channel
      
      dcomplex P_pp( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PP-channel
      dcomplex P_ph( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PH-channel
      dcomplex P_xph( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for XPH-channel

      dcomplex P_pp_t( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PP-channel TRIPLET
      dcomplex P_pp_s( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;        ///< Return P function for PP-channel SINGLET
      dcomplex P_ph_m( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PH-channel MAGNETIC
      dcomplex P_ph_d( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return P function for PH-channel DENSITY

      dcomplex R_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return R function for PP-channel
      dcomplex R_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return R function for PH-channel
      dcomplex R_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const;	///< Return R function for XPH-channel

      dcomplex Gval( const idx_1p_t& Gvec, double Lam ) const; 

      MatQN GMat( const idx_1p_mat_t& Gvec, double Lam ) const; 

      inline dcomplex Sig( const idx_1p_t& idx ) const
      {
	 return Sig( idx(0), idx(1), idx(2), idx(3) ); 
      }

      inline MatQN SigMat( const idx_1p_mat_t& idx ) const 
      {
	 return SigMat( idx(0), idx(1) ); 
      }

      inline dcomplex vertx( const idx_2p_t& idx ) const
      {
	 return vertx( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }
 
      inline dcomplex vertx_pp( const idx_phi_t& idx ) const
      {
	 return vertx_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex vertx_ph( const idx_phi_t& idx ) const
      {
	 return vertx_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex vertx_xph( const idx_phi_t& idx ) const
      {
	 return vertx_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }
      
      inline dcomplex genchi_pp( const idx_phi_t& idx ) const
      {
	 return genchi_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex genchi_ph( const idx_phi_t& idx ) const
      {
	 return genchi_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex genchi_xph( const idx_phi_t& idx ) const
      {
	 return genchi_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }
      
      inline dcomplex lambda( const idx_2p_t& idx ) const
      {
	 return lambda( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex phi_pp( const idx_phi_t& idx ) const
      {
	 return phi_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex phi_ph( const idx_phi_t& idx ) const
      {
	 return phi_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex phi_xph( const idx_phi_t& idx ) const
      {
	 return phi_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex chi_pp( const idx_chi_t& idx ) const
      {
	 return chi_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5) ); 
      }

      inline dcomplex chi_ph( const idx_chi_t& idx ) const
      {
	 return chi_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5) ); 
      }

      inline dcomplex chi_xph( const idx_chi_t& idx ) const
      {
	 return chi_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5) ); 
      }

      inline dcomplex R_pp( const idx_phi_t& idx ) const
      {
	 return R_pp( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex R_ph( const idx_phi_t& idx ) const
      {
	 return R_ph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

      inline dcomplex R_xph( const idx_phi_t& idx ) const
      {
	 return R_xph( idx(0), idx(1), idx(2), idx(3), idx(4), idx(5), idx(6), idx(7), idx(8), idx(9) ); 
      }

}; 

//namespace boost {
//   namespace numeric {
//      namespace odeint {
//	 template<>
//	    struct vector_space_norm_inf< state_t >
//	    {
//	       typedef double result_type;
//	       double operator()( const state_t &state_vec ) const
//	       {
//		  using namespace std; 
//		  return norm( state_vec ); 
//	       }
//	    };
//      }
//   }
//}

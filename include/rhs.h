
/*******************************************************************************************//** @file
 *  		
 * 	file: 		rhs.h
 * 	contents:  	Functor specifying the calculation of the rhs of the flow equation
 * 
 ****************************************************************************************************/


#pragma once
#include <state.h>
#include <grid.h>
#include <symmetry_group.h>
#include <symmetries.h>

class rhs_t 		///< Functor to specify the rhs calculation for both the self-energy and the vertex
{
   public:
      void operator() ( const state_t& state_vec, state_t &dfdl, const double Lam ); 	///< Overload call operator to calculate full rhs

      /********************* RHS evaluation ********************/
      static dcomplex eval_rhs_Sig( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam ); 		///< Calculate rhs of self-energy flow for specific index
      
      static dcomplex eval_diag_pp( const idx_phi_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam );		///< Calculate particle-particle diagram for the flow
      static dcomplex eval_diag_ph( const idx_phi_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam );		///< Calculate particle-hole diagram for the flow
      static dcomplex eval_diag_xph( const idx_phi_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam );		///< Calculate particle-hole diagram for the flow
      
      static void phi_pp_inverse( const state_t& state_vec, gf_phi_t& gf_phi_pp, const gf_1p_mat_t& Gvec, double Lam ); 
      static void phi_ph_xph_inverse( const state_t& state_vec, gf_phi_t& gf_phi_ph, gf_phi_t& gf_phi_xph, const gf_1p_mat_t& Gvec, double Lam ); 

   private:
      static symmetry_grp_t<4> symm_grp_sig; 

      static symmetry_grp_t<10> symm_grp_phi_pp; 
      static symmetry_grp_t<10> symm_grp_phi_ph; 
      static symmetry_grp_t<10> symm_grp_phi_xph; 
      
      static symmetry_grp_t<8> symm_grp_P_pp; 
      static symmetry_grp_t<8> symm_grp_P_ph; 
      static symmetry_grp_t<8> symm_grp_P_xph; 
      
      static symmetry_grp_t<6> symm_grp_chi_pp; 
      static symmetry_grp_t<6> symm_grp_chi_ph; 
      static symmetry_grp_t<6> symm_grp_chi_xph; 
      
      static gf<double, 1> weight_vec; 
      static gf<double, 2> weight_vec_2d; 

};



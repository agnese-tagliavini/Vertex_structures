
/************************************************************************************************//**
 *  		
 * 	file: 		rhs.cpp
 * 	contents:   	See rhs.h
 * 
 ****************************************************************************************************/


#include <rhs.h>
#include <frg.h>
#include <params.h>
#include <observables.h>
#include <cmath>
#include <translate.h>
#include <mymath.h>
#include <Eigen/LU>

using namespace std; 
using boost::bind;

// --- Create static objects: symmetry groups, weight vector

symmetry_grp_t<4> rhs_t::symm_grp_sig( gf_1p_t(), { cconj_sig } ); 

symmetry_grp_t<10> rhs_t::symm_grp_phi_pp( gf_phi_t(), { } ); 
symmetry_grp_t<10> rhs_t::symm_grp_phi_ph( gf_phi_t(), { } ); 
symmetry_grp_t<10> rhs_t::symm_grp_phi_xph( gf_phi_t(), { } ); 


symmetry_grp_t<8> rhs_t::symm_grp_P_pp( gf_P_t(),{
//      hmirror_P_pp
//      , cconj_timerev_P_pp
      } ); 

symmetry_grp_t<8> rhs_t::symm_grp_P_ph( gf_P_t(),{ 
//      hmirror_timerev_P_ph
//      , vmirror_P_ph
      } ); 

symmetry_grp_t<8> rhs_t::symm_grp_P_xph( gf_P_t(),{ 
//      cconj_P_xph, 
//      timerev_P_xph
      } ); 

symmetry_grp_t<6> rhs_t::symm_grp_chi_pp( gf_chi_t(), { 
//      cconj_chi_pp
      } ); 
symmetry_grp_t<6> rhs_t::symm_grp_chi_ph( gf_chi_t(), { 
//      hmirror_chi_ph
//      , cconj_chi_ph
      } ); 

symmetry_grp_t<6> rhs_t::symm_grp_chi_xph( gf_chi_t(), { 
//      hmirror_chi_xph
//     , cconj_chi_xph
//      , timerev_chi_xph 
      } ); 

//symmetry_grp_t<8> rhs_t::symm_grp_P_pp( gf_P_t(), { } ); 
//symmetry_grp_t<8> rhs_t::symm_grp_P_xph( gf_P_t(), { } ); 
//symmetry_grp_t<8> rhs_t::symm_grp_P_ph( gf_P_t(), { } ); 
//
//symmetry_grp_t<6> rhs_t::symm_grp_chi_pp( gf_chi_t(), { } ); 
//symmetry_grp_t<6> rhs_t::symm_grp_chi_ph( gf_chi_t(), {  } ); 
//symmetry_grp_t<6> rhs_t::symm_grp_chi_xph( gf_chi_t(), { } ); 


gf<double, 1> rhs_t::weight_vec( generate_weights( POS_INT_RANGE-TAIL_LENGTH, TAIL_LENGTH, FIT_ORDER ) ); 
gf<double, 2> rhs_t::weight_vec_2d( generate_2d_weights( POS_INT_RANGE-TAIL_LENGTH, TAIL_LENGTH, FIT_ORDER ) ); 


void rhs_t::operator() ( const state_t& state_vec, state_t& dfdl, const double Lam )
{

   // Precalculate Green function and single-scale propagator on frequency grid

   gf_1p_mat_t Gvec( POS_1P_RANGE ); 
   Gvec.init( bind( &state_t::GMat, boost::cref(state_vec), _1, Lam ) ); // Initialize big Green function vector 

   /*********************  Open MP parallelized RHS calculation  ********************/

   cout << " ... Rhs calculation: ";

#ifndef READIN 

   cout << " ... Sig " << flush;
   symm_grp_sig.init( dfdl.gf_Sig(), [&state_vec, &Gvec, &Lam]( const idx_1p_t& idx ){ return eval_rhs_Sig( idx, state_vec, Gvec, Lam ); } );

#endif

   cout << " ... chi " << endl;
   symm_grp_chi_pp.init( dfdl.gf_chi_pp(), [&state_vec, &Gvec, Lam]( const idx_chi_t& idx ){ return eval_diag_chi_pp( idx, state_vec, Gvec, Lam ); } ); 
   symm_grp_chi_ph.init( dfdl.gf_chi_ph(), [&state_vec, &Gvec, Lam]( const idx_chi_t& idx ){ return eval_diag_chi_ph( idx, state_vec, Gvec, Lam ); } ); 
   symm_grp_chi_ph.init( dfdl.gf_chi_xph(), [&state_vec, &Gvec, Lam]( const idx_chi_t& idx ){ return eval_diag_chi_xph( idx, state_vec, Gvec, Lam ); } ); 

   cout << " ... P " << endl << endl;

   symm_grp_P_pp.init( dfdl.gf_P_pp(), [&state_vec, &dfdl, &Gvec, Lam]( const idx_P_t& idx ){ return eval_diag_P_pp( idx, state_vec, Gvec, Lam ) - dfdl.gf_chi_pp()( iP_to_ichi( idx ) ); } ); 
   symm_grp_P_ph.init( dfdl.gf_P_ph(), [&state_vec, &dfdl, &Gvec, Lam]( const idx_P_t& idx ){ return eval_diag_P_ph( idx, state_vec, Gvec, Lam ) - dfdl.gf_chi_ph()( iP_to_ichi( idx ) ); } ); 
   symm_grp_P_ph.init( dfdl.gf_P_xph(), [&state_vec, &dfdl, &Gvec, Lam]( const idx_P_t& idx ){ return eval_diag_P_xph( idx, state_vec, Gvec, Lam ) - dfdl.gf_chi_xph()( iP_to_ichi( idx ) ); } ); 
   

}

dcomplex rhs_t::eval_rhs_Sig( const idx_1p_t& se_idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

   int w = se_idx( I1P::w ); 
   int k = se_idx( I1P::k ); 
   int s_in = se_idx( I1P::s_in ); 
   int s_out = se_idx( I1P::s_out ); 

   for( int s1 = 0; s1 < QN_COUNT; ++s1 )
      for( int s1p = 0; s1p < QN_COUNT; ++s1p )
	 for( int s2 = 0; s2 < QN_COUNT; ++s2 )
	    for( int s2p = 0; s2p < QN_COUNT; ++s2p )
	       for( int s3 = 0; s3 < QN_COUNT; ++s3 )
		  for( int s3p = 0; s3p < QN_COUNT; ++s3p )
		     for( int k1 = 0; k1 < PATCH_COUNT; ++k1 )
			for( int w1 = -POS_INT_RANGE; w1 < POS_INT_RANGE; ++w1 )
			   for( int k2 = 0; k2 < PATCH_COUNT; ++k2 )
			      for( int w2 = -POS_INT_RANGE; w2 < POS_INT_RANGE; ++w2 )
			      {
				 val += vert_bare( s_in, s3, s1p, s2p ) * Gvec[w1][k1]( s1p, s1 ) * Gvec[w2][k2]( s2p, s2 ) 
				    * weight_vec_2d[w1][w2]  // Do not use 1d fitting here!!
				    * Gvec[ w1 + w2 - w ][ add_k( k1, dif_k( k2, k ) ) ]( s3p, s3 ) 
				    * state_vec.vertx( w1, w2, w, k1, k2, k, s1, s2, s_out, s3p ); 
			      }

   return - val / ( BETA * PATCH_COUNT ) / ( BETA * PATCH_COUNT );  // Minus sign from loop
}

dcomplex rhs_t::eval_diag_chi_pp( const idx_chi_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

#ifdef FORCED_ZEROS
   if( forced_zero_check( idx ) )							// check wether element should be forced to zero
      return val;
#endif

   // Introduce help variables
   int W = idx( ICHI::W );

   int K = idx( ICHI::K );

   int s1_in = idx( ICHI::s1_in );
   int s2_in = idx( ICHI::s2_in );
   int s1_out = idx( ICHI::s1_out );
   int s2_out = idx( ICHI::s2_out );
      
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	 for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	    for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	          for( int k = 0; k < PATCH_COUNT; ++k ){
		     int K_minus_k = dif_k( K, k ); 
		     for( int w = -POS_INT_RANGE ; w < POS_INT_RANGE; ++w ){
		        for( int s5 = 0; s5 < QN_COUNT; ++s5 )
                           for( int s5p = 0; s5p < QN_COUNT; ++s5p )
                   	      for( int s6 = 0; s6 < QN_COUNT; ++s6 )
                   	         for( int s6p = 0; s6p < QN_COUNT; ++s6p )
                                       for( int kp = 0; kp < PATCH_COUNT; ++kp ){
		     		          int K_minus_kp = dif_k( K, kp ); 
				       for( int wp = -POS_INT_RANGE ; wp < POS_INT_RANGE; ++wp )
		  			{
					     // Particle particle channel
						val += (state_vec.genchi_pp(W, w, wp, K, k, kp, s3, s4, s5p, s6p )
						      +state_vec.genchi_pp(W, w, -wp -1, K, k, kp, s3, s4, s5p, s6p ) + 
						      2*state_vec.genchi_0_pp(W, w, wp, K, k, kp, s3, s4, s5p, s6p))
						     * weight_vec_2d[w][wp];  // Do not use 1d fitting here!!


		  			}
				       }
				    }
				 }

      val *= 0.5*vert_bare(s1_in, s2_in, s1_out, s2_out)*vert_bare(s1_in, s2_in, s1_out, s2_out)/BETA/BETA/PATCH_COUNT/PATCH_COUNT;
      return val;
}

dcomplex rhs_t::eval_diag_P_pp( const idx_P_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

#ifdef FORCED_ZEROS
   if( forced_zero_check( idx ) )							// check wether element should be forced to zero
      return val;
#endif

   // Introduce help variables
   int W = idx( IP::W );
   int w_in = idx( IP::w);

   int K = idx( IP::K );
   int k_in = idx( IP::k );

   int s1_in = idx( IP::s1_in );
   int s2_in = idx( IP::s2_in );
   int s1_out = idx( IP::s1_out );
   int s2_out = idx( IP::s2_out );

      
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	 for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	    for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
		  for( int k = 0; k < PATCH_COUNT; ++k )
		  {

		     int K_minus_k = dif_k( K, k ); 

		     // Particle particle channel
		     for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
		     {
			val += state_vec.genchi_pp(W, w_in, w, K, k_in, k, s1_in, s2_in, s3p, s4p )
			   * weight_vec[w]; //READJUST for varying summation range
		     }

		  }

   val *= vert_bare(0,0,0,0)/ BETA / PATCH_COUNT * 1./Gvec[ w_in + div2_ceil( W ) ][0]( 0,0 ) * 1./Gvec[ div2_floor( W ) - w_in - 1 ][0]( 0, 0 );
  return val; 
}

dcomplex rhs_t::eval_diag_chi_ph( const idx_chi_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

#ifdef FORCED_ZEROS
   if( forced_zero_check( idx ) )							// check wether element should be forced to zero
      return val;
#endif

   // Introduce help variables
   int W = idx( ICHI::W );

   int K = idx( ICHI::K );

   int s1_in = idx( ICHI::s1_in );
   int s2_in = idx( ICHI::s2_in );
   int s1_out = idx( ICHI::s1_out );
   int s2_out = idx( ICHI::s2_out );

      
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	 for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	    for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
		  for( int s5 = 0; s5 < QN_COUNT; ++s5 )
                     for( int s5p = 0; s5p < QN_COUNT; ++s5p )
                   	for( int s6 = 0; s6 < QN_COUNT; ++s6 )
                   	   for( int s6p = 0; s6p < QN_COUNT; ++s6p )
			      for( int k = 0; k < PATCH_COUNT; ++k ){
		                 int K_plus_k = add_k( K, k ); 
                                 for( int kp = 0; kp < PATCH_COUNT; ++kp ){
		                   int K_plus_kp = add_k( K, kp ); 
				    for( int w = -POS_INT_RANGE ; w < POS_INT_RANGE; ++w ){
				       for( int wp = -POS_INT_RANGE ; wp < POS_INT_RANGE; ++wp )
					     // Particle hole channel
					     {
						val += 
						    state_vec.genchi_ph( W, w, wp, K, k, kp, s3, s5, s4p, s6p )
                                 	            * weight_vec_2d[w][wp];  // Do not use 1d fitting here!! 
					     }
					}
				    }
			      }

   val *= vert_bare(0,0,0,0)*vert_bare(0,0,0,0)*1.0/BETA/BETA/PATCH_COUNT/PATCH_COUNT;  // Two interanal loops-> plus sign in front
   return val; 
}

dcomplex rhs_t::eval_diag_P_ph( const idx_P_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

#ifdef FORCED_ZEROS
   if( forced_zero_check( idx ) )							// check wether element should be forced to zero
      return val;
#endif

   // Introduce help variables
   int W = idx( IP::W );
   int w_in = idx( IP::w );

   int K = idx( IP::K );
   int k_in = idx( IP::k );

   int s1_in = idx( IP::s1_in );
   int s2_in = idx( IP::s2_in );
   int s1_out = idx( IP::s1_out );
   int s2_out = idx( IP::s2_out );

      
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	 for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	    for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
		  for( int k = 0; k < PATCH_COUNT; ++k )
		  {

		     int K_plus_k = add_k( K, k ); 

		     // Particle hole channel
		     for( int w = -POS_INT_RANGE; w < POS_INT_RANGE; ++w )
		     {
			val += weight_vec[w] * //READJUST for varying summation range
			   (state_vec.genchi_ph( W, w_in, w, K, k_in, k, s1_in, s3, s1_out, s4p )-
			    state_vec.genchi_xph( W, w_in, w, K, k_in, k, s1_in, s3, s1_out, s4p ));
		     }


		  }

   val *= -vert_bare(0,0,0,0)/ BETA / PATCH_COUNT * 1./Gvec[ w_in + div2_ceil( W ) ][0](0) * 1./Gvec[ -div2_floor( W ) + w_in ][0]( 0, 0 );  // minus sign for the internal loop 
   return val-vert_bare(0,0,0,0); 
}

dcomplex rhs_t::eval_diag_chi_xph( const idx_chi_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

#ifdef FORCED_ZEROS
   if( forced_zero_check( idx ) )							// check wether element should be forced to zero
      return val;
#endif

   // Introduce help variables
   int W = idx( ICHI::W );

   int K = idx( ICHI::K );

   int s1_in = idx( ICHI::s1_in );
   int s2_in = idx( ICHI::s2_in );
   int s1_out = idx( ICHI::s1_out );
   int s2_out = idx( ICHI::s2_out );

      
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	 for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	    for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
		  for( int s5 = 0; s5 < QN_COUNT; ++s5 )
                     for( int s5p = 0; s5p < QN_COUNT; ++s5p )
                   	for( int s6 = 0; s6 < QN_COUNT; ++s6 )
                   	   for( int s6p = 0; s6p < QN_COUNT; ++s6p )
			      for( int k = 0; k < PATCH_COUNT; ++k ){
		                 int K_plus_k = add_k( K, k ); 
                                 for( int kp = 0; kp < PATCH_COUNT; ++kp ){
		                   int K_plus_kp = add_k( K, kp ); 
				    for( int w = -POS_INT_RANGE ; w < POS_INT_RANGE; ++w ){
				       for( int wp = -POS_INT_RANGE ; wp < POS_INT_RANGE; ++wp )
				       {
					val += 
					   state_vec.genchi_xph( W, wp, w, K, kp, k, s3, s5, s4p, s6p )
					   * weight_vec_2d[w][wp];
				       }
				    }
				 }
			      }


   val *= 1.0/BETA/BETA/PATCH_COUNT/PATCH_COUNT * vert_bare(0,0,0,0) * vert_bare(0,0,0,0);
   return val; 
}
dcomplex rhs_t::eval_diag_P_xph( const idx_P_t& idx, const state_t& state_vec, const gf_1p_mat_t& Gvec, double Lam )
{
   dcomplex val( 0.0, 0.0 );

#ifdef FORCED_ZEROS
   if( forced_zero_check( idx ) )							// check wether element should be forced to zero
      return val;
#endif

   // Introduce help variables
   int W = idx( IP::W );
   int w_in = idx( IP::w );

   int K = idx( IP::K );
   int k_in = idx( IP::k );
   
   int s1_in = idx( IP::s1_in );
   int s2_in = idx( IP::s2_in );
   int s1_out = idx( IP::s1_out );
   int s2_out = idx( IP::s2_out );

      
   for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	 for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	    for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
		  for( int k = 0; k < PATCH_COUNT; ++k )
		  {
		     int K_plus_k = add_k( K, k ); 

		     // Particle hole exchange channel, this way of connecting has no loop and no swap
		     for( int w = -POS_INT_RANGE ; w < POS_INT_RANGE ; ++w )
		     {
			val += 
			   state_vec.genchi_xph( W, w_in, w, K, k_in, k, s1_in, s3, s2_out, s4p ) *  
			   weight_vec[w]; 
		     }

		  }    

   val *= vert_bare(0,0,0,0)/ BETA / PATCH_COUNT * 1./Gvec[ w_in + div2_ceil( W ) ][0](0,0 ) * 1./Gvec[-div2_floor( W ) + w_in ][0]( 0, 0 ); 
   return val-vert_bare(0,0,0,0); 
}


/*********************************************************************************************************************************************************
 *
 * 				BETHE-SALPETER EQUATIONS INVERTED WITH CORRECTION TO THE FINITE MATRIX INVERSION:
 * 				1) AGNESE'S CORRECTIONS 2) STEFAN'S CORRECTIONS 3) BIG MATRIX WITH ASYMPTOTICS USES(TODO:CHECK NOTATION-> NON TESTED)
 *
 ***********************************************************************************************************************************************************/ 				

#ifdef METHOD2 // AGNESE'S METHOD FOR THE INVERSION OF THE BETHE-SALPETER EQUATIONS -> VERTEX ASYMPTOTICS USED 
void rhs_t::phi_pp_inverse( const state_t& state_vec, gf_phi_t& gf_phi_pp, const gf_1p_mat_t& Gvec, double Lam )
{
   using gf_mat_t = gf<dcomplex, 8>; 
   using gf_mat_corr_t = gf<dcomplex, 6>; 
#ifndef SINGLETHREADED 
#pragma omp parallel for schedule( dynamic )
#endif
   for( int W = -POS_BFREQ_COUNT_PHI; W < POS_BFREQ_COUNT_PHI + 1; ++W )
   //for( int W = 20; W < 20 + 1; ++W )
   {
      const int MAT_DIM = 2 * POS_INV_RANGE * PATCH_COUNT * QN_COUNT * QN_COUNT; 

      gf_mat_t gam_t( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] ); 
      gf_mat_t gam_s(gam_t);
      gf_mat_t gam_0(gam_t);
      gf_mat_t corr_gamma_s(gam_t);
      
      for( int K = 0; K < PATCH_COUNT; ++K )
      {
	 // Initialize gam_t with genchi_t_minus_0_pp

	 gam_t.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchit_minus_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 
	 gam_s.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchis_plus_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } );
	 //gam_s.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } );
	 gam_0.init([W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } );
   
	 // Build 'Identity + weight_vec * G * G * F' matrixx
	 corr_gamma_s.init( [W,K,&Gvec,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       dcomplex val( 0.0, 0.0 );

	       int w_in = idx(0);
	       int k_in = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s1_in= idx(2);
	       int s2_in = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);

	       for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	       for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       {
	       val += Gvec[div2_floor( W )-w_in-1 ][dif_k(K,k_in)]( s2_in,s4 )*Gvec[ div2_ceil( W )+w_in ][k_in]( s1_in,s3 )*
		  (2.0 * vert_bare(s3,s4,s3p,s4p)+state_vec.P_pp_s(W,w_in,K,k_in,s3,s4,s3p,s4p)+state_vec.chi_pp_s(W,K,s3,s4,s3p,s4p))*
		  0.5*vert_bare(s3p,s4p,s1_out,s2_out)*FUNC_PP(W,POS_INV_RANGE); //Check sign of UINT
	       
	       }

	       
	       val *= -1.0 / BETA / BETA; 

	       return val; 
	       } );
         
	 MapXcd( gam_s.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_s.data(), MAT_DIM, MAT_DIM ).inverse() * dcomplex(-4.0,0.0);
#ifdef CORRECTIONS
	 MapXcd( corr_gamma_s.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_s.data(), MAT_DIM, MAT_DIM ) * MapXcd( corr_gamma_s.data(), MAT_DIM, MAT_DIM );
#endif	 
	 
	 MapXcd( gam_t.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_t.data(), MAT_DIM, MAT_DIM ).inverse() * dcomplex(-4.0,0.0);
	 MapXcd( gam_0.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_0.data(), MAT_DIM, MAT_DIM ).inverse() * dcomplex(2.0,0.0);
         
	 
	 gam_s += gam_0;
	 gam_t -= gam_0;
#ifdef CORRECTIONS
	 gam_s -= corr_gamma_s;
#endif
	 gam_s *= 1.0 * BETA * BETA; 
         gam_t *= 1.0 * BETA * BETA; 

	 for( int w_in = -POS_FFREQ_COUNT_PHI; w_in < POS_FFREQ_COUNT_PHI; ++w_in )
	    for( int w_out = -POS_FFREQ_COUNT_PHI; w_out < POS_FFREQ_COUNT_PHI; ++w_out )
	       for( int k_in = 0; k_in < PATCH_COUNT; ++k_in )
		  for( int k_out = 0; k_out < PATCH_COUNT; ++k_out )
		     for( int s1 = 0; s1 < QN_COUNT; ++s1 )
			for( int s2 = 0; s2 < QN_COUNT; ++s2 )
			   for( int s1p = 0; s1p < QN_COUNT; ++s1p )
			      for( int s2p = 0; s2p < QN_COUNT; ++s2p ){
				//int Wb = W - 20;
				 gf_phi_pp[W][w_in][w_out][K][k_in][k_out][s1][s2][s1p][s2p] = 
				    state_vec.vertx_pp( W, w_in, w_out, K, k_in, k_out, s1, s2, s1p, s2p ) - 
				    0.5*(gam_s[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]+gam_t[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]); //GAMMA UPDO PP = (GAMMA_S + GAMMA_T)/2 
			      }
      }
   }
}

void rhs_t::phi_ph_xph_inverse( const state_t& state_vec, gf_phi_t& gf_phi_ph, gf_phi_t& gf_phi_xph, const gf_1p_mat_t& Gvec, double Lam )
{
   using gf_mat_t = gf<dcomplex, 8>; 
#ifndef SINGLETHREADED 
#pragma omp parallel for schedule( dynamic )
#endif
   for( int W = -POS_BFREQ_COUNT_PHI; W < POS_BFREQ_COUNT_PHI + 1; ++W )
   //for( int W = 20; W < 20 + 1; ++W )
   {

      const int MAT_DIM = 2 * POS_INV_RANGE * PATCH_COUNT * QN_COUNT * QN_COUNT; 

      gf_mat_t gam_m( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] ); 
      gf_mat_t gam_0(gam_m);
      gf_mat_t gam_d(gam_m);
      gf_mat_t corr_gamma_m(gam_m);
      gf_mat_t corr_gamma_d(gam_m);
      
      for( int K = 0; K < PATCH_COUNT; ++K )
      {
	 // Initialize gam_m gam_d and gam_0
	 gam_m.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_m( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 
	 gam_0.init([W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_0_ph( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7)) ; } );
	 gam_d.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_d( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } );
	 
	 // Corrections to the density channel
	 corr_gamma_d.init( [W,K,&Gvec,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       dcomplex val( 0.0, 0.0 );

	       int w_in = idx(0);
	       int k_in = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s1_in= idx(2);
	       int s2_in = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);


	       for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	       for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       {
	       val += Gvec[div2_ceil( W )+w_in ][add_k(K,k_in)]( s2_in,s4 )*Gvec[-div2_floor( W )+w_in ][k_in]( s1_in,s3 )* 
		  (vert_bare(s3,s4,s3p,s4p)+state_vec.P_ph_d(W,w_in,K,k_in,s3,s4,s3p,s4p)+state_vec.chi_ph_d(W,K,s3,s4,s3p,s4p))*
		  vert_bare(s3p,s4p,s1_out,s2_out)*FUNC_PH(W, POS_INV_RANGE); //check sign of UINT
	       }
	       
	       val *= -1.0 / BETA / BETA; 

	       return val; 
	       } ); 

	 //Corrections to the magnetic channel
	 corr_gamma_m.init( [W,K,&Gvec,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       dcomplex val( 0.0, 0.0 );

	       int w_in = idx(0);
	       int k_in = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s1_in= idx(2);
	       int s2_in = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);


	       for( int s3 = 0; s3 < QN_COUNT; ++s3 )
	       for( int s4 = 0; s4 < QN_COUNT; ++s4 )
	       for( int s3p = 0; s3p < QN_COUNT; ++s3p )
	       for( int s4p = 0; s4p < QN_COUNT; ++s4p )
	       {
	       val += Gvec[div2_ceil( W )+w_in ][add_k(K,k_in)]( s2_in,s4 )*Gvec[-div2_floor( W )+w_in ][k_in]( s1_in,s3 )* 
		  (-vert_bare(s3,s4,s3p,s4p)+state_vec.P_ph_m(W,w_in,K,k_in,s3,s4,s3p,s4p)+state_vec.chi_ph_m(W,K,s3,s4,s3p,s4p))*
		  -1.0 * vert_bare(s3p,s4p,s1_out,s2_out)*FUNC_PH(W, POS_INV_RANGE); //check sign of UINT
	       }
	       
	       val *=-1.0 / BETA / BETA; 

	       return val; 
	       } ); 
#ifdef CORRECTIONS   
	 MapXcd( corr_gamma_m.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_m.data(), MAT_DIM, MAT_DIM ).inverse() * (MapXcd( corr_gamma_m.data(), MAT_DIM, MAT_DIM ));
	 MapXcd( corr_gamma_d.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_d.data(), MAT_DIM, MAT_DIM ).inverse() * (MapXcd( corr_gamma_d.data(), MAT_DIM, MAT_DIM ));
#endif
	 MapXcd( gam_0.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_0.data(), MAT_DIM, MAT_DIM ).inverse(); 
         MapXcd( gam_d.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_d.data(), MAT_DIM, MAT_DIM ).inverse();
	 MapXcd( gam_m.data(), MAT_DIM, MAT_DIM ) = MapXcd( gam_m.data(), MAT_DIM, MAT_DIM ).inverse();
	 
	 gam_d += gam_0;
	 gam_m += gam_0;
#ifdef CORRECTIONS	 
	 gam_d -= corr_gamma_d;
	 gam_m -= corr_gamma_m;
#endif
	 gam_d *= -BETA*BETA;
	 gam_m *= -BETA*BETA;
	 
         
	 for( int w_in = -POS_FFREQ_COUNT_PHI; w_in < POS_FFREQ_COUNT_PHI; ++w_in )
	    for( int w_out = -POS_FFREQ_COUNT_PHI; w_out < POS_FFREQ_COUNT_PHI; ++w_out )
	       for( int k_in = 0; k_in < PATCH_COUNT; ++k_in )
		  for( int k_out = 0; k_out < PATCH_COUNT; ++k_out )
		     for( int s1 = 0; s1 < QN_COUNT; ++s1 )
			for( int s2 = 0; s2 < QN_COUNT; ++s2 )
			   for( int s1p = 0; s1p < QN_COUNT; ++s1p )
			      for( int s2p = 0; s2p < QN_COUNT; ++s2p ){
				//int Wb = W - 20;
				 gf_phi_ph[W][w_in][w_out][K][k_in][k_out][s1][s2][s1p][s2p] = 
				    state_vec.vertx_ph( W, w_in, w_out, K, k_in, k_out, s1, s2, s1p, s2p ) 
				   - 0.5*(gam_d[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]-gam_m[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]); //GAMMA UPDO PH = (GAMMA_D - GAMMA_M)/2 
      
				 gf_phi_xph[W][w_in][w_out][K][k_in][k_out][s1][s2][s1p][s2p] = 
				    state_vec.vertx_xph( W, w_in, w_out, K, k_in, k_out, s1, s2, s1p, s2p ) 
				   + gam_m[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]; //GAMMA UPDO XPH = - GAMMA_M 
			      }
      }
   }

}

#endif

#ifdef METHOD1 // STEFAN'S METHOD FOR THE INVERSION OF THE BETHE-SALPETER EQUATIONS-> GAMMA'S ASYMPTOTICS USED


void rhs_t::phi_pp_inverse( const state_t& state_vec, gf_phi_t& gf_phi_pp, const gf_1p_mat_t& Gvec, double Lam )
{
   gf<double, 2> weight_vec_2d_asy = generate_2d_weights_asy( POS_ASY_RANGE + POS_INV_RANGE-TAIL_LENGTH_ASY, TAIL_LENGTH_ASY, FIT_ORDER_ASY ); 
   using gf_mat_t = gf<dcomplex, 8>; 
  
#ifndef SINGLETHREADED 
#pragma omp parallel for schedule( dynamic )
#endif
   for( int W = -POS_BFREQ_COUNT_PHI; W < POS_BFREQ_COUNT_PHI + 1; ++W )

   //for( int W = 20; W < 20 + 1; ++W )
   {
      const int MAT_DIM = 2 * POS_INV_RANGE * PATCH_COUNT * QN_COUNT * QN_COUNT;
      const int MAT_11_DIM = 2 * POS_ASY_RANGE * PATCH_COUNT * QN_COUNT * QN_COUNT;
      
      gf_mat_t gam_s_00_pp( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      gf_mat_t gam_t_00_pp( gam_s_00_pp ); 
      gf_mat_t gam_0_00_pp( gam_s_00_pp);
      gf_mat_t gam_11_pp( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      gf_mat_t gam_0_11_pp( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_01_pp( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_10_pp( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_corr_pp( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      
      for( int K = 0; K < PATCH_COUNT; ++K )
      {
	 // Initialize gam_t, gam_s in all the matrix blocks
	 gam_t_00_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchit_minus_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 

	 gam_s_00_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchis_plus_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 
	 
	 gam_0_00_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_0_pp( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 
	 
	 //gam_01_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return 2.0*vert_bare( idx(2), idx(3), idx(6), idx(7) ); } );
      	 //
	 //gam_10_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return 2.0*vert_bare( idx(2), idx(3), idx(6), idx(7) ); } );
      	 //
         //gam_corr_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return dcomplex(0.0,0.0); } );
#ifdef CORRECTIONS
	 gam_0_11_pp.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       int w = idx(0);
	       int k = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s3p = idx(2);
	       int s4p = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);

	       int K_minus_k = dif_k( K, k ); 

	       if ((w < 0) && (w_out < 0)){

	          w = w-POS_INV_RANGE;
		  w_out = w_out-POS_INV_RANGE;
	       	  return state_vec.genchi_0_pp( W, w, w_out, K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) );
	       }

	       else if((w >= 0) && (w_out >= 0)){
		   
		   w = w+POS_INV_RANGE;
		   w_out = w_out+POS_INV_RANGE;

	       	  return state_vec.genchi_0_pp( W, w, w_out, K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) );
	       
	       }
               else
		  return dcomplex(0.0,0.0);
   	} );

	 gam_11_pp.init( [W,K,&Gvec,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 

	       int w = idx(0);
	       int k = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s3p = idx(2);
	       int s4p = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);

	       int K_minus_k = dif_k( K, k ); 

	       if ((w < 0) && (w_out < 0)){

	          w = w-POS_INV_RANGE;
		  w_out = w_out-POS_INV_RANGE;
	          
	          return -1./BETA/BETA/2.0 * (state_vec.chi_ph_s(-w-w_out-(W+10000)%2-1,dif_k(K_minus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_s(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) + 2.0 * vert_bare(s3p,s4p,s1_out, s2_out)) * state_vec.genchi_0_pp( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out ) +
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       
	       }

	       else if((w >= 0) && (w_out >= 0)){
		   
		   w = w+POS_INV_RANGE;
		   w_out = w_out+POS_INV_RANGE;

	          return -1./BETA/BETA/2.0 * (state_vec.chi_ph_s(-w-w_out-(W+10000)%2-1,dif_k(K_minus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_s(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out)+ 2.0 * vert_bare(s3p,s4p,s1_out, s2_out)) * state_vec.genchi_0_pp( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out ) +
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       
	       }


	       else if ((w < 0) && (w_out >= 0)){
	          
		   w = w-POS_INV_RANGE;
		   w_out = w_out+POS_INV_RANGE;

	          return -1./BETA/BETA/2.0 * (state_vec.chi_ph_s(-w-w_out-(W+10000)%2-1,dif_k(K_minus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_s(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out)+ 2.0 * vert_bare(s3p,s4p,s1_out, s2_out)) * state_vec.genchi_0_pp( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out ) +
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }

	       else if ((w >= 0) && (w_out < 0)){
	          
		   w = w+POS_INV_RANGE;
		   w_out = w_out-POS_INV_RANGE;

	          return -1./BETA/BETA/2.0 * (state_vec.chi_ph_s(-w-w_out-(W+10000)%2-1,dif_k(K_minus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_s(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out)+ 2.0 * vert_bare(s3p,s4p,s1_out, s2_out)) * state_vec.genchi_0_pp( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out ) +
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }
	       else
		  return dcomplex(0.0,0.0);

	       } ); 

	 
	 MapXcd(  gam_11_pp.data(), MAT_11_DIM, MAT_11_DIM ) = MapXcd(  gam_11_pp.data(), MAT_11_DIM, MAT_11_DIM ).inverse();
         //gam_11_pp *= gam_0_11_pp;
	 //MapXcd(  gam_corr_pp.data(), MAT_11_DIM, MAT_DIM ) = MapXcd(  gam_11_pp.data(), MAT_11_DIM, MAT_11_DIM ) * MapXcd(  gam_10_pp.data(), MAT_11_DIM, MAT_DIM );

	 dcomplex corr_inv(0.0,0.0); 
	 for (int w1=-POS_ASY_RANGE; w1<POS_ASY_RANGE; ++w1){
	    for (int w2=-POS_ASY_RANGE; w2<POS_ASY_RANGE; ++w2){
	    corr_inv +=  0.5/BETA/BETA/BETA/BETA * (2*vert_bare(0,0,0,0)) * (2*vert_bare(0,0,0,0)) * 
	                 gam_11_pp[w1][0][0][0][w2][0][0][0] *  gam_0_11_pp[w2][0][0][0][w2][0][0][0]
	        	 * weight_vec_2d_asy[w1][w2];
	    }
	 }
	 cout << "CORR DONE" << endl;
#endif

         MapXcd(  gam_t_00_pp.data(), MAT_DIM, MAT_DIM ) = MapXcd(  gam_t_00_pp.data(), MAT_DIM, MAT_DIM ).inverse();
	 MapXcd(  gam_s_00_pp.data(), MAT_DIM, MAT_DIM ) = MapXcd(  gam_s_00_pp.data(), MAT_DIM, MAT_DIM ).inverse();
         MapXcd(  gam_0_00_pp.data(), MAT_DIM, MAT_DIM ) = MapXcd(  gam_0_00_pp.data(), MAT_DIM, MAT_DIM ).inverse();
	 
	 MapXcd(  gam_t_00_pp.data(), MAT_DIM, MAT_DIM ) *= dcomplex(-4.0,0.0);
	 MapXcd(  gam_s_00_pp.data(), MAT_DIM, MAT_DIM ) *= dcomplex(-4.0,0.0);
	 MapXcd(  gam_0_00_pp.data(), MAT_DIM, MAT_DIM ) *= dcomplex(2.0,0.0);

         
	 gam_t_00_pp -= gam_0_00_pp;
	 gam_s_00_pp += gam_0_00_pp;

#ifdef CORRECTIONS
	 //MapXcd( gam_s_00_pp.data(), MAT_DIM, MAT_DIM ) -= 1./BETA/BETA/BETA/BETA/2.0*MapXcd(  gam_01_pp.data(), MAT_DIM, MAT_11_DIM ) * MapXcd(  gam_corr_pp.data(), MAT_11_DIM, MAT_DIM ); 
	 gam_s_00_pp -= corr_inv;
#endif         

	 gam_t_00_pp *= BETA*BETA; 
	 gam_s_00_pp *= BETA*BETA;

	 for( int w_in = -POS_FFREQ_COUNT_PHI; w_in < POS_FFREQ_COUNT_PHI; ++w_in )
	    for( int w_out = -POS_FFREQ_COUNT_PHI; w_out < POS_FFREQ_COUNT_PHI; ++w_out )
	       for( int k_in = 0; k_in < PATCH_COUNT; ++k_in )
		  for( int k_out = 0; k_out < PATCH_COUNT; ++k_out )
		     for( int s1 = 0; s1 < QN_COUNT; ++s1 )
			for( int s2 = 0; s2 < QN_COUNT; ++s2 )
			   for( int s1p = 0; s1p < QN_COUNT; ++s1p )
			      for( int s2p = 0; s2p < QN_COUNT; ++s2p ){
				 //int Wb = W - 20;
				 gf_phi_pp[W][w_in][w_out][K][k_in][k_out][s1][s2][s1p][s2p] = 
				    state_vec.vertx_pp( W, w_in, w_out, K, k_in, k_out, s1, s2, s1p, s2p ) - 
				    0.5*(gam_s_00_pp[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]+gam_t_00_pp[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]); 
					  }
	}
   }
}

void rhs_t::phi_ph_xph_inverse( const state_t& state_vec, gf_phi_t& gf_phi_ph, gf_phi_t& gf_phi_xph, const gf_1p_mat_t& Gvec, double Lam )
{
   gf<double, 2> weight_vec_2d_asy = generate_2d_weights_asy( POS_ASY_RANGE+POS_INV_RANGE-TAIL_LENGTH_ASY, TAIL_LENGTH_ASY, FIT_ORDER_ASY ) ; 
   using gf_mat_t = gf<dcomplex, 8>; 
#ifndef SINGLETHREADED 
#pragma omp parallel for schedule( dynamic )
#endif
   for( int W = -POS_BFREQ_COUNT_PHI; W < POS_BFREQ_COUNT_PHI + 1; ++W )
     
   //for( int W = 20; W < 20 + 1; ++W )
   {
      const int MAT_DIM = 2 * POS_INV_RANGE * PATCH_COUNT * QN_COUNT * QN_COUNT;
      const int MAT_11_DIM = 2 * POS_ASY_RANGE * PATCH_COUNT * QN_COUNT * QN_COUNT;

      gf_mat_t gam_d_00_ph( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      gf_mat_t gam_m_00_ph( gam_d_00_ph ); 
      gf_mat_t gam_0_00_ph( gam_d_00_ph);
      
      gf_mat_t gam_d_11_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      gf_mat_t gam_m_11_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      gf_mat_t gam_0_11_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_d_01_ph( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_d_10_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_m_01_ph( boost::extents[ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_m_10_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_corr_m_10_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] );
      //gf_mat_t gam_corr_d_10_ph( boost::extents[ffreq(POS_ASY_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT][ffreq(POS_INV_RANGE)][PATCH_COUNT][QN_COUNT][QN_COUNT] ); 
      
      for( int K = 0; K < PATCH_COUNT; ++K )
      {
	 // Initialize gam_pp with full vertex
	 gam_d_00_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_d( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 

	 gam_m_00_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_m( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } ); 
	 
	 gam_0_00_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return state_vec.genchi_0_ph( W, idx(0), idx(4), K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) ); } );
	 
	 //gam_d_01_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return vert_bare( idx(2), idx(3), idx(6), idx(7) ); } );
      	 //gam_d_10_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return vert_bare( idx(2), idx(3), idx(6), idx(7) ); } );
      	 //gam_m_01_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return -vert_bare( idx(2), idx(3), idx(6), idx(7) ); } );
      	 //gam_m_10_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx ){ return -vert_bare( idx(2), idx(3), idx(6), idx(7) ); } );
       	 
#ifdef CORRECTIONS	 
	 gam_0_11_ph.init( [W,K,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       int w = idx(0);
	       int k = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s3p = idx(2);
	       int s4p = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);

	       int K_plus_k = add_k( K, k ); 

	       if ((w < 0) && (w_out<0)){

	          w = w-POS_INV_RANGE;
		  w_out = w_out-POS_INV_RANGE;
	       
		  return state_vec.genchi_0_ph( W, w, w_out, K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) );
	          
	       }

	       if((w >= 0) && (w_out>=0)){
		   w = w+POS_INV_RANGE;
		   w_out = w+POS_INV_RANGE;

		  return state_vec.genchi_0_ph( W, w, w_out, K, idx(1), idx(5), idx(2), idx(3), idx(6), idx(7) );
	       
	       }
	       else
		  return dcomplex(0.0,0.0);
   } );
	       
	       

	 // Build 'Identity + weight_vec * G * G * F' matrixx
	 gam_d_11_ph.init( [W,K,&Gvec,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       dcomplex val( 0.0, 0.0 );

	       int w = idx(0);
	       int k = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s3p = idx(2);
	       int s4p = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);

	       int K_plus_k = add_k( K, k ); 

	       if ((w < 0) && (w_out<0)){

	          w = w-POS_INV_RANGE;
		  w_out = w_out-POS_INV_RANGE;
	          
	          return 1./BETA/BETA*(state_vec.chi_pp_d(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_d(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) + vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }

	       if((w >= 0) && (w_out>=0)){
		   w = w+POS_INV_RANGE;
		   w_out = w+POS_INV_RANGE;

	          return 1./BETA/BETA*(state_vec.chi_pp_d(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_d(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) + vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       
	       }


	       if ((w < 0) && (w_out >= 0)){
	          
		   w = w-POS_INV_RANGE;
		   w_out = w+POS_INV_RANGE;

	          return 1./BETA/BETA*(state_vec.chi_pp_d(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_d(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) + vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }

	       if ((w >= 0) && (w_out < 0)){
	          
		   w = w+POS_INV_RANGE;
		   w_out = w-POS_INV_RANGE;

	          return 1./BETA/BETA*(state_vec.chi_pp_d(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_d(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) + vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }
	       else
		  return dcomplex(0.0,0.0);
	       } ); 


	 gam_m_11_ph.init( [W,K,&Gvec,&state_vec]( const gf_mat_t::idx_t& idx )
	       { 
	       dcomplex val( 0.0, 0.0 );

	       int w = idx(0);
	       int k = idx(1);

	       int w_out = idx(4);
	       int k_out = idx(5);

	       int s3p = idx(2);
	       int s4p = idx(3);
	       int s1_out = idx(6);
	       int s2_out = idx(7);

	       int K_plus_k = add_k( K, k ); 

	       if ((w < 0) && (w_out<0)){

	          w = w-POS_INV_RANGE;
		  w_out = w_out-POS_INV_RANGE;
	          
	          return 1./BETA/BETA*(state_vec.chi_pp_m(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_m(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) - vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }

	       if((w >= 0) && (w_out>=0)){
		   w = w+POS_INV_RANGE;
		   w_out = w_out+POS_INV_RANGE;

	          return 1./BETA/BETA*(state_vec.chi_pp_m(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_m(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) - vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       
	       }


	       if ((w < 0) && (w_out >= 0)){
	          
		   w = w-POS_INV_RANGE;
		   w_out = w_out+POS_INV_RANGE;

	          return 1./BETA/BETA*(state_vec.chi_pp_m(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_m(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) - vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }

	       if ((w >= 0) && (w_out < 0)){
	          
		   w = w+POS_INV_RANGE;
		   w_out = w_out-POS_INV_RANGE;

	          return 1./BETA/BETA*(state_vec.chi_pp_m(w+w_out+(W+10000)%2+1,add_k(K_plus_k,k_out),s3p,s4p,s1_out,s2_out) + 
		         state_vec.chi_xph_m(w_out-w,dif_k(k_out,k),s3p,s4p,s1_out,s2_out) - vert_bare(s3p,s4p,s1_out, s2_out))* state_vec.genchi_0_ph( W, w, w_out, K, k, k_out, s3p, s4p, s1_out, s2_out )+
		         1.0 * ( idx(0) == idx(4) ) * ( idx(1) == idx(5) ) * ( idx(2) == idx(6) ) * ( idx(3) == idx(7) );
	       }
	       else
		  return dcomplex(0.0,0.0);
	       } ); 

         
	 MapXcd(  gam_d_11_ph.data(), MAT_11_DIM, MAT_11_DIM ) = MapXcd(  gam_d_11_ph.data(), MAT_11_DIM, MAT_11_DIM ).inverse();
	
	 MapXcd(  gam_m_11_ph.data(), MAT_11_DIM, MAT_11_DIM ) = MapXcd(  gam_m_11_ph.data(), MAT_11_DIM, MAT_11_DIM ).inverse();

	 //gam_d_11_ph *= gam_0_11_ph;
	 //gam_m_11_ph *= gam_0_11_ph;
	 //  
	 //MapXcd(  gam_corr_d_10_ph.data(), MAT_11_DIM, MAT_DIM ) = MapXcd(  gam_d_11_ph.data(), MAT_11_DIM, MAT_11_DIM ) * MapXcd(  gam_d_10_ph.data(), MAT_11_DIM, MAT_DIM );
	 //MapXcd(  gam_corr_m_10_ph.data(), MAT_11_DIM, MAT_DIM ) = MapXcd(  gam_m_11_ph.data(), MAT_11_DIM, MAT_11_DIM ) * MapXcd(  gam_m_10_ph.data(), MAT_11_DIM, MAT_DIM ); 

	 dcomplex corr_inv_d(0.0,0.0);
	 dcomplex corr_inv_m(0.0,0.0);
	 
	 for (int w1=-POS_ASY_RANGE; w1<POS_ASY_RANGE; ++w1){
	    for (int w2=-POS_ASY_RANGE; w2<POS_ASY_RANGE; ++w2){
	    
	       corr_inv_d +=  1.0/BETA/BETA/BETA/BETA * (vert_bare(0,0,0,0)) * (vert_bare(0,0,0,0)) * 
	                 gam_d_11_ph[w1][0][0][0][w2][0][0][0] *  gam_0_11_ph[w2][0][0][0][w2][0][0][0]
	        	 * weight_vec_2d_asy[w1][w2];
	    
	    corr_inv_m +=  1.0/BETA/BETA/BETA/BETA * (-vert_bare(0,0,0,0)) * (-vert_bare(0,0,0,0)) * 
	                 gam_m_11_ph[w1][0][0][0][w2][0][0][0] *  gam_0_11_ph[w2][0][0][0][w2][0][0][0]
	        	 * weight_vec_2d_asy[w1][w2];
	    }
	 }

#endif         
	 MapXcd(  gam_d_00_ph.data(), MAT_DIM, MAT_DIM ) = MapXcd(  gam_d_00_ph.data(), MAT_DIM, MAT_DIM ).inverse();
	 MapXcd(  gam_m_00_ph.data(), MAT_DIM, MAT_DIM ) = MapXcd(  gam_m_00_ph.data(), MAT_DIM, MAT_DIM ).inverse();
         MapXcd(  gam_0_00_ph.data(), MAT_DIM, MAT_DIM ) = MapXcd(  gam_0_00_ph.data(), MAT_DIM, MAT_DIM ).inverse();
	 

	 gam_d_00_ph += gam_0_00_ph;
	 gam_m_00_ph += gam_0_00_ph;

#ifdef CORRECTIONS	 
	 //MapXcd(  gam_d_00_ph.data(), MAT_DIM, MAT_DIM ) -= 1./BETA/BETA/BETA/BETA* MapXcd(  gam_d_01_ph.data(), MAT_DIM, MAT_11_DIM )*MapXcd(  gam_corr_d_10_ph.data(), MAT_11_DIM, MAT_DIM ); 
	 //MapXcd(  gam_m_00_ph.data(), MAT_DIM, MAT_DIM ) -= 1./BETA/BETA/BETA/BETA* MapXcd(  gam_m_01_ph.data(), MAT_DIM, MAT_11_DIM )*MapXcd(  gam_corr_m_10_ph.data(), MAT_11_DIM, MAT_DIM );
	 
	 gam_d_00_ph -= 2*corr_inv_d;
         gam_m_00_ph -= 2*corr_inv_m;
#endif
	 gam_d_00_ph *= -BETA*BETA; 
	 gam_m_00_ph *= -BETA*BETA;

	 //copy this part to gf_phi_pp
	 for( int w_in = -POS_FFREQ_COUNT_PHI; w_in < POS_FFREQ_COUNT_PHI; ++w_in )
	    for( int w_out = -POS_FFREQ_COUNT_PHI; w_out < POS_FFREQ_COUNT_PHI; ++w_out )
	       for( int k_in = 0; k_in < PATCH_COUNT; ++k_in )
		  for( int k_out = 0; k_out < PATCH_COUNT; ++k_out )
		     for( int s1 = 0; s1 < QN_COUNT; ++s1 )
			for( int s2 = 0; s2 < QN_COUNT; ++s2 )
			   for( int s1p = 0; s1p < QN_COUNT; ++s1p )
			      for( int s2p = 0; s2p < QN_COUNT; ++s2p ){
				 //int Wb = W - 20;
				 gf_phi_ph[W][w_in][w_out][K][k_in][k_out][s1][s2][s1p][s2p] = 
				    state_vec.vertx_ph( W, w_in, w_out, K, k_in, k_out, s1, s2, s1p, s2p ) - 
				    0.5*(gam_d_00_ph[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]-gam_m_00_ph[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p]);
					  
				 gf_phi_xph[W][w_in][w_out][K][k_in][k_out][s1][s2][s1p][s2p] = 
				    state_vec.vertx_xph( W, w_in, w_out, K, k_in, k_out, s1, s2, s1p, s2p ) + gam_m_00_ph[w_in][k_in][s1][s2][w_out][k_out][s1p][s2p];
					  } 
      }

   }
}


#endif


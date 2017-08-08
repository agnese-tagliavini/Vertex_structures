
/******************************************************************************************//** @file
 *  		
 * 	file: 		state.cpp
 * 	contents:  	See state.h
 * 
 ****************************************************************************************************/


#include <state.h>
#include <mymath.h>
#include <frg.h>

using namespace std; 
/********************* Interfacing gf containers  ********************/

/*************************************************************************************************************************************************
 *
 * 						SELF-ENERGY: NOT USED IF THE SE IS GIVEN AS INPUT
 *
 ************************************************************************************************************************************************/ 						

dcomplex state_t::Sig( int w, int k, int s_in, int s_out ) const
{
   if ( w < -POS_FFREQ_COUNT_SIG )
      return -1.0 * POS_FFREQ_COUNT_SIG / w * gf_Sig()[-POS_FFREQ_COUNT_SIG][k][s_in][s_out]; 

   if ( w > POS_FFREQ_COUNT_SIG - 1 )
      return 1.0 * ( POS_FFREQ_COUNT_SIG - 1 ) / w * gf_Sig()[POS_FFREQ_COUNT_SIG-1][k][s_in][s_out]; 

   //if ( w < -POS_FFREQ_COUNT_SIG || w > POS_FFREQ_COUNT_SIG - 1 ) 
   //   return 0.0; 

   return gf_Sig()[w][k][s_in][s_out]; 
}

MatQN state_t::SigMat( int w, int k ) const
{
   if ( w < -POS_FFREQ_COUNT_SIG )
      return -1.0 * POS_FFREQ_COUNT_SIG / w * Eigen::Map<const MatQN>( &(gf_Sig()[-POS_FFREQ_COUNT_SIG][k][0][0]) ); 
	//return MatQN::Zero();

   if ( w > POS_FFREQ_COUNT_SIG - 1 )
      return 1.0 * ( POS_FFREQ_COUNT_SIG - 1 ) / w * Eigen::Map<const MatQN>( &(gf_Sig()[POS_FFREQ_COUNT_SIG-1][k][0][0]) ); 
	//return MatQN::Zero();

   //if ( w < -POS_FFREQ_COUNT_SIG || w > POS_FFREQ_COUNT_SIG - 1 ) 
   //   return MatQN::Zero(); 

   MatQN SigMatrix = Eigen::Map<const MatQN>( &(gf_Sig()[w][k][0][0]) );  

   return SigMatrix; 
}


/***************************************************************************************************************************************************
 *
 * 						FULL VERTEX IN PURELY FERMIONIC NOTATION
 *
 *************************************************************************************************************************************************/ 						

dcomplex state_t::vertx( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

#ifdef READIN
   
   return vertx_ph( w1_out - w1_in, w1_in + div2_floor( w1_out - w1_in ),w2_in - div2_ceil( w1_out - w1_in ),dif_k(k1_out, k1_in), k1_in, add_k( k1_in, dif_k( k2_in, k1_out ) ), s1_in, s2_in, s1_out, s2_out); // Here the purely fermionic notation is recovered from the ph channel (the other channels would be good as well!)-> It's nor really used if the SE is not calculated using the Schwringer-Dyson equation!!!
 
#endif

   int W_pp  = w1_in + w2_in + 1;
   int W_ph  = w1_out - w1_in;
   int W_xph = w2_in - w1_out;

   int K_pp = add_k( k1_in, k2_in ); 
   int K_ph = dif_k( k1_out, k1_in ); 
   int K_xph = dif_k( k2_in, k1_out ); 

   int k2_out = add_k( k1_in, dif_k( k2_in, k1_out ) ); // calculate k2_out 

   return phi_pp_outside( W_pp, w1_in - div2_ceil( W_pp ), -w1_out-1+ div2_floor( W_pp ), K_pp, k1_in, k1_out, s1_in, s2_in, s1_out, s2_out ) + 
      phi_ph_outside( W_ph, w1_in + div2_floor( W_ph ), w2_in - div2_ceil( W_ph ), K_ph, k1_in, k2_out, s1_in, s2_in, s1_out, s2_out ) +
      phi_xph_outside( W_xph, w1_in + div2_floor( W_xph ), w2_in - div2_ceil( W_xph ), K_xph, k1_in, k1_out, s1_in, s2_in, s1_out, s2_out ) +
      vert_bare( s1_in, s2_in, s1_out, s2_out );
}


/*****************************************************************************************************************************************************************************
 *
 * 						FULL VERTEX IN THE DIFFERENT CHANNEL-DEPENDENT NOTATIONS
 *
 ***************************************************************************************************************************************************************************/

//PP
dcomplex state_t::vertx_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

#ifdef READIN
    //Return exact solution if available
   if (  W  >= -POS_BOS_VERT_COUNT_EXACT_SMALL && W  <= POS_BOS_VERT_COUNT_EXACT_SMALL &&
	 w_in  >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_in  <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 &&
	 w_out >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_out <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 )
      return vert_exact_pp[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; //INPUT DATA FROM POMEROL-> WARNING: DIFFERENT NOTATION!!
#endif

#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

   int W_xph  = w_out - w_in;
   int W_ph  = - ( W + 100000 ) % 2 - w_out - w_in - 1;

   int K_xph = dif_k( k_out, k_in ); 
   int K_ph = dif_k( K, add_k( k_in, k_out ) ); 

   return phi_pp_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) + 
      phi_xph_outside( W_xph, w_in + div2_ceil( W ) + div2_floor( W_xph ), div2_floor( W ) - w_out + div2_floor( W_xph ) - 1, K_xph, k_in, dif_k( K, k_out ), s1_in, s2_in, s1_out, s2_out ) +
      phi_ph_outside( W_ph, w_in + div2_ceil( W ) + div2_floor( W_ph ), w_out + div2_ceil( W ) + div2_floor( W_ph ), K_ph, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) +
      vert_bare( s1_in, s2_in, s1_out, s2_out );
}

//PH
dcomplex state_t::vertx_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

#ifdef READIN
    //Return exact solution if available
   if (  W  >= -POS_BOS_VERT_COUNT_EXACT_SMALL && W  <= POS_BOS_VERT_COUNT_EXACT_SMALL &&
	 w_in  >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_in  <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 &&
	 w_out >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_out <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 )
      return vert_exact_ph[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
#endif

#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif
   int W_pp  = w_in + w_out + ( W + 100000 ) % 2 + 1;
   int W_xph  = w_out - w_in;

   int K_pp = add_k( K, add_k( k_out, k_in ) ); 
   int K_xph = dif_k( k_out, k_in ); 

   return phi_pp_outside( W_pp, w_in - div2_floor( W ) - div2_ceil( W_pp ), w_out - div2_floor( W ) - div2_ceil( W_pp ), K_pp, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) + 
      phi_ph_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) +
      phi_xph_outside( W_xph, w_in - div2_floor( W ) + div2_floor( W_xph ), w_out + div2_ceil( W ) - div2_ceil( W_xph ), K_xph, k_in, add_k( K, k_in ), s1_in, s2_in, s1_out, s2_out ) +
      vert_bare( s1_in, s2_in, s1_out, s2_out );
}

//XPH
dcomplex state_t::vertx_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

#ifdef READIN
    //Return exact solution if available
   if (  W  >= -POS_BOS_VERT_COUNT_EXACT_SMALL && W  <= POS_BOS_VERT_COUNT_EXACT_SMALL &&
	 w_in  >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_in  <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 &&
	 w_out >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_out <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 )
      return vert_exact_xph[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
#endif

#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

   int W_pp  = w_in + w_out + ( W + 100000 ) % 2 + 1;
   int W_ph  = w_out - w_in;

   int K_pp = add_k( K, add_k( k_out, k_in ) ); 
   int K_ph = dif_k( k_out, k_in ); 

   return phi_pp_outside( W_pp, w_in - div2_floor( W ) - div2_ceil( W_pp ), -w_out -1 + div2_floor( W ) + div2_floor( W_pp ), K_pp, k_in, add_k( K,k_in ), s1_in, s2_in, s1_out, s2_out ) + 
      phi_ph_outside( W_ph, w_in - div2_floor( W ) + div2_floor( W_ph ), w_out + div2_ceil( W ) - div2_ceil( W_ph ), K_ph, k_in, add_k( K, k_in ), s1_in, s2_in, s1_out, s2_out ) +
      phi_xph_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) +
      vert_bare( s1_in, s2_in, s1_out, s2_out );
}

dcomplex state_t::vertx_upup( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return vertx( w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out ) - vertx( w2_in, w1_in, w1_out, k2_in, k1_in, k1_out, s2_in, s1_in, s1_out, s2_out ); 
}

dcomplex state_t::vertx_ph_upup( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return vertx_ph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) - vertx_xph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ); 
}

/*************************************************************************************************************************************************************************
 *
 * 						 GENERALIZED SUSCEPTIBILITY -> FOR THE INVERSION OF THE BETHE-SALPETER EQUATIONS
 *
 ************************************************************************************************************************************************************************/ 						 
//PP
dcomplex state_t::genchi_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif
#ifdef READIN
    //Return exact solution if available
   if (  W  >= -POS_BOS_VERT_COUNT_EXACT_SMALL && W  <= POS_BOS_VERT_COUNT_EXACT_SMALL &&
         w_in  >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_in  <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 &&
         w_out >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_out <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 ){
      return genchi_exact_pp[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; //INPUT DATA FROM POMEROL-> WARNING: DIFFERENT NOTATION!!
     }
   
#endif

   return G(w_val(div2_floor(W)-w_in-1),mom_grid[dif_k(K,k_in)].first,mom_grid[dif_k(K,k_in)].second,Lam,SigMat(div2_floor(W)-w_in-1,dif_k(K,k_in)))(s2_in,s2_in)*G(w_val(w_in+div2_ceil(W)),mom_grid[k_in].first,mom_grid[k_in].second,Lam,SigMat(w_in+div2_ceil(W),k_in))(s1_in,s1_in)*vertx_pp(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out)* G(w_val(div2_floor(W)-w_out-1),mom_grid[dif_k(K,k_out)].first,mom_grid[dif_k(K,k_out)].second,Lam,SigMat(div2_floor(W)-w_out-1,dif_k(K,k_out)))(s1_out,s1_out)*G(w_val(w_out+div2_ceil(W)),mom_grid[k_out].first,mom_grid[k_out].second,Lam,SigMat(w_out+div2_ceil(W),k_out))(s2_out,s2_out);


}

//PH
dcomplex state_t::genchi_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif
#ifdef READIN
    //Return exact solution if available
   if (  W  >= -POS_BOS_VERT_COUNT_EXACT_SMALL && W  <= POS_BOS_VERT_COUNT_EXACT_SMALL  &&
	 w_in  >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_in  <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 &&
	 w_out >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_out <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 )
      return genchi_exact_ph[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
#endif

   return G(w_val(-div2_floor(W)+w_in),mom_grid[k_in].first,mom_grid[k_in].second,Lam,SigMat(-div2_floor(W)+w_in,k_in))(s1_in,s1_in)*G(w_val(w_out+div2_ceil(W)),mom_grid[add_k(K,k_out)].first,mom_grid[add_k(K,k_out)].second,Lam,SigMat(w_out+div2_ceil(W),add_k(K,k_out)))(s2_in,s2_in)*vertx_ph(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out)* G(w_val(div2_ceil(W)+w_in),mom_grid[add_k(K,k_in)].first,mom_grid[add_k(K,k_in)].second,Lam,SigMat(div2_ceil(W)+w_in,add_k(K,k_in)))(s1_out,s1_out)*G(w_val(w_out-div2_floor(W)),mom_grid[k_out].first,mom_grid[k_out].second,Lam,SigMat(w_out-div2_floor(W),k_out))(s2_out,s2_out);
}

//XPH
dcomplex state_t::genchi_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif
#ifdef READIN
    //Return exact solution if available
   if (  W  >= -POS_BOS_VERT_COUNT_EXACT_SMALL && W  <= POS_BOS_VERT_COUNT_EXACT_SMALL &&
	 w_in  >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_in  <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 &&
	 w_out >= -POS_FERM_VERT_COUNT_EXACT_SMALL && w_out <= POS_FERM_VERT_COUNT_EXACT_SMALL - 1 )
      return genchi_exact_xph[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
#endif
   if((w_in == w_out) && (k_in == k_out) && (s1_in == s1_out) && (s2_in == s2_out)){
   
      return BETA*G(w_val(-div2_floor(W)+w_in),mom_grid[k_in].first,mom_grid[k_in].second,Lam,SigMat(-div2_floor(W)+w_in,k_in))(s1_in,s1_out)*G(w_val(w_in+div2_ceil(W)),mom_grid[add_k(K,k_in)].first,mom_grid[add_k(K,k_in)].second,Lam,SigMat(w_in+div2_ceil(W),add_k(K,k_in)))(s2_in,s2_out) + G(w_val(-div2_floor(W)+w_in),mom_grid[k_in].first,mom_grid[k_in].second,Lam,SigMat(-div2_floor(W)+w_in,k_in))(s1_in,s1_in)*G(w_val(w_out+div2_ceil(W)),mom_grid[add_k(K,k_out)].first,mom_grid[add_k(K,k_out)].second,Lam,SigMat(w_out+div2_ceil(W),add_k(K,k_out)))(s2_in,s2_in)*vertx_xph(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out)* G(w_val(div2_ceil(W)+w_in),mom_grid[add_k(K,k_in)].first,mom_grid[add_k(K,k_in)].second,Lam,SigMat(div2_ceil(W)+w_in,add_k(K,k_in)))(s2_out,s2_out)*G(w_val(w_out-div2_floor(W)),mom_grid[k_out].first,mom_grid[k_out].second,Lam,SigMat(w_out-div2_floor(W),k_out))(s1_out,s1_out);
   }
      return G(w_val(-div2_floor(W)+w_in),mom_grid[k_in].first,mom_grid[k_in].second,Lam,SigMat(-div2_floor(W)+w_in,k_in))(s1_in,s1_in)*G(w_val(w_out+div2_ceil(W)),mom_grid[add_k(K,k_out)].first,mom_grid[add_k(K,k_out)].second,Lam,SigMat(w_out+div2_ceil(W),add_k(K,k_out)))(s2_in,s2_in)*vertx_xph(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out)* G(w_val(div2_ceil(W)+w_in),mom_grid[add_k(K,k_in)].first,mom_grid[add_k(K,k_in)].second,Lam,SigMat(div2_ceil(W)+w_in,add_k(K,k_in)))(s2_out,s2_out)*G(w_val(w_out-div2_floor(W)),mom_grid[k_out].first,mom_grid[k_out].second,Lam,SigMat(w_out-div2_floor(W),k_out))(s1_out,s1_out);
}

//BARE BUBBLE PP
dcomplex state_t::genchi_0_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if((w_in == w_out) && (k_in == k_out) && (s1_in == s1_out) && (s2_in == s2_out)){
      return BETA*G(w_val(w_in + div2_ceil(W)),mom_grid[k_in].first, mom_grid[k_in].second, Lam, SigMat(w_in + div2_ceil(W),k_in))(s1_in,s1_out)*G(w_val(div2_floor(W) - w_out - 1),mom_grid[dif_k(K,k_out)].first, mom_grid[dif_k(K,k_out)].second, Lam, SigMat(div2_floor(W) - w_out - 1, dif_k(K,k_out)))(s2_in,s2_out);
   }
  return dcomplex(0.0,0.0); 
}

//BARE BUBBLE PH
dcomplex state_t::genchi_0_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if((w_in == w_out) && (k_in == k_out) && (s1_in == s1_out) && (s2_in == s2_out)){
      return BETA*G(w_val(w_in + div2_ceil(W)),mom_grid[add_k(K,k_in)].first, mom_grid[add_k(K,k_in)].second, Lam, SigMat(w_in + div2_ceil(W),add_k(K,k_in)))(s1_in,s1_out)*G(w_val(w_out - div2_floor(W)),mom_grid[k_out].first, mom_grid[k_out].second, Lam, SigMat(w_out - div2_floor(W), k_out))(s2_in,s2_out);
   }
  return dcomplex(0.0,0.0); 
}

/*******************************************************************************************************************************************************************************************
 *
 * 				FUNCTIONS USEFUL FOR THE INVERSION OF THE BETHE SALPETER EQUATIONS
 *
 *******************************************************************************************************************************************************************************************/ 				
dcomplex state_t::genchis_plus_0_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   //cout << " W, w_in, w_out" << W << w_in << w_out << "genchi0_pp" << genchi_0_pp(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out) << endl;
   return genchi_pp(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out)+genchi_pp(W, w_in, -w_out-1-( W + 100000 ) % 2, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out) + 2.*genchi_0_pp(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out) ; 
}


dcomplex state_t::genchit_minus_0_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return genchi_pp(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out)-genchi_pp(W, w_in, -w_out-1-( W + 100000 ) % 2, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out) - 2.*genchi_0_pp(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out) ; 
}

dcomplex state_t::genchi_d( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return 2.0 * genchi_ph(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out)-genchi_xph(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out); 
}

dcomplex state_t::genchi_m( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return -genchi_xph(W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out);
}


/***********************************************************************************************************************************************************
 *
 * 					 TWO-PARTICLE IRREDUCIBLE VERTEX -> GAMMMAS (FOR THE SELF-CONSISTENCY!!) 
 *
 **********************************************************************************************************************************************************/ 					 

dcomplex state_t::gamma_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

   return vertx_pp(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out) - phi_pp(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out);

   //return vertx_pp( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) -
   //phi_pp( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ); 

   //int W_ph  = w_out - w_in;
   //int W_xph  = - ( W + 100000 ) % 2 - w_out - w_in - 1;

   //int K_ph = dif_k( k_out, k_in ); 
   //int K_xph = dif_k( K, add_k( k_in, k_out ) ); 

   //return phi_ph( W_ph, w_in + div2_ceil( W ) + div2_floor( W_ph ), div2_floor( W ) - w_out + div2_floor( W_ph ) - 1, K_ph, k_in, dif_k( K, k_out ), s1_in, s2_in, s1_out, s2_out ) +
      //phi_xph( W_xph, w_in + div2_ceil( W ) + div2_floor( W_xph ), w_out + div2_ceil( W ) + div2_floor( W_xph ), K_xph, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) +
      //vert_bare( s1_in, s2_in, s1_out, s2_out );
}

dcomplex state_t::gamma_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

   return vertx_ph(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out) - phi_ph(W,w_in,w_out,K,k_in,k_out,s1_in,s2_in,s1_out,s2_out);

   //return vertx_ph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) -
   //phi_ph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ); 

   //int W_pp  = w_in + w_out + ( W + 100000 ) % 2 + 1;
   //int W_xph  = w_out - w_in;

   //int K_pp = add_k( K, add_k( k_out, k_in ) ); 
   //int K_xph = dif_k( k_out, k_in ); 

   //return phi_pp( W_pp, w_in - div2_floor( W ) - div2_ceil( W_pp ), w_in + div2_ceil( W ) - div2_ceil( W_pp ), K_pp, k_in, add_k( K, k_in ), s1_in, s2_in, s1_out, s2_out ) + 
      //phi_xph( W_xph, w_in - div2_floor( W ) + div2_floor( W_xph ), w_out + div2_ceil( W ) - div2_ceil( W_xph ), K_xph, k_in, add_k( K, k_in ), s1_in, s2_in, s1_out, s2_out ) +
      //vert_bare( s1_in, s2_in, s1_out, s2_out );
}

dcomplex state_t::gamma_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

   return vertx_xph(  W,  w_in,  w_out,  K,  k_in,  k_out,  s1_in,  s2_in,  s1_out,  s2_out) - phi_xph(  W,  w_in,  w_out,  K,  k_in,  k_out,  s1_in,  s2_in,  s1_out, s2_out);

   //return vertx_xph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) -
   //phi_xph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ); 

   //int W_pp  = w_in + w_out + ( W + 100000 ) % 2 + 1;
   //int W_ph  = w_out - w_in;

   //int K_pp = add_k( K, add_k( k_out, k_in ) ); 
   //int K_ph = dif_k( k_out, k_in ); 

   //return phi_pp( W_pp, w_in - div2_floor( W ) - div2_ceil( W_pp ), w_out - div2_floor( W ) - div2_ceil( W_pp ), K_pp, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) + 
      //phi_ph( W_ph, w_in - div2_floor( W ) + div2_floor( W_ph ), w_out + div2_ceil( W ) - div2_ceil( W_ph ), K_ph, k_in, add_k( K, k_in ), s1_in, s2_in, s1_out, s2_out ) +
      //vert_bare( s1_in, s2_in, s1_out, s2_out );
}

dcomplex state_t::gamma_ph_upup( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return gamma_ph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) - gamma_xph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ); //SU2 SYMMETRY USED 
}

// FULLY IRREDUCIBLE VERTEX -> LAMBDA
dcomplex state_t::lambda( int w1_in, int w2_in, int w1_out, int k1_in, int k2_in, int k1_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
#ifdef FORCED_ZEROS
   if( forced_zero_check( s1_in, s2_in, s1_out, s2_out ) )							// check wether element should be forced to zero
      return dcomplex( 0.0, 0.0 );
#endif

   int W_pp  = w1_in + w2_in + 1;
   int W_ph  = w1_out - w1_in;
   int W_xph = w2_in - w1_out;

   int K_pp = add_k( k1_in, k2_in ); 
   int K_ph = dif_k( k1_out, k1_in ); 
   int K_xph = dif_k( k2_in, k1_out ); 

   int k2_out = add_k( k1_in, dif_k( k2_in, k1_out ) ); // calculate k2_out 

   return vertx( w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out ) - 
      phi_pp( W_pp, w1_in - div2_ceil( W_pp ), -w1_out-1 + div2_floor( W_pp ), K_pp, k1_in, k1_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_ph( W_ph, w1_in + div2_floor( W_ph ), w2_in - div2_ceil( W_ph ), K_ph, k1_in, k2_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_xph( W_xph, w1_in + div2_floor( W_xph ), w2_in - div2_ceil( W_xph ), K_xph, k1_in, k1_out, s1_in, s2_in, s1_out, s2_out ); 
}

/**************************************************************************************************************************************
 *
 * 				 TWO-PARTICLE REDUCUIBLE VERTICES-> PHI
 *
 ************************************************************************************************************************************/ 				 

//PP
dcomplex state_t::phi_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return phi_pp_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out );

   return gf_phi_pp()[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
}

//ASYMPTOTIC PP
dcomplex state_t::phi_pp_outside( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return chi_pp( W, K, s1_in, s2_in, s1_out, s2_out ) + 
      P_pp( W, w_in, K, k_in, s1_in, s2_in, s1_out, s2_out ) + 
      P_pp( W, -w_out-1-(W+10000)%2, K, k_out, s1_out, s2_out, s1_in, s2_in );  // time reversal symmetry used
}
//PH
dcomplex state_t::phi_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return phi_ph_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out );

   return gf_phi_ph()[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
}

//ASYMPTOTIC PH
dcomplex state_t::phi_ph_outside( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return chi_ph( W, K, s1_in, s2_in, s1_out, s2_out ) +
      P_ph( W, w_in, K, k_in, s1_in, s2_in, s1_out, s2_out ) + 
      P_ph( -W, w_out, neg_k(K), add_k(k_out,K), s2_in, s1_in, s2_out, s1_out );  // flip diagram horizontally
}

//XPH
dcomplex state_t::phi_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return phi_xph_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out );

   return gf_phi_xph()[W][w_in][w_out][K][k_in][k_out][s1_in][s2_in][s1_out][s2_out]; 
}

//ASYMPTOTIC XPH
dcomplex state_t::phi_xph_outside( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return chi_xph( W, K, s1_in, s2_in, s1_out, s2_out ) + 
      P_xph( W, w_in, K, k_in, s1_in, s2_in, s1_out, s2_out ) + 
      P_xph( -W, w_out, neg_k(K), add_k(k_out,K), s2_in, s1_in, s2_out, s1_out );  // flip diagram horizontally
}


/*****************************************************************************************************************
 * 				CHI
 ****************************************************************************************************************/

//PP
dcomplex state_t::chi_pp( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   

   if ( W < -POS_BFREQ_COUNT_CHI || W > POS_BFREQ_COUNT_CHI ) 
      return 1.0 * POS_BFREQ_COUNT_CHI * POS_BFREQ_COUNT_CHI / W / W * gf_chi_pp()[sgn(W)*POS_BFREQ_COUNT_CHI][K][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;
   return gf_chi_pp()[W][K][s1_in][s2_in][s1_out][s2_out]; 
}

//PH
dcomplex state_t::chi_ph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{

   if ( W < -POS_BFREQ_COUNT_CHI || W > POS_BFREQ_COUNT_CHI ) 
      return 1.0 * POS_BFREQ_COUNT_CHI * POS_BFREQ_COUNT_CHI * POS_BFREQ_COUNT_CHI * POS_BFREQ_COUNT_CHI / W / W / W / W * gf_chi_ph()[sgn(W)*POS_BFREQ_COUNT_CHI][K][s1_in][s2_in][s1_out][s2_out]; 
      //  return 0.0;
   return gf_chi_ph()[W][K][s1_in][s2_in][s1_out][s2_out]; 
}

//XPH
dcomplex state_t::chi_xph( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   

   if ( W < -POS_BFREQ_COUNT_CHI || W > POS_BFREQ_COUNT_CHI ) 
      return 1.0 * POS_BFREQ_COUNT_CHI * POS_BFREQ_COUNT_CHI / W / W * gf_chi_xph()[sgn(W)*POS_BFREQ_COUNT_CHI][K][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;
   
   return gf_chi_xph()[W][K][s1_in][s2_in][s1_out][s2_out]; 
}

/************************************************************************************************************
 *
 * MAGNETIC, DENSITY, SINGLET, TRIPLET SPIN COMBINATION USEFUL FOR BETHE SALPETER EQUATIONS
 * 
 * **********************************************************************************************************/
// MAGNETIC UPUP-UPDO

dcomplex state_t::chi_pp_m( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return -chi_pp(W,K,s1_in,s2_in,s1_out,s2_out); 
}

dcomplex state_t::chi_ph_m( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return -chi_xph(W,K,s1_in,s2_in,s1_out,s2_out); 
}


dcomplex state_t::chi_xph_m( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return -chi_ph(W,K,s1_in,s2_in,s1_out,s2_out); 
}

// DENSITY UPUP+UPDO

dcomplex state_t::chi_pp_d( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return chi_pp(W,K,s1_in,s2_in,s1_out,s2_out); 
}

dcomplex state_t::chi_ph_d( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return 2.*chi_ph(W,K,s1_in,s2_in,s1_out,s2_out)-chi_xph(W,K,s1_in,s2_in,s1_out,s2_out); 
}

dcomplex state_t::chi_xph_d( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return 2.*chi_xph(W,K,s1_in,s2_in,s1_out,s2_out)-chi_ph(W,K,s1_in,s2_in,s1_out,s2_out); 
}

// SINGLET UPDO-XUPDO

dcomplex state_t::chi_pp_s( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return 2.*chi_pp(W,K,s1_in,s2_in,s1_out,s2_out); 
}

dcomplex state_t::chi_ph_s( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return chi_ph(W,K,s1_in,s2_in,s1_out,s2_out)+chi_xph(W,K,s1_in,s2_in,s1_out,s2_out); 
}

dcomplex state_t::chi_xph_s( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return chi_xph(W,K,s1_in,s2_in,s1_out,s2_out)+chi_ph(W,K,s1_in,s2_in,s1_out,s2_out); 
}


// TRIPLET UPDO+XUPDO

dcomplex state_t::chi_pp_t( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return dcomplex(0.0,0.0); 
}

dcomplex state_t::chi_ph_t( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return chi_ph(W,K,s1_in,s2_in,s1_out,s2_out)-chi_xph(W,K,s1_in,s2_in,s1_out,s2_out); 
}

dcomplex state_t::chi_xph_t( int W, int K, int s1_in, int s2_in, int s1_out, int s2_out ) const
{   
   return chi_xph(W,K,s1_in,s2_in,s1_out,s2_out)-chi_ph(W,K,s1_in,s2_in,s1_out,s2_out); 
}

/************************************************************************************************************
 *
 * 			KERNEL 2 (P-function)
 *
 ***********************************************************************************************************/ 			

//PP
dcomplex state_t::P_pp( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_P || W > POS_BFREQ_COUNT_P )
      return 0.0; 
      
   if( w < -POS_FFREQ_COUNT_P ) // Check validity out of PH symmetry! ( imaginary part )
      return ( - 0.25 * W * W + POS_FFREQ_COUNT_P * POS_FFREQ_COUNT_P ) / ( - 0.25 * W * W + w * w ) * gf_P_pp()[W][-POS_FFREQ_COUNT_P][K][k][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;

   if( w > POS_FFREQ_COUNT_P - 1 )
      return ( - 0.25 * W * W + ( POS_FFREQ_COUNT_P - 1 ) * ( POS_FFREQ_COUNT_P - 1 ) ) / ( - 0.25 * W * W + w * w ) * gf_P_pp()[W][POS_FFREQ_COUNT_P-1][K][k][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;

   return gf_P_pp()[W][w][K][k][s1_in][s2_in][s1_out][s2_out]; 
}

//PH
dcomplex state_t::P_ph( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_P || W > POS_BFREQ_COUNT_P )
      return 0.0; 
      
   if( w < -POS_FFREQ_COUNT_P ) // Check validity out of PH symmetry! ( imaginary part )
     return ( - 0.25 * W * W + POS_FFREQ_COUNT_P * POS_FFREQ_COUNT_P ) / ( - 0.25 * W * W + w * w ) * gf_P_ph()[W][-POS_FFREQ_COUNT_P][K][k][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;
   if( w > POS_FFREQ_COUNT_P - 1 )
      return ( - 0.25 * W * W + ( POS_FFREQ_COUNT_P - 1 ) * ( POS_FFREQ_COUNT_P - 1 ) ) / ( - 0.25 * W * W + w * w ) * gf_P_ph()[W][POS_FFREQ_COUNT_P-1][K][k][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;

   return gf_P_ph()[W][w][K][k][s1_in][s2_in][s1_out][s2_out]; 
}

//XPH
dcomplex state_t::P_xph( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_P || W > POS_BFREQ_COUNT_P )
      return 0.0; 
      
   if( w < -POS_FFREQ_COUNT_P ) // Check validity out of PH symmetry! ( imaginary part )
      return ( - 0.25 * W * W + POS_FFREQ_COUNT_P * POS_FFREQ_COUNT_P ) / ( - 0.25 * W * W + w * w ) * gf_P_xph()[W][-POS_FFREQ_COUNT_P][K][k][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;
   if( w > POS_FFREQ_COUNT_P - 1 )
      return ( - 0.25 * W * W + ( POS_FFREQ_COUNT_P - 1 ) * ( POS_FFREQ_COUNT_P - 1 ) ) / ( - 0.25 * W * W + w * w ) * gf_P_xph()[W][POS_FFREQ_COUNT_P-1][K][k][s1_in][s2_in][s1_out][s2_out]; 
	//return 0.0;

   return gf_P_xph()[W][w][K][k][s1_in][s2_in][s1_out][s2_out]; 
}

/************************************************************************************************************
 *
 * MAGNETIC, DENSITY, SINGLET, TRIPLET SPIN COMBINATION USEFUL FOR BETHE SALPETER EQUATIONS
 * 
 * **********************************************************************************************************/
// MAGNETIC UPUP-UPDO

dcomplex state_t::P_ph_m( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return -P_xph(W,w,K,k,s1_in,s2_in,s1_out,s2_out); 
}

//DENSITY UPUP+UPDO

dcomplex state_t::P_ph_d( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return 2.*P_ph(W,w,K,k,s1_in,s2_in,s1_out,s2_out)-P_xph(W,w,K,k,s1_in,s2_in,s1_out,s2_out); 
}

//SINGLET UPDO+XUPDO

dcomplex state_t::P_pp_s( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return 2.*P_pp(W,w,K,k,s1_in,s2_in,s1_out,s2_out); 
}

//TRIPLET UPDO-XUPDO

dcomplex state_t::P_pp_t( int W, int w, int K, int k, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   return 0.0; 
}

/******************************************************************************************************************
 *
 * 			REST FUNCTIONS
 *
 ****************************************************************************************************************/


//PP
dcomplex state_t::R_pp( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return 0.0;

   return phi_pp( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_pp_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out );  
}

//PH
dcomplex state_t::R_ph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return 0.0;

   return phi_ph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_ph_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ); 
}

//XPH
dcomplex state_t::R_xph( int W, int w_in, int w_out, int K, int k_in, int k_out, int s1_in, int s2_in, int s1_out, int s2_out ) const
{
   if ( W < -POS_BFREQ_COUNT_PHI || W > POS_BFREQ_COUNT_PHI ||
	 w_in < -POS_FFREQ_COUNT_PHI || w_in > POS_FFREQ_COUNT_PHI - 1 ||
	 w_out < -POS_FFREQ_COUNT_PHI || w_out > POS_FFREQ_COUNT_PHI - 1 )
      return 0.0;

   return phi_xph( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out ) -
      phi_xph_outside( W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out );  
}

//GREEN's FUNCTIONS -> IN THIS CASE NOT NECESSARY SINCE THE LAMBDA IS SET TO 1 AND THE SE IS FIXED (INPUT VALUE)

dcomplex state_t::Gval( const idx_1p_t& idx, double Lam ) const
{
   return G( w_val( idx(0) ), mom_grid[ idx(1) ].first, mom_grid[ idx(1) ].second, Lam, SigMat( idx(0), idx(1) ) )( idx(2), idx(3) ); 
}

MatQN state_t::GMat( const idx_1p_mat_t& idx, double Lam ) const
{
   return G( w_val( idx(0) ), mom_grid[ idx(1) ].first, mom_grid[ idx(1) ].second, Lam, SigMat( idx(0), idx(1) ) ); 
}

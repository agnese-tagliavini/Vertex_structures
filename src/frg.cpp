
/************************************************************************************************//**
 *  		
 * 	file: 		frg.cpp
 * 	contents:  	for further documentation see frg.h
 * 
 ****************************************************************************************************/


#include <frg.h>
#include <params.h>
#include <Eigen/Eigenvalues>
#include <cmath>
#include <mymath.h>
#include <H5Tools.h>
#include <iostream>

using namespace Eigen;
using namespace std;

/********************* Exact data ********************/

#ifdef READIN

const int POS_FERM_VERT_COUNT_EXACT = 40;
const int POS_BOS_VERT_COUNT_EXACT = 60;

gf_phi_t vert_exact_pp(POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT); 
gf_phi_t vert_exact_ph(POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT);
gf_phi_t vert_exact_xph(POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT);


gf_phi_t genchi_exact_pp(POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT); 
gf_phi_t genchi_exact_ph(POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT);
gf_phi_t genchi_exact_xph(POS_FERM_VERT_COUNT_EXACT, POS_BOS_VERT_COUNT_EXACT);


const int POS_SIG_COUNT_EXACT = 160; 
gf_1p_t Sig_exact(POS_SIG_COUNT_EXACT); 

void read_exact()
{
   using namespace H5;
   H5File input_file( "dat/U1.0_beta20.0_FFREQ_40_BFREQ_60.h5", H5F_ACC_RDONLY );
   cout << "Got file" << endl;

   Group vert_group =  input_file.openGroup("/VERT");
   
   cout << "Got vert" << endl;
   Group vert_group_pp = vert_group.openGroup("PP");
   cout << "Got PP" << endl; 
   read( vert_exact_pp, vert_group_pp,"_F_UPDO");
   cout << "Got vert pp data" << endl; 
   
   Group vert_group_ph = vert_group.openGroup("PH");
   read( vert_exact_ph, vert_group_ph,"_F_UPDO");
   cout << "Got vert ph data" << endl;

   Group vert_group_xph = vert_group.openGroup("XPH");
   read( vert_exact_xph, vert_group_xph,"_F_UPDO");
   cout << "Got vert xph data" << endl;

   Group genchi_group =  input_file.openGroup("/GENCHI"); 
   
   Group genchi_group_pp =  genchi_group.openGroup("PP");
   read( genchi_exact_pp, genchi_group_pp, "_GENCHI_UPDO" );
   cout << "Genchi PP" << endl;
   
   Group genchi_group_ph =  genchi_group.openGroup("PH");
   read( genchi_exact_ph, genchi_group_ph, "_GENCHI_UPDO" );
   cout << "Genchi PH" << endl;
   
   Group genchi_group_xph =  genchi_group.openGroup("XPH");
   read( genchi_exact_xph, genchi_group_xph, "_GENCHI_UPDO" );
   cout << "Genchi XPH" << endl;
   
   Group Sig_group =  input_file.openGroup("/Sig"); 
   read( Sig_exact, Sig_group , "");

   cout << "SE" << endl;
}

#endif

/********************* Hybridization function ********************/

#ifdef ED_BATH

MatQN Gam( double w )
{
   dcomplex val = 0.0; 
   for( int ed_idx = 0; ed_idx < energies.size(); ++ed_idx )
      val += hybridizations[ed_idx]*hybridizations[ed_idx] /( I*w - energies[ed_idx]) ; 

   MatQN Gam; 
   Gam <<
      val;
   return Gam; 
}

#elif QMC_BATH

MatQN Gam( double w )
{
   dcomplex val = 0.0;
   val = I / 2.0 * ( w - sgn( w ) * sqrt( 4.0 + w*w ) ); 

   MatQN Gam; 
   Gam <<
      val; 
   return Gam; 
}

#else // Wide band 

MatQN Gam( double w )
{
   double sq = sqrt( w*w + DEL*DEL );
   double dt = DEL / sq * DD;

   MatQN Gam; 
   Gam << 
      -I*w/sq; 
   return Gam; 
}

#endif


/*********************  INTERACTION FlOW  ********************/

#ifdef INT_FLOW

MatQN G( double w, double kx, double ky, double Lam, const MatQN& selfEn )
{
   return Lam * ( G0inv( w, kx, ky ) - Lam * selfEn ).inverse();
}

dcomplex asympt_GG_pp( int W_int, double Lam )
{
   double PIR = POS_INT_RANGE + abs(W_int/2); 
   double W = W_int; 

   if( W_int == 0 )
      return (BETA*Lam*Lam)/(2.0*pow(PI,2.0)*PIR); 

   return 0.0;
}

#endif

/************************ Common to all implemented cutoff schemes ********************************/

#ifdef NO_MOMENTA

// SIAM / SQDJJ, G0 scale independent for interaction cutoff, Matrix in Nambu basis
MatQN G0inv( double w, double kx, double ky )
{
   MatQN ginv; 
   ginv <<		
      I*w - EPS - B;

   return ginv;
   //return ginv - Gam( w ); 
}

#else

// Hubbard model, G0 scale independent for interaction cutoff, Matrix in Nambu basis
MatQN G0inv( double w, double kx, double ky )
{
   double cos_kx  = cos( kx );
   double cos_ky  = cos( ky );

   MatQN G0inv; 
   G0inv <<		
      I*w + 2.0 * ( cos_kx + cos_ky ) + 4.0 * cos_kx * cos_ky * T_PRIME - B - MU;

   return G0inv; 
}

#endif

dcomplex asympt_GG_ph( int W, double Lam )
{
   return -asympt_GG_pp( W, Lam ); 
}

dcomplex FUNC_PH( int W_int, int POS_INV )
{
   int W = W_int; 

   if( W_int == 0 )
      return -BETA*BETA/(PI*PI)*2.*(1.0/(2.0+4.0*POS_INV)); 

   //return 0.0;
   return BETA*BETA/(PI*PI)*(1.0/(4*W))*log((1.0+2.0*POS_INV-2.0*div2_floor(W))*(1.0+2.0*POS_INV-2.0*div2_ceil(W))/((1.0+2.0*POS_INV+2.0*div2_ceil(W))*(1.0+2.0*POS_INV+2.0*div2_floor(W))));
      
}

dcomplex FUNC_PP( int W_int,int POS_INV )
{
     return -FUNC_PH( W_int, POS_INV ); 
}

// ---- Initial values 

dcomplex vert_bare( int s1_in, int s2_in, int s1_out, int s2_out ) // bare interaction vertex ( Nambu basis )
{
   return -UINT;
}

dcomplex vert_bare( const idx_2p_t& idx ) // bare interaction vertex
{
   return vert_bare( idx( I2P::s1_in ), idx( I2P::s2_in ), idx( I2P::s1_out ), idx( I2P::s2_out ) ); // return bare interaction 
}

dcomplex Sig_init( const idx_1p_t& idx ) // initial values for Sigma
{
#ifdef READIN
   if ((idx(0) >= -POS_SIG_COUNT_EXACT) && (idx(0) <= POS_SIG_COUNT_EXACT - 1) )
      return Sig_exact[idx(0)][idx(1)][idx(2)][idx(3)]; 
#endif

   return 0.0;
}

dcomplex phi_init( const idx_phi_t& idx ) // initial values for phi function
{
   return 0.0;
}

dcomplex P_init( const idx_P_t& idx ) // initial values for P function
{
   return 0.0;
}

dcomplex chi_init( const idx_chi_t& idx ) // inital values for chi function
{
   return 0.0;
}
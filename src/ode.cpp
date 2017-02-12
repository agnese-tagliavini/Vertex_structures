
/************************************************************************************************//**
 *  		
 * 	file: 		ode.cpp
 * 	contents:   	main routine, sets up and runs ODE solver, then outputs results
 * 
 ****************************************************************************************************/


#include <ode.h>
#include <boost/numeric/odeint.hpp>
#include <rhs.h>
#include <params.h>
#include <frg.h>
#include <output.h>
#include <observables.h>
#include <state.h>
#include <observer.h>
#include <H5Tools.h>

using namespace boost::numeric::odeint;
using namespace H5; 
using namespace std; 

int main ( int argc, char * argv[])
{
   // ------ Set output format

   cout << scientific << setprecision(4); 

   // ------ Read in parameters from file

   if( argc == readIn_lst.size() + 1 )
   {
      for( int i = 0; i < readIn_lst.size(); i++ )
	 *readIn_lst[i] = atof( argv[i+1]);
      update_dep_params(); 
   }
   else
   {
      cout << "Wrong amount of arguments, using defaults." << endl << endl; 
   }

   // ------ Initialize rhs object and observer list

   cout << " Initializing ODE solver... " << endl << endl;

   rhs_t rhs;
   
   //vector< pair< string, vector<double> >  >   flow_obs_lst;		///< Vector to contain the values of observables tracked during the flow

#ifdef READIN
   // ------ Reading exact data as specified in fRG.cpp
   cout << "Reading from file:" << endl;
   read_exact(); 
#endif

   // ------ Write initial values

   cout << " Writing initial conditions... " << endl << endl; 
   state_t state_vec; 
   
   // Write initial values for Sig, phi's, P's and chi's
   state_vec.gf_Sig().init( Sig_init ); 

   state_vec.gf_phi_pp().init( phi_init ); 
   state_vec.gf_phi_ph().init( phi_init ); 
   state_vec.gf_phi_xph().init( phi_init ); 

   state_vec.gf_P_pp().init( P_init ); 
   state_vec.gf_P_ph().init( P_init ); 
   state_vec.gf_P_xph().init( P_init ); 

   state_vec.gf_chi_pp().init( chi_init ); 
   state_vec.gf_chi_ph().init( chi_init ); 
   state_vec.gf_chi_xph().init( chi_init ); 

   cout << " Starting scale-dependent PARQUET solver... "  << endl << endl;
 //bool success = true; 

   vector<double> Lam_list = { 1.0 }; 
   bool success; 

   for( double Lam: Lam_list )
   {
      cout << " Starting self-consistency cycle at scale " << Lam << endl << endl; 

      state_t state_vec_old;
      double diff; 
      double damping = 0.0; 
      int count = 0; 

      do
      {
         count++; 
	 cout << " Iteration " << count << endl; 
	 state_vec_old = state_vec; 

	 // --- Update state_vec
	 rhs( state_vec_old, state_vec, Lam ); 
	 cout << " ... Damping with factor " << damping << endl;
	 state_vec = state_vec_old * damping + state_vec * ( 1.0 - damping ); 
	 cout << " ... Calculating difference ... ";
	 diff = norm( state_vec - state_vec_old ); 
	 cout << diff << endl << endl; 

	 // Output intermediate files
	 if( count == 1 || count % 20 == 0 )
	    write_all( "log/iter_" + to_string(count) + ".h5", state_vec ); 

      } while( diff > 1e-12 && diff < MAX_COUPLING ); 

      if( diff >= MAX_COUPLING )
      {
	 cout << " DIVERGENT! " << endl << endl; 
	 success = false; 
      }
      else
      {
	 cout << " converged after " << count << " iterations " << endl << endl; 
	 success = true; 
      }

      cout << " Calculating phi functions by inverse ... " << endl;

      gf_1p_mat_t Gvec( POS_1P_RANGE ); 
      Gvec.init( bind( &state_t::GMat, boost::cref(state_vec), _1, Lam ) ); // Initialize big Green function vector 

      rhs_t::phi_pp_inverse( state_vec, state_vec.gf_phi_pp(), Gvec, LAM_FIN ); 
      rhs_t::phi_ph_xph_inverse( state_vec, state_vec.gf_phi_ph(), state_vec.gf_phi_xph(), Gvec, LAM_FIN ); 

   }

 // ------ Write output file

   H5std_string FILE_NAME("dat/dat");
   vector<pair<string, double>> fname_params = {{ "U", UINT }, { "Beta", BETA }, { "PFCB", POS_BFREQ_COUNT_CHI }}; //{ "PFC", POS_FFREQ_COUNT_P }}; //, { "Eps", EPS }};
   for( auto parStr: fname_params )
   {
      FILE_NAME.append("_" + parStr.first );
      string valStr = to_string( ( long double ) parStr.second );
      valStr.erase( valStr.find_last_not_of('0') + 1, string::npos ); // delete trailing zeros
      valStr.erase( valStr.find_last_not_of('.') + 1, string::npos ); // delete trailing dot
      replace( valStr.begin(), valStr.end(), '.', 'p');	// replace dot with p
      FILE_NAME.append( valStr ); 
   }

#ifdef NO_MOMENTA

#ifdef QMC_BATH
   FILE_NAME.append("_QMCB"); 
#elif ED_BATH
   FILE_NAME.append("_EDB"); 
#endif

#else
   FILE_NAME.append("_HUB"); 
#endif

   //FILE_NAME.append( "_" + FLOW_SCHEME_ABBREV ); 

   FILE_NAME.append("_PARQ"); 

#ifndef READIN
   FILE_NAME.append("_APPR"); 
#endif

   FILE_NAME.append("_SU2"); 

#ifdef METHOD2
   FILE_NAME.append("_METH2"); 
#elif METHOD1
   FILE_NAME.append("_METH1");
#endif
  
   if( !success ) FILE_NAME.append("_DIVERGENT"); 

   FILE_NAME.append(".h5"); 

   write_all( FILE_NAME, state_vec ); 

   return 0;

}

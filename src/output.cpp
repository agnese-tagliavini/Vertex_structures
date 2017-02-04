
/************************************************************************************************//**
 *  		
 * 	file: 		output.cpp
 * 	contents:   	See output.h
 * 
 ****************************************************************************************************/


#include <output.h>
#include <params.h>
#include <const.h>
#include <frg.h>
#include <observables.h>
#include <grid.h>
#include <ctime>
#include <H5Tools.h>

//using namespace boost;
using namespace H5; 
using namespace std; 


void init( H5File& file )
{
   time_t curr_time = time( nullptr );
   char* time_str = asctime( localtime(&curr_time )); 
      
   DataSpace dsp = DataSpace( H5S_SCALAR );

   StrType strdatatype( PredType::C_S1, H5T_VARIABLE ); // String type with variable length

   hid_t time_att_id = H5Acreate2( file.getId(), "DATE_TIME", strdatatype.getId(), dsp.getId(), H5P_DEFAULT, H5P_DEFAULT ); // Adds attribute to root group of file - find c++ equivalent! 

   Attribute time_att( time_att_id );
   time_att.write( strdatatype, &time_str );
}

void write_config( H5File& file )
{
   Group group( file.createGroup("/Config"));

   vector<pair<string, double>> config_scalar_list = 
   { 
      { "POS_FFREQ_COUNT_SIG", POS_FFREQ_COUNT_SIG }, 
      { "POS_FFREQ_COUNT_PHI", POS_FFREQ_COUNT_PHI }, 
      { "POS_BFREQ_COUNT_PHI", POS_BFREQ_COUNT_PHI }, 
      { "POS_FFREQ_COUNT_P", POS_FFREQ_COUNT_P }, 
      { "POS_BFREQ_COUNT_P", POS_BFREQ_COUNT_P }, 
      { "POS_BFREQ_COUNT_CHI", POS_BFREQ_COUNT_CHI }, 
      { "POS_INT_RANGE", POS_INT_RANGE }, 
      { "PATCH_COUNT", PATCH_COUNT }, 
      { "QN_COUNT", QN_COUNT }, 
      { "LAM_START", LAM_START }, 
      { "LAM_FIN", LAM_FIN }, 
      { "INIT_STEP", INIT_STEP }, 
      { "ERR_ABS", ERR_ABS }, 
      { "ERR_REL", ERR_REL }, 
      { "MAX_COUPLING", MAX_COUPLING }
   };

   for( auto conf : config_scalar_list )
      write( conf.second, group, conf.first); 
   
   vector<pair<string, string>> config_text_list = 
   {
      { "ERR_STEPPER", ERR_STEPPER_STRING }, 
      { "FLOW_SCHEME", FLOW_SCHEME_STRING }, 
      { "BATH_DOS", BATH_DOS_STRING }
   }; 
   
   for( auto conf : config_text_list )
      write( conf.second, group, conf.first); 

   DataSpace config_dsp = DataSpace ( H5S_SCALAR );

   int katanin_int =
#ifdef KATANIN
      1;
#else
      0; 
#endif
   write( katanin_int, group, "KATANIN" ); 
   
   int forced_zeros_int =
#ifdef FORCED_ZEROS
      1;
#else
      0; 
#endif
   write( forced_zeros_int, group, "FORCED_ZEROS" ); 
}

void write_params( H5File& file )
{

   Group group( file.createGroup("/Params"));

#ifdef NO_MOMENTA
   vector<pair<string, double>> par_lst = 
   { 
      { "UINT", UINT }, 
      { "BETA", BETA },
      { "B", B }, 
      { "GAM_L", GAM_L }, 
      { "DEL", DEL }, 
      { "EPS", EPS }, 
      { "PHI", PHI }
   };
#else
   vector<pair<string, double>> par_lst = 
   { 
      { "UINT", UINT }, 
      { "BETA", BETA },
      { "B", B }, 
      { "MU", MU }, 
      { "T_PRIME", T_PRIME } 
   };
#endif

   for( auto par : par_lst )
      write( par.second, group, par.first); 

}

void write_vert_tensor( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/Vert") );

   gf_2p_t vert_plot( POS_PLOT_RANGE_VERT ); 
   vert_plot.init( boost::bind( &state_t::vertx, boost::ref(state_vec), _1 ) ); 

   write( vert_plot, group, "" ); 
   write( F_Grid( POS_PLOT_RANGE_VERT, 2.0*PI / BETA ), group );
}

void write_lambda_tensor( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/Lambda") );

   gf_2p_t lambda_plot( POS_PLOT_RANGE_VERT ); 
   lambda_plot.init( boost::bind( &state_t::lambda, boost::ref(state_vec), _1 ) ); 

   write( lambda_plot, group, "" ); 
   write( F_Grid( POS_PLOT_RANGE_VERT, 2.0*PI / BETA ), group );
}

void write_Sig_tensor( H5::H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/Sig") );
   write( state_vec.gf_Sig(), group ); 
   write( F_Grid( POS_FFREQ_COUNT_SIG, 2.0*PI / BETA ), group );
}

//void write_vert_func( H5File& file, const state_t& state_vec )
//{
//   Group group( file.createGroup("/vert_func") );
//
//   gf_phi_t gf_vert_pp_plot( POS_PLOT_RANGE_PHI ); 
//   gf_phi_t gf_vert_ph_plot( POS_PLOT_RANGE_PHI ); 
//   gf_phi_t gf_vert_xph_plot( POS_PLOT_RANGE_PHI ); 
//
//   gf_vert_pp_plot.init( bind( &state_t::vertx_pp, boost::cref(state_vec), _1 ) ); 
//   gf_vert_ph_plot.init( bind( &state_t::vertx_ph, boost::cref(state_vec), _1 ) ); 
//   gf_vert_xph_plot.init(bind( &state_t::vertx_xph, boost::cref(state_vec), _1 ) ); 
//
//   write( gf_vert_pp_plot, group, "_PP" ); 
//   write( gf_vert_ph_plot, group, "_PH" ); 
//   write( gf_vert_xph_plot, group, "_XPH" ); 
//
//   write( Bos_Grid( POS_BFREQ_COUNT_PHI, 2.0*PI / BETA ), group );
//   write( F_Grid( POS_PLOT_RANGE_PHI, 2.0*PI / BETA ), group );
//}
void write_phi_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/phi_func") );

   gf_phi_t gf_phi_pp_plot( POS_PLOT_RANGE_PHI ); 
   gf_phi_t gf_phi_ph_plot( POS_PLOT_RANGE_PHI ); 
   gf_phi_t gf_phi_xph_plot( POS_PLOT_RANGE_PHI ); 

   gf_phi_pp_plot.init( bind( &state_t::phi_pp, boost::cref(state_vec), _1 ) ); 
   gf_phi_ph_plot.init( bind( &state_t::phi_ph, boost::cref(state_vec), _1 ) ); 
   gf_phi_xph_plot.init( bind( &state_t::phi_xph, boost::cref(state_vec), _1 ) ); 

   write( gf_phi_pp_plot, group, "_PP" ); 
   write( gf_phi_ph_plot, group, "_PH" ); 
   write( gf_phi_xph_plot, group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_PHI, 2.0*PI / BETA ), group );
   write( F_Grid( POS_PLOT_RANGE_PHI, 2.0*PI / BETA ), group );
}

void write_chi_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/chi_func") );

   write( state_vec.gf_chi_pp(), group, "_PP" ); 
   write( state_vec.gf_chi_ph(), group, "_PH" ); 
   write( state_vec.gf_chi_xph(), group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_CHI, 2.0*PI / BETA ), group );
}

void write_P_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/P_func") );

   write( state_vec.gf_P_pp(), group, "_PP" ); 
   write( state_vec.gf_P_ph(), group, "_PH" ); 
   write( state_vec.gf_P_xph(), group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_P, 2.0*PI / BETA ), group );
   write( F_Grid( POS_FFREQ_COUNT_P, 2.0*PI / BETA ), group );
}

void write_R_func( H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/R_func") );

   gf_phi_t R_pp_plot( POS_PLOT_RANGE_PHI ); 
   gf_phi_t R_ph_plot( POS_PLOT_RANGE_PHI ); 
   gf_phi_t R_xph_plot( POS_PLOT_RANGE_PHI ); 

   R_pp_plot.init( boost::bind( &state_t::R_pp, boost::ref(state_vec), _1 ) ); 
   R_ph_plot.init( boost::bind( &state_t::R_ph, boost::ref(state_vec), _1 ) ); 
   R_xph_plot.init( boost::bind( &state_t::R_xph, boost::ref(state_vec), _1 ) ); 

   write( R_pp_plot, group, "_PP" ); 
   write( R_ph_plot, group, "_PH" ); 
   write( R_xph_plot, group, "_XPH" ); 

   write( Bos_Grid( POS_BFREQ_COUNT_PHI, 2.0*PI / BETA ), group );
   write( F_Grid( POS_PLOT_RANGE_PHI, 2.0*PI / BETA ), group );
}

void write_Giw_tensor( H5::H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/Giw") );

   gf_1p_t Gvec(POS_1P_RANGE);
   Gvec.init( bind( &state_t::Gval, boost::cref(state_vec), _1, LAM_FIN ) ); // Initialize big Green function vector 

   write( Gvec, group ); 
   write( F_Grid( POS_1P_RANGE, 2.0 * PI / BETA ), group );  
}

void write_observables( H5::H5File& file, const state_t& state_vec )
{
   Group group( file.createGroup("/Observables") );

   vector<pair<string, double>> obs_lst = 
   { 
      { "EFF_MASS", eff_mass( state_vec ) }, 
      { "ERR_EFF_MASS", err_eff_mass( state_vec ) }, 
   };

   // Add all possible fillings
   for( int s_in = 0; s_in < QN_COUNT; s_in++)
      for( int s_out = 0; s_out < QN_COUNT; s_out++)
	 obs_lst.push_back( { "FILLING_" + to_string( ( long long ) s_out ) + to_string( ( long long ) s_in ), filling( s_in, s_out, state_vec ) } );
   
   for( auto obs : obs_lst )
      write( obs.second, group, obs.first ); 
}

//void write_flow_observables( vector<  pair< string, vector<double> >  >&   flow_obs_lst, H5::H5File& file )
//{
   //Group group( file.createGroup("/Flow_obs") );

   //for( auto flow_obs : flow_obs_lst )
      //write( flow_obs.second, group, flow_obs.first ); 
//}

void write_all( std::string file_name, const state_t& state_vec )
{
   try
   {
      H5File file( file_name, H5F_ACC_TRUNC );
      init( file );		
      write_config( file );		
      write_params( file );
      write_Sig_tensor( file, state_vec );
      write_vert_tensor( file, state_vec ); 
      write_lambda_tensor( file, state_vec ); 
//      write_vert_func( file, state_vec );
      write_phi_func( file, state_vec ); 
      write_chi_func( file, state_vec ); 
      write_P_func( file, state_vec ); 
      write_R_func( file, state_vec ); 
      write_Giw_tensor( file, state_vec );
      write_observables( file, state_vec );
      file.flush( H5F_SCOPE_LOCAL ); 
      file.close(); 
   }
   // catch failure caused by the H5File operations
   catch( FileIException error )
   {
      error.printError();
   }
   // catch failure caused by the DataSet operations
   catch( DataSetIException error )
   {
      error.printError();
   }
   // catch failure caused by the DataSpace operations
   catch( DataSpaceIException error )
   {
      error.printError();
   }
   // catch failure caused by the DataSpace operations
   catch( DataTypeIException error )
   {
      error.printError();
   }
}

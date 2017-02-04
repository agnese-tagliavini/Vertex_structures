
/************************************************************************************************//**
 *  		
 * 	file: 		observer.cpp
 * 	contents:   	See observer.h
 * 
 ****************************************************************************************************/


#include <observer.h>
#include <observables.h>

//using namespace std;

bool observer_t::operator() ( const state_t& state_vec, const double Lam )
{

   // Write the observables to be tracked during the flow
   auto flow_obs_lst_itr = flow_obs_lst.begin(); 
   ( *flow_obs_lst_itr ).second.push_back( Lam ); 
   
   // Calculate maximal coupling
   double max_cpl = get_max_cpl( state_vec ); 
   ++flow_obs_lst_itr; 
   ( *flow_obs_lst_itr ).second.push_back( max_cpl ); 

   // Calculate further observables to be tracked
   for( auto flow_obs_func : flow_obs_func_lst )
   {
      ++flow_obs_lst_itr; 
      ( *flow_obs_lst_itr ).second.push_back( flow_obs_func.second( state_vec ) ); 
   }

   return ( max_cpl < MAX_COUPLING ) ? 1 : 0; 

}

observer_t::observer_t( std::vector<   std::pair< std::string, std::vector<double> >    >&   flow_obs_lst_ ): flow_obs_lst( flow_obs_lst_ )
{

   // Track lambda, abs_max_cpl 
   flow_obs_lst.push_back( { "LAM", {} } ); 
   flow_obs_lst.push_back( { "ABS_MAX_CPL", {} } ); 

   // Fill rest of flow observables into list, COMPARE observables.cpp
   for( auto flow_obs_func : flow_obs_func_lst )
      flow_obs_lst.push_back( { flow_obs_func.first, {} } ); 

}


/*******************************************************************************************//** @file
 *  		
 * 	file: 		observer.h
 * 	contents:  	functor to track observables during the flow
 * 
 ****************************************************************************************************/


#pragma once

#include <state.h>

class observer_t 												///< Functor to specify the rhs calculation for both the self-energy and the vertex
{
   public:
      bool operator() ( const state_t& state_vec, const double Lam ); 								///< Overload call operator to calculate full rhs
      observer_t( std::vector<   std::pair< std::string, std::vector<double> >    >&   flow_obs_lst_ );		///< Constructor taking the reference of an rhs object

   private:
      std::vector<   std::pair< std::string, std::vector<double> >    >&   flow_obs_lst;					///< Vector to contain the values of observables tracked during the flow
};

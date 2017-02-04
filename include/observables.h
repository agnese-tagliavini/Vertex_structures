
/*******************************************************************************************//** @file
 *  		
 * 	file: 		observalbes.h
 * 	contents:  	Define functions to calculate observables from results
 * 
 ****************************************************************************************************/


#pragma once

#include <state.h>

typedef double (*obs_func_t )( const state_t& ); 	///< Symmetry function that acts on an idx_2p_t object and alters it 
typedef std::vector< std::pair<  std::string, obs_func_t > >	flow_obs_func_lst_t; 

extern flow_obs_func_lst_t flow_obs_func_lst; 	///< List to contain all observable functions to be tracked during the flow

double eff_mass( const state_t& state_vec ); 		///< Return effective mass based on self-energy slope at frequency zero
double err_eff_mass( const state_t& state_vec ); 	///< Return an error estimate on the effective mass by using the first and third matsubara frequency
double filling( int s_in, int s_out, const state_t& state_vec ); 	///< Calculate the average occupation number summing the diagonal part of all quantum numbers
double jos_curr( const state_t& state_vec ); 	///< Calculate the average occupation number summing the diagonal part of all quantum numbers
double get_max_cpl( const state_t& state_vec ); 	///< Find the coupling with the largest absolute value. Returns the absolute value and the position

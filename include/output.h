
/*******************************************************************************************//** @file
 *  		
 * 	file: 		output.h
 * 	contents:  	Functions to write output to hdf5 file
 * 
 ****************************************************************************************************/


#pragma once

#include <H5Cpp.h>
#include <rhs.h>
#include <string>

void init( H5::H5File& file );			///< Initialize with time and date
void write_config( H5::H5File& file );		///< Write configuration to output file
void write_params( H5::H5File& file );		///< Write parameters to output file
//void write_vert_tensor( H5::H5File& file, const state_t& state_vec );  	///<	Write 1PI two-particle vertex to output file
void write_vert_func( H5::H5File& file, const state_t& state_vec ); 		///< 	Write phi functions to output file
void write_genchi_func( H5::H5File& file, const state_t& state_vec ); 		///< 	Write phi functions to output file
void write_lambda_tensor( H5::H5File& file, const state_t& state_vec );  	///<	Write 2PI two-particle vertex to output file
void write_Sig_tensor( H5::H5File& file, const state_t& state_vec );		///<	Write self-energy to output file
void write_phi_func( H5::H5File& file, const state_t& state_vec ); 		///< 	Write phi functions to output file
void write_chi_func( H5::H5File& file, const state_t& state_vec );  	///<	Write Karrasch functions to output file
void write_P_func( H5::H5File& file, const state_t& state_vec ); 		///<	Write P functions to output file
void write_R_func( H5::H5File& file, const state_t& state_vec ); 		///<	Write R functions to output file
void write_Giw_tensor( H5::H5File& file, const state_t& state_vec );	///<	Write Green function to output file
void write_observables( H5::H5File& file, const state_t& state_vec );	///<	Write observables to output file
//void write_flow_observables( std::vector<   std::pair< std::string, std::vector<double> >    >&   flow_obs_lst, H5::H5File& file );	///<	Write observables tracked during the flow

void write_all( std::string file_name, const state_t& state_vec ); 

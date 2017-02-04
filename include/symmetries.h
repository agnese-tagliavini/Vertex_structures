
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symmetries.h
 * 	contents:  	List of symmetry functions
 * 
 ****************************************************************************************************/

#pragma once

#include <def.h>
#include <grid.h>
#include <symmetry_group.h>
#include <translate.h>

// Symmetries - Two-particle ( Copied from Vertex code )
template<typename notation> operation exch_in( idx_2p_t& idx );			///< Exchange ingoing lines
template<typename notation> operation exch_out( idx_2p_t& idx );		///< Exchange outgoing lines
template<typename notation> operation compl_conj( idx_2p_t& idx );		///< Complex conjugation 
template<typename notation> operation time_rev( idx_2p_t& idx );		///< Time reversal symmetry
template<typename notation> operation particle_hole( idx_2p_t& idx );		///< Particle hole symmetry in nambu notation
template<typename notation> operation spin_symm( idx_2p_t& idx );		///< Spin symmetry in nambu notation

template<typename notation> operation rot_k( idx_2p_t& idx );			///< Rotate all momenta by 90 degrees
template<typename notation> operation mirror_vert( idx_2p_t& idx );		///< Mirror all momenta vertically
template<typename notation> operation mirror_diag( idx_2p_t& idx );		///< Mirror all momenta diagonally

// Translation of Two-particle symmetries to P ( problematic )
using symm_func_t = operation (*)( idx_2p_t& idx ); 
template<symm_func_t symm_func> operation symm_P( idx_P_t& idx )
{
   idx_phi_t idx_phi = iP_to_iphi( idx ); 
   operation op = symm_func( idx_phi ); 
   idx = iphi_to_iP( idx_phi ); 
}

// symmetries from symmetries.cpp , lattice/system independent
// Two-particle - Mixed Notation ( analytic derivation )
operation hmirror_phi_pp( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_phi_ph( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_phi_xph( idx_phi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_phi_pp( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_phi_ph( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_phi_xph( idx_phi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_phi_pp( idx_phi_t& idx);	///< time reversal
operation timerev_phi_ph( idx_phi_t& idx);	///< time reversal
operation timerev_phi_xph( idx_phi_t& idx);	///< time reversal

// Specific symmetries of P functions - Observed
operation hmirror_P_pp( idx_P_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_timerev_P_ph( idx_P_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_timerev_P_pp( idx_P_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation vmirror_P_ph( idx_P_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_P_xph( idx_P_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_P_xph( idx_P_t& idx);	///< time reversal

// Specific symmetries of chi functions - Observed

operation hmirror_chi_ph( idx_chi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation hmirror_chi_xph( idx_chi_t& idx );	///< horizontal mirroring - exchange of both incoming/outcoming legs
operation cconj_chi_pp( idx_chi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_chi_ph( idx_chi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation cconj_chi_xph( idx_chi_t& idx );	///< complex conjugate -exchange incoming with outgoing
operation timerev_chi_xph( idx_chi_t& idx);	///< time reversal

// Specific symmetries of Self energy
operation cconj_sig( idx_1p_t& idx );	///< Change sign of bosonic frequency + compl conj

// -- Tools
inline void freq_sign_change( int& idx )
{
   idx = - idx - 1;
}

inline void flip_spin( int& idx )
{   
   idx = !idx;
}

inline void mirror_mom_vert( int& idx )
{
   idx = mirror_mom_vert_arr[idx];
}

inline void mirror_mom_diag( int& idx )
{
   idx = mirror_mom_diag_arr[idx];
}

inline void mirror_mom_pipi( int& idx )
{   
   idx = mirror_mom_pipi_arr[idx]; 
}

#include <symmetries_impl.h> 	// contains implementations of template functions


/******************************************************************************************//** @file
 *  		
 * 	file: 		def.h
 * 	contents:  	Definition of used correlation function containers ( wrapper around gf container )
 * 
 ****************************************************************************************************/

#pragma once

#define EIGEN_STACK_ALLOCATION_LIMIT 0
#include <complex>
#include <Eigen/Core>
#include <vector>
#include <boost/multi_array.hpp>

#include <const.h>
#include <gf.h>

using dcomplex = std::complex<double>;						///< Complex double type
using MatQN = Eigen::Matrix<dcomplex, QN_COUNT, QN_COUNT, Eigen::RowMajor>;	///< Complex matrix representing the discrete quantum number structure
using MatQNQN = Eigen::Matrix<dcomplex, QN_COUNT*QN_COUNT, QN_COUNT*QN_COUNT, Eigen::RowMajor>;	///< Complex matrix representing the discrete quantum number structure of two-particle function

using MapXcd = Eigen::Map< Eigen::Matrix< dcomplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >  >; 
using MapXd = Eigen::Map< Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >  >; 
using MapXdvec = Eigen::Map< Eigen::Matrix< double, Eigen::Dynamic, 1 > >; 

using MapArrXcd = Eigen::Map< Eigen::Array< dcomplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >  >; 
using MapArrXd =  Eigen::Map< Eigen::Array< double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor >  >; 

#define INSERT_COPY_AND_ASSIGN(X) 					\
X( const X & gf_obj ):    						\
   base_t( gf_obj )							\
{}       								\
X( X && gf_obj ):							\
   base_t( std::move(gf_obj) )						\
{}      								\
X & operator=( const X & gf_obj )					\
{									\
   base_t::operator=( gf_obj ); 					\
   return *this; 							\
} 									\
X & operator=( X && gf_obj )						\
{									\
   base_t::operator=( std::move( gf_obj) ); 				\
   return *this; 							\
} 


// Container and index types
class gf_1p_mat_t : public gf< MatQN, 2 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< MatQN, 2 >; 

      gf_1p_mat_t( int pos_freq_count_= POS_FFREQ_COUNT_SIG ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][PATCH_COUNT] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_mat_t)
}; 
using idx_1p_mat_t = gf_1p_mat_t::idx_t;   

enum class I1P{ w, k, s_in, s_out }; 
class gf_1p_t : public  gf< dcomplex, 4 > 		///< Container type for one-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_1p_t( int pos_freq_count_ = POS_FFREQ_COUNT_SIG ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][PATCH_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_t)
}; 
using idx_1p_t = gf_1p_t::idx_t; 

enum class I2P{ w1_in, w2_in, w1_out, k1_in, k2_in, k1_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_2p_t : public gf< dcomplex, 10 > 		///< Container type for two-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 10 >; 

      gf_2p_t( int pos_freq_count_ = POS_FFREQ_COUNT_PHI ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][ffreq(pos_freq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][PATCH_COUNT][PATCH_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_2p_t)
}; 
using idx_2p_t = gf_2p_t::idx_t; 

enum class IPHI{ W, w_in, w_out, K, k_in, k_out, s1_in, s2_in, s1_out, s2_out }; 
class gf_phi_t : public gf< dcomplex, 10 > 		///< Container type for two-particle correlation functions, holds ind_cpl_t
{
   public:
      using base_t = gf< dcomplex, 10 >;

      gf_phi_t( int pos_freq_count_ = POS_FFREQ_COUNT_PHI, int pos_bfreq_count_ = POS_BFREQ_COUNT_PHI ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][PATCH_COUNT][PATCH_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_phi_t)
}; 
using idx_phi_t = gf_phi_t::idx_t; 

enum class IP{ W, w, K, k, s1_in, s2_in, s1_out, s2_out }; 
class gf_P_t : public gf< dcomplex, 8 > 		///< Container type for two-particle correlation functions, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 8 >; 

      gf_P_t(int pos_freq_count_ = POS_FFREQ_COUNT_P, int pos_bfreq_count_ = POS_BFREQ_COUNT_P ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][ffreq(pos_freq_count_)]
	       [PATCH_COUNT][PATCH_COUNT]
	       [QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_P_t)
}; 
using idx_P_t = gf_P_t::idx_t; 

enum class ICHI{ W, K, s1_in, s2_in, s1_out, s2_out }; 
class gf_chi_t : public gf< dcomplex, 6 > 		///< Container type for two-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 6 >; 

      gf_chi_t( int pos_bfreq_count = POS_BFREQ_COUNT_CHI ):
	 base_t( boost::extents[bfreq(pos_bfreq_count)][PATCH_COUNT][QN_COUNT][QN_COUNT][QN_COUNT][QN_COUNT] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_chi_t)
}; 
using idx_chi_t = gf_chi_t::idx_t; 

#ifdef FORCED_ZEROS
bool forced_zero_check( const idx_2p_t& idx ); 	///< Check if coupling should be artificially forced to zero
bool forced_zero_check( int s1_in, int s2_in, int s1_out, int s2_out ); 
#endif

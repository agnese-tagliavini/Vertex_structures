
/******************************************************************************************//** @file
 *  		
 * 	file: 		mymath.h
 * 	contents:  	Contains useful mathematical functions / routines
 * 
 ****************************************************************************************************/

#pragma once

#include <def.h>
#include <gsl/gsl_errno.h>
#include <complex>
#include <cmath>
#include<boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace boost::numeric::ublas;

// (m * t) [i] [j] = m [i] [j] * t

template<class E1, class T2> typename matrix_binary_scalar2_traits<E1, T2, scalar_multiplies<typename E1::value_type, T2> >::result_type operator * (const matrix_expression<E1> &e1,const T2 &e2);

// Matrix-matrix multiplication

//template<class T1, class E1, class T2, class E2>
//struct matrix_matrix_binary_traits {
//   typedef unknown_orientation_tag dispatch_category;
//   typedef typename promote_traits<T1, T2>::promote_type promote_type;
//   typedef matrix_matrix_binary<typename E1::const_closure_type, typename E2::const_closure_type, matrix_matrix_prod<T1, T2, promote_type> > expression_type;
//   typedef expression_type result_type;
//};

template<class E1, class E2> typename matrix_matrix_binary_traits<typename E1::value_type, E1, typename E2::value_type, E2>::result_type prod (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2, unknown_orientation_tag);

// Dispatcher

template<class E1, class E2> typename matrix_matrix_binary_traits<typename type_traits<typename E1::value_type>::precision_type, E1, typename type_traits<typename E2::value_type>::precision_type, E2>::result_type prec_prod (const matrix_expression<E1> &e1, const matrix_expression<E2> &e2, unknown_orientation_tag);

//void multiplyMatrices(dcomplex firstMatrix[][], dcomplex secondMatrix[][], dcomplex mult[][], int rowFirst, int columnFirst, int rowSecond, int columnSecond);
// ---- Useful functions

inline int div2_ceil( int W )
{
   return (W - 1000000)/2 + 500000; 
}

inline int div2_floor( int W )
{
   return (W + 1000000)/2 - 500000; 
}

double Theta( double x );			///< Theta step function
double sgn( double x );			///< Sign function with sgn( 0.0 )==0
double kdel( double a, double b );	///< Kronecker delta for two double values


// --- Self Consistency ---

typedef void(*sc_Func )( const std::complex<double>[], std::complex<double>[], void * q ); ///< General function for self_consistency
double dist( const dcomplex a[], const dcomplex b[], int len );	///< Arbitrary distance measure on array of length len
/**
 * 	Perform self consistency loop using general function type defined above, till the distance measure dist(..) between two
 * 	steps fullfills the errorbound err	
 */
int self_consistency( std::complex<double> init[], std::complex<double> fin[], int len, sc_Func fw, double err, void * q );


// --- ARRAY INTEGRATION -------

/**
 *	Trapezoidal integration of function y( x ) with array of x and array of y values as input ( both of length len )
 *	Nonequidistant grid for x possible
 */
std::complex<double> arrInt( double x[], std::complex<double> y[], int len );
/**
 *	Trapezoidal integration of function y( x ) with array of y values as input
 *	Equidistant grid with spacing h for x grid assumed
 */
std::complex<double> arrInt( double h, std::complex<double> y[], int len );
/**
 *	Extended Simpson integration of function y( x ) with array of x and array of y values as input ( both of length len )
 *	Nonequidistant grid for x possible
 */
std::complex<double> arrIntSimpson( double h, std::complex<double> y[], int len );
/**
 *	Extended Simpson integration of function y( x ) with array of y values as input
 *	Equidistant grid with spacing h for x grid assumed
 */
std::complex<double> arrIntSimpson( double h, std::complex<double> y[], int lenx, int leny );




// --- Integrations using ODE solver ---------

typedef int(*intFuncP )( double, const double[], double[], void*);
typedef std::complex<double>(*intFuncPCPLX )( double, void*);

double integrate( intFuncP fw, void* q, double err );
double integrate( intFuncP fw, double err );

std::complex<double> integrateCPLX( intFuncPCPLX fw, void* q, double err );
int realFunc( double w, const double y[], double g[], void* q );
int imagFunc( double w, const double y[], double g[], void* q );

struct intCPLXList{
    intFuncPCPLX fw;
    void* q;
    intCPLXList( intFuncPCPLX fw_val, void* q_val ){fw = fw_val; q = q_val;}
};


// --- INTEGRATION QAGS ---------

typedef double(*intFuncP_QAGS )( double, void*);

double integrate_QAGS( intFuncP_QAGS fw, void* q, double err );
double integrate_QAGS( intFuncP_QAGS fw, double err );

std::complex<double> integrateCPLX( intFuncPCPLX fw, void* q, double err );
double realFunc_QAGS( double w, void* q );
double imagFunc_QAGS( double w, void* q );

// -- Oscilatory

double integrate_QAGS_COS( intFuncP_QAGS fw, void* q, double err );
double integrate_QAGS_SIN( intFuncP_QAGS fw, void* q, double err );

// generate_sum_function provides the fitting of a finite matsubara sum in order to get the extension of the sum up to infinity
// The Least-Square fitting uses the following fitting function : sum_{i=0}^{fit_order} a_i (x)^{-i}
// This requires the minimization of the chisquare which has been done in a semi-analytic way (see notes by Georg Roehringer)
// generate_sum_func returns sumfit of a given single-variable function 
std::vector<double> generate_tail_weights( int iMin, int tail_length, int fit_order ); 
gf<double, 1> generate_weights( int iMin, int tail_length, int fit_order ); 
gf<double, 2> generate_2d_weights( int iMin, int tail_length, int fit_order ); 

template< typename value_t >
boost::function< value_t ( boost::function< value_t ( int i ) > ) > generate_sum_func( int iMin, int tail_length, int fit_order )
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   return [tail_weights,iMin,tail_length]( boost::function< value_t ( int i ) > func )
   {
      value_t val; 
      for( int i = -iMin; i <= iMin; ++i )
	 val += func(i); 

      for( int i = 1; i < tail_length; ++i )
	 val += ( func(iMin + i) + func(-iMin -i) ) * tail_weights[i]; 

      return val; 
   }; 
}

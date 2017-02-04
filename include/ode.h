
/*******************************************************************************************//** @file
 *  		
 * 	file: 		ode.h
 * 	contents:  	Adaptive stepper implementation that aborts for large coupling ( see const.h )
 * 			Also define norm for state_t
 * 
 ****************************************************************************************************/


#pragma once

#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>
#include <boost/numeric/odeint.hpp>
#include <def.h>

namespace boost {
namespace numeric {
namespace odeint {
namespace detail {

 /*
 * Changed integrate adaptive to consider abort
 */
template< class Stepper, class System, class State, class Time, class Observer >
bool integrate_adaptive_check(
	Stepper stepper, System system, State &start_state,
	Time &start_time, Time end_time, Time &dt,
	Observer observer //, controlled_stepper_tag
)
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;

    const size_t max_attempts = 1000;
    const char *error_string = "Integrate adaptive : Maximal number of iterations reached. A step size could not be found.";
    size_t count = 0;
    while( less_with_sign( start_time, end_time, dt ) )
    {
	if( !obs( start_state, start_time ) )
	{
	   std::cout << std::endl << " !!!! Vertex exceeded MAX_COUPLING at scale " << start_time << ", ODE solver stopped !!!! " << std::endl << std::endl; 
	   return 0; 
	}	   
	if( less_with_sign( end_time, start_time + dt, dt ) )
	{
	    dt = end_time - start_time;
	}

	size_t trials = 0;
	controlled_step_result res = success;
	do
	{
	    res = stepper.try_step( system, start_state, start_time, dt );
	    ++trials;
	}
	while( ( res == fail ) && ( trials < max_attempts ) );
	if( trials == max_attempts ) throw std::overflow_error( error_string );

	++count;
    }
    obs( start_state, start_time );
    std::cout << " ODE solver finished with " << count << " integration steps." << std::endl << std::endl; 
    return 1;
}

} // Namespace detail

template< class Stepper, class System, class State, class Time, class Observer >
bool integrate_adaptive_check(
        Stepper stepper, System system, State &start_state,
        Time start_time, Time end_time, Time dt,
        Observer observer )
{
    return detail::integrate_adaptive_check(
            stepper, system, start_state,
            start_time, end_time, dt,
            observer ); //, typename Stepper::stepper_category() );
}

} // Namespace odeint
} // Namespace numeric
} // Namespace boost

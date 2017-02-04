
/******************************************************************************************//** @file
 *  		
 * 	file: 		def.cpp
 * 	contents:  	General typedefs and constants
 * 
 ****************************************************************************************************/


#include <def.h>

#ifdef FORCED_ZEROS
bool forced_zero_check( const idx_2p_t& idx )
{
   return false;    
}

bool forced_zero_check( int s1_in, int s2_in, int s1_out, int s2_out )
{
   return false;    
}
#endif

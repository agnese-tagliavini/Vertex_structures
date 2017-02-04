
/************************************************************************************************//**
 *  		
 * 	file: 		symmetries.cpp
 * 	contents:   	see symmetries.h
 * 
 ****************************************************************************************************/

#include <symmetries.h>
#include <grid.h>

// -- Phi
operation hmirror_phi_pp( idx_phi_t& idx )
{
   idx(0) *= 1; 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) *=1;
   idx(4) = dif_k(idx(3),idx(4));
   idx(5) = dif_k(idx(3),idx(5));
   swap(idx(6),idx(7));
   swap(idx(8),idx(9));
   return operation( false, false ); 
}

operation hmirror_phi_ph( idx_phi_t& idx )
{
   idx(0) *= -1;
   swap(idx(1),idx(2)); 
   idx(3) =neg_k(idx(3));
   swap(idx(4),idx(5));
   idx(4) = dif_k(idx(4),idx(3));
   idx(5) = dif_k(idx(5),idx(3));
   swap(idx(6),idx(7));
   swap(idx(8),idx(9));
   return operation( false, false ); 
}

operation hmirror_phi_xph( idx_phi_t& idx )
{
   idx(0) *= -1;
   swap(idx(1),idx(2)); 
   idx(3) = neg_k(idx(3));
   swap(idx(4),idx(5));
   idx(4) = dif_k(idx(4),idx(3));
   idx(5) = dif_k(idx(5),idx(3));
   swap(idx(6),idx(7));
   swap(idx(8),idx(9));
   return operation( false, false ); 
}

operation cconj_phi_pp( idx_phi_t& idx )
{
   idx(0) *= -1; 
   swap(idx(1),idx(2));
   idx(3) *=1;
   swap(idx(4),idx(5));
   idx(4) = dif_k(idx(3),idx(4));
   idx(5) = dif_k(idx(3),idx(5));
   swap(idx(6),idx(9));
   swap(idx(7),idx(8));
   return operation( false, true ); 
}

operation cconj_phi_ph( idx_phi_t& idx )
{
   idx(0) *= -1;
   swap(idx(1),idx(2)); 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) *=1;
   swap(idx(4),idx(5));
   swap(idx(6),idx(9));
   swap(idx(7),idx(8));
   return operation( false, true ); 
}

operation cconj_phi_xph( idx_phi_t& idx )
{
   idx(0) *= 1;
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) = neg_k(idx(3));
   idx(4) = dif_k(idx(4),idx(3));
   idx(5) = dif_k(idx(5),idx(3));
   
   swap(idx(6),idx(9));
   swap(idx(7),idx(8));
   return operation( false, true ); 
}

operation timerev_phi_pp( idx_phi_t& idx )
{
   idx(0) *= 1; 
   swap(idx(1),idx(2));
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) =-idx(2)-1-abs(idx(0)%2);
   idx(3) *=1;
   swap(idx(4),idx(5));
   idx(4) = dif_k(idx(3),idx(4));
   idx(5) = dif_k(idx(3),idx(5));
   swap(idx(6),idx(9));
   swap(idx(7),idx(8));
   return operation( false, false ); 
}

operation timerev_phi_ph( idx_phi_t& idx )
{
   idx(0) *= 1;
   swap(idx(1),idx(2)); 
   idx(3) *=1;
   swap(idx(4),idx(5));
   swap(idx(6),idx(9));
   swap(idx(7),idx(8));
   return operation( false, false ); 
}

operation timerev_phi_xph( idx_phi_t& idx )
{
   idx(0) *= -1;
   idx(1) *= 1;
   idx(2) *= 1;
   idx(3) = neg_k(idx(3));
   idx(4) = dif_k(idx(4),idx(3));
   idx(5) = dif_k(idx(5),idx(3));
   swap(idx(6),idx(9));
   swap(idx(7),idx(8));
   return operation( false, false ); 
}

// --- P
operation hmirror_P_pp( idx_P_t& idx )
{
   idx(0) *= 1; 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) *=1;
   idx(3) = dif_k(idx(2),idx(3));
   swap(idx(4),idx(5));
   swap(idx(6),idx(7));
   return operation( false, false ); 
}

operation cconj_timerev_P_pp( idx_P_t& idx )
{
   idx(0) *= -1; 
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) *=1;
   idx(3) = dif_k(idx(2),idx(3));
   return operation( false, true ); // Quantum numbers don't change after cconj and time reversal
}

operation hmirror_timerev_P_ph( idx_P_t& idx )
{
   idx(0) *= -1;
   idx(1) *= 1;
   idx(2) =neg_k(idx(2));
   idx(3) = dif_k(idx(3),idx(2));
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   return operation( false, false ); 
}

operation vmirror_P_ph( idx_P_t& idx )
{
   idx(0) *=  1;
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) = neg_k(idx(2));
   idx(3) = dif_k(idx(3),idx(2));
   swap(idx(4),idx(6));
   swap(idx(5),idx(7));
   return operation( false, true ); 
}

operation cconj_P_xph( idx_P_t& idx )
{
   idx(0) *= 1;
   idx(1) =-idx(1)-1-abs(idx(0)%2);
   idx(2) = neg_k(idx(2));
   idx(3) = dif_k(idx(3),idx(2));
   swap(idx(4),idx(7));
   swap(idx(5),idx(6));
   return operation( false, true ); 
}

operation timerev_P_xph( idx_P_t& idx )
{
   idx(0) *= -1;
   idx(1) *= 1;
   idx(2) = neg_k(idx(2));
   idx(3) = dif_k(idx(3),idx(2));
   swap(idx(4),idx(7));
   swap(idx(5),idx(6));
   return operation( false, false ); 
}

// --- Chi

operation hmirror_chi_ph( idx_chi_t& idx )
{
   idx(0) *= -1;
   idx(1) =neg_k(idx(1));
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   return operation( false, false ); 
}

operation hmirror_chi_xph( idx_chi_t& idx )
{
   idx(0) *= -1;
   idx(1) = neg_k(idx(1));
   swap(idx(2),idx(3));
   swap(idx(4),idx(5));
   return operation( false, false ); 
}

operation cconj_chi_pp( idx_chi_t& idx )
{
   idx(0) *= -1; 
   idx(1) *=  1;
   swap(idx(2),idx(5));
   swap(idx(3),idx(4));
   return operation( false, true ); 
}

operation cconj_chi_ph( idx_chi_t& idx )
{
   idx(0) *= -1;
   idx(1) *=  1;
   swap(idx(2),idx(5));
   swap(idx(3),idx(4));
   return operation( false, true ); 
}

operation cconj_chi_xph( idx_chi_t& idx )
{
   idx(0) *= 1;
   idx(1) = neg_k(idx(1));
   swap(idx(2),idx(5));
   swap(idx(3),idx(4));
   return operation( false, true ); 
}

operation timerev_chi_xph( idx_chi_t& idx )
{
   idx(0) *= -1;
   idx(1) = neg_k(idx(1));
   swap(idx(2),idx(5));
   swap(idx(3),idx(4));
   return operation( false, false ); 
}

// --- Sig
operation cconj_sig( idx_1p_t& idx )
{
   idx(0) = -idx(0) -1;
   idx(1) *= 1;  
   swap(idx(2),idx(3));
   return operation( false, true ); 
}

// =========== One particle ( Copied from old vertex code )

//operation spin_symm( idx_1p_t& idx )
//{
   //freq_sign_change( idx(0) );
   
   //swap( idx(2), idx(3) );

   //flip_spin( idx(2) );
   //flip_spin( idx(3) );

   //if ( ( idx(2) + !( idx(3) ) ) % 2  != 0 ) // if uneven amount of creation operators minus sign, see masterthesis
      //return operation( true, false );

   //return operation( false, false );
//}

//operation compl_conj( idx_1p_t& idx )
//{
   //freq_sign_change( idx.w );

   //// Changing momenta sign is unnecessary since equal to rotating twice by 90 degrees

   //swap( idx.s_in, idx.s_out );
   //return operation( false, true );
//}

//operation time_rev( idx_1p_t& idx )
//{
   //swap( idx.s_in, idx.s_out );
   //return operation( false, false );
//}

//operation particle_hole( idx_1p_t& idx )
//{
//#ifndef NO_MOMENTA
   //mirror_mom_pipi( idx.k );
//#endif
   //return operation( false, true );
//}

//operation spin_symm( idx_1p_t& idx )
//{
   //freq_sign_change( idx.w );
   
   //swap( idx.s_in, idx.s_out );

   //flip_spin( idx.s_in );
   //flip_spin( idx.s_out );

   //if ( ( idx.s_in + !( idx.s_out ) ) % 2  != 0 ) // if uneven amount of creation operators minus sign, see masterthesis
      //return operation( true, false );

   //return operation( false, false );
//}

//operation rot_k( idx_1p_t& idx )
//{
   //idx.k = rot_k_idx_arr[idx.k];
   //return operation( false, false );
//}

//operation mirror_vert( idx_1p_t& idx )
//{
   //mirror_mom_vert( idx.k );
   //return operation( false, false );
//}


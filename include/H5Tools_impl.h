
/*******************************************************************************************//** @file
 *  		
 * 	file: 		H5Tools_impl.h
 * 	contents:  	Implementation of templated hdf5 write-functions
 * 
 ****************************************************************************************************/

#include <cassert>
#include <iostream>
#include <def.h>

template< size_t ndims > void write( const boost::multi_array<double, ndims>& my_arr, H5::Group& group, const std::string& dataset_name )
{
   using namespace H5; 

   // Establish H5 DataSpace from my_arr 
   auto& dims = reinterpret_cast<const boost::array<hsize_t, ndims>&>(*my_arr.shape()); 
   hsize_t hndims(ndims); 		
   DataSpace dataspace( hndims, dims.data() );

   // Modify dataset creation with property list
   DSetCreatPropList plist;

#ifdef COMPRESS
   std::vector<hsize_t> chunk_dims;	// chunk dimensions
   for( int i = 0; i < ndims; ++i )
      chunk_dims.push_back( std::min( int(dims[i]) , 20 ) ); 

   plist.setChunk( ndims, chunk_dims.data() );

   // Set ZLIB (DEFLATE) Compression using level 6.
   // To use SZIP compression comment out this line.
   plist.setDeflate(6);

   // Uncomment these lines to set SZIP Compression
   //unsigned szip_options_mask = H5_SZIP_NN_OPTION_MASK;
   //unsigned szip_pixels_per_block = 16;
   //plist.setSzip(szip_options_mask, szip_pixels_per_block);
#endif
		     
   DataSet dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, dataspace, plist );

   dataset.write( my_arr.data(), PredType::NATIVE_DOUBLE );
}

template< size_t ndims > void write( const boost::multi_array<dcomplex, ndims>& my_arr, H5::Group& group, const std::string& dataset_name )
{
   using namespace std; 
   
   // Create array for real and imaginary part with equal size
   boost::multi_array<double, ndims> arr_reals, arr_imags; 
   auto& dims = reinterpret_cast<const boost::array<size_t, ndims>&>(*my_arr.shape()); 
   arr_reals.resize( dims ); 
   arr_imags.resize( dims ); 

   // Copy real and imaginary part into arrays
   transform( my_arr.data(), my_arr.data() + my_arr.num_elements(), arr_reals.data(), [](dcomplex a){ return real(a); } ); 
   transform( my_arr.data(), my_arr.data() + my_arr.num_elements(), arr_imags.data(), [](dcomplex a){ return imag(a); } ); 
   
   // Write arrays seperately as double arrays
   write( arr_reals, group, string("RE") + dataset_name ); 
   write( arr_imags, group, string("IM") + dataset_name ); 
}

template< size_t ndims > void read( boost::multi_array<double, ndims>& my_arr, H5::Group& group, const std::string& dataset_name )
{
   using namespace H5; 

   boost::array< hsize_t, ndims > dims; 

   
   //Open dataset and dataspace
   DataSet dataset = group.openDataSet( dataset_name );
   DataSpace dataspace = dataset.getSpace();


   int hndims = dataspace.getSimpleExtentDims( dims.data(), NULL );
//   std:: cout << "HNDIMS:" << hndims <<"NDIMS" << ndims << std::endl;
   assert( ndims == hndims );  // enforce equal dimensions
     
//
//   for(int i=0;i<hndims; ++i)
//     std:: cout << "EACH DIM:" << dims[i] << my_arr.shape()[i] << std::endl;
//
//   std:: cout << "NDIMS:" << ndims << std::endl;
//   // Resize my_arr if necessary
//   std:: cout << "NDIMS:" << my_arr.shape() << dims.begin() << std::endl; 

   if( !std::equal( my_arr.shape(), my_arr.shape() + ndims, dims.begin() ) )
   {
      std::cout << "Reading (double) required resize! Does not consider index-bases!" << std::endl; 
      my_arr.resize( dims ); 
   }
   
   dataset.read( my_arr.data(),PredType::NATIVE_DOUBLE ); 
}

template< size_t ndims > void read( boost::multi_array<dcomplex, ndims>& my_arr, H5::Group& group, const std::string& dataset_name )
{
   using namespace std; 
   
   // Create array for real and imaginary part with equal size
   boost::multi_array<double,ndims> arr_reals, arr_imags; 

   // Read real and imaginary part
   read( arr_reals, group, string("RE") + dataset_name ); 
   read( arr_imags, group, string("IM") + dataset_name ); 
    
   // Resize my_arr if necessary
   auto& dims = reinterpret_cast<const boost::array<size_t, ndims>&>(*arr_reals.shape()); 
   if( !std::equal( my_arr.shape(), my_arr.shape() + ndims, dims.begin() ) )
   {
      std::cout << "Reading (dcomplex) required resize! Does not consider index-bases!" << std::endl; 
      my_arr.resize( dims ); 
   }

   // Merge real and imaginary part into complex array
   transform( arr_reals.data(), arr_reals.data() + arr_reals.num_elements(), arr_imags.data(), my_arr.data(), []( double a, double b ){ return dcomplex( a, b ); } ); 
}

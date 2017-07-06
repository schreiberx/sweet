#ifndef _CENCAP_HPP_
#define _CENCAP_HPP_

#include <iomanip>
#include "PlaneDataVars.hpp"
#include "PlaneDataCtx.hpp"

/*
  "Encap" functions called from Fortran to manipulate PlaneData
*/


extern "C"
{
  // instantiates and returns the sweet data encapsulated object
  void c_sweet_data_create(
			   PlaneDataCtx *i_ctx, 
			   int i_level,
			   PlaneDataVars **o_Y, 
			   int *o_size
			   );

  // calls the destructor of the sweet data encapsulated object
  void c_sweet_data_destroy(
			    PlaneDataVars *i_Y
			    );

  // sets the value of the sweet data encapsulated object
  void c_sweet_data_setval(
			   PlaneDataVars *io_Y,  
			   double i_val
			   );

  // copies i_src into o_dst
  void c_sweet_data_copy(
			 PlaneDataVars *i_src,    
			 PlaneDataVars *o_dst
			 );

  // computes the norm of the sweet data encapsulated object
  void c_sweet_data_norm(
			 PlaneDataVars *io_Y,     
			 double* o_val
			 );

  // packs all the values contained in the sweet data object into a flat array
  void c_sweet_data_pack(
			 PlaneDataVars *io_Y,
			 double** o_flat_data_ptr
			 );
  // unpacks the flat array into the sweet data object 
  void c_sweet_data_unpack(
			   double** i_flat_data_ptr,
			   PlaneDataVars *io_Y
			   );
  
  // computes io_Y = i_a * i_X + io_Y
  void c_sweet_data_saxpy(
			  double i_a,          
			  PlaneDataVars *i_X,
			  PlaneDataVars *io_Y      
			  );

  // prints the data to the terminal
  void c_sweet_data_eprint(
			   PlaneDataVars *i_Y
			   );
}

#endif

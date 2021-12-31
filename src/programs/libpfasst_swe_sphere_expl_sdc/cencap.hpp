#ifndef _CENCAP_HPP_
#define _CENCAP_HPP_

#include <iomanip>
#include "../libpfasst_interface/SphereDataVars.hpp"
#include "SphereDataCtxSDC.hpp"

/*
  "Encap" functions called from Fortran to manipulate SphereData
*/


extern "C"
{
  // instantiates and returns the sweet data encapsulated object
  void c_sweet_data_create(
			   SphereDataCtxSDC *i_ctx, 
			   int i_level,
			   SphereDataVars **o_Y, 
			   int *o_size
			   );

  // calls the destructor of the sweet data encapsulated object
  void c_sweet_data_destroy(
			    SphereDataVars *i_Y
			    );

  // sets the value of the sweet data encapsulated object
  void c_sweet_data_setval(
			   SphereDataVars *io_Y,  
			   double i_val
			   );

  // copies i_src into o_dst
  void c_sweet_data_copy(
			 SphereDataVars *i_src,    
			 SphereDataVars *o_dst
			 );

  // computes the norm of the sweet data encapsulated object
  void c_sweet_data_norm(
			 SphereDataVars *io_Y,     
			 double* o_val
			 );

  // packs all the values contained in the sweet data object into a flat array
  void c_sweet_data_pack(
			 SphereDataVars *io_Y,
			 double** o_flat_data_ptr
			 );
  // unpacks the flat array into the sweet data object 
  void c_sweet_data_unpack(
			   double** i_flat_data_ptr,
			   SphereDataVars *io_Y
			   );
  
  // computes io_Y = i_a * i_X + io_Y
  void c_sweet_data_saxpy(
			  double i_a,          
			  SphereDataVars *i_X,
			  SphereDataVars *io_Y      
			  );

  // prints the data to the terminal
  void c_sweet_data_eprint(
			   SphereDataVars *i_Y
			   );
}

#endif

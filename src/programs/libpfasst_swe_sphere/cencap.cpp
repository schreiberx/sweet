#include <iomanip>
#include <cstddef>
#include "cencap.hpp"

extern "C"
{
  /*
    "Encap" functions called from Fortran to manipulate SphereData
   */
  
  // instantiates and returns the sweet data encapsulated object
  void c_sweet_data_create(
			   SphereDataCtx   *i_ctx, 
			   int             i_level, 
			   SphereDataVars **o_Y, 
			   int            *o_size
			   ) 
  {
    SphereDataConfig *Y_config = i_ctx->get_sphere_data_config(i_level);
    //Y_config->physical_array_data_number_of_elements = 
    // Y_config->physical_res[0]*Y_config->physical_res[1];

    // create the SphereDataVars object
    *o_Y  = new SphereDataVars(
			      Y_config,
			      i_level
			      );
    
    SphereData& phi  = (*o_Y)->get_phi();
    SphereData& vort = (*o_Y)->get_vort();
    SphereData& div  = (*o_Y)->get_div();

    // initialize the SphereData vectors
    phi.physical_set_zero();
    vort.physical_set_zero();
    div.physical_set_zero();

    // return the size of the number of elements 
    *o_size = (phi.sphereDataConfig->physical_array_data_number_of_elements 
	       + vort.sphereDataConfig->physical_array_data_number_of_elements 
	       + div.sphereDataConfig->physical_array_data_number_of_elements);
  }

  // calls the destructor of the sweet data encapsulated object
  void c_sweet_data_destroy(
			    SphereDataVars *i_Y
			    ) 
  {
    delete i_Y; // call the sweet object destructor
  }

  // sets the value of the sweet data encapsulated object
  void c_sweet_data_setval(
			   SphereDataVars *io_Y, 
			   double i_val
			   ) 
  {
    SphereData& phi  = io_Y->get_phi();
    SphereData& vort = io_Y->get_vort();
    SphereData& div  = io_Y->get_div();

    // set the SphereData vectors to i_val
    phi.physical_set_all_value(i_val);
    vort.physical_set_all_value(i_val);
    div.physical_set_all_value(i_val);
  }

  // copies i_src into o_dst
  void c_sweet_data_copy(SphereDataVars *i_src, 
			 SphereDataVars *o_dst) 
  {
    const SphereData& phi_src  = i_src->get_phi();
    const SphereData& vort_src = i_src->get_vort();
    const SphereData& div_src  = i_src->get_div();

    SphereData&       phi_dst  = o_dst->get_phi();
    SphereData&       vort_dst = o_dst->get_vort();
    SphereData&       div_dst  = o_dst->get_div();

    phi_dst = phi_src;
    vort_dst = vort_src;
    div_dst = div_src;
  }
  
  // computes the norm of the sweet data encapsulated object
  void c_sweet_data_norm(
			 SphereDataVars *i_Y, 
			 double *o_val
			 ) 
  {
    const SphereData& phi  = i_Y->get_phi();
    const SphereData& vort = i_Y->get_vort();
    const SphereData& div  = i_Y->get_div();

    *o_val = phi.physical_reduce_max_abs();
    const double vort_max = vort.physical_reduce_max_abs();
    const double div_max  = div.physical_reduce_max_abs();

    // L-infinity norm
    if (vort_max > *o_val) 
      *o_val = vort_max;
    if (div_max > *o_val)
      *o_val = div_max;
  }

  // packs all the values contained in the sweet data object into a flat array
  void c_sweet_data_pack(
			 SphereDataVars *io_Y, 
			 double **o_flat_data_ptr
			 ) 
  {
    SphereData& phi  = io_Y->get_phi();
    SphereData& vort = io_Y->get_vort();
    SphereData& div  = io_Y->get_div();

    // make sure that the physical data is up to date
    phi.request_data_physical();  
    vort.request_data_physical();  
    div.request_data_physical();  

    // allocate the flat data array
    const int n_elems = (phi.sphereDataConfig->physical_array_data_number_of_elements 
		      +  vort.sphereDataConfig->physical_array_data_number_of_elements
		      +  div.sphereDataConfig->physical_array_data_number_of_elements);
    io_Y->allocate_flat_data_array(n_elems);
    double*& flat_data_array = io_Y->get_flat_data_array();
    
    int j = 0;

    // phi
    for (int i = 0; i < phi.sphereDataConfig->physical_array_data_number_of_elements; ++i) 
      flat_data_array[j++] = phi.physical_space_data[i];

    // vort
    for (int i = 0; i < vort.sphereDataConfig->physical_array_data_number_of_elements; ++i)
      flat_data_array[j++] = vort.physical_space_data[i];

    // div
    for (int i = 0; i < div.sphereDataConfig->physical_array_data_number_of_elements; ++i)
      flat_data_array[j++] = div.physical_space_data[i]; 

    // return the pointer to the array
    *o_flat_data_ptr = flat_data_array;
  } 

  // unpacks the flat array into the sweet data object 
  void c_sweet_data_unpack(
			   double **i_flat_data_ptr, 
			   SphereDataVars *o_Y
			   ) 
  {
    int j = 0;

    // copy the values into physical_space_data array

    // phi
    SphereData& phi = o_Y->get_phi();
    for (int i = 0; i < phi.sphereDataConfig->physical_array_data_number_of_elements; ++i) 
      phi.physical_space_data[i] = i_flat_data_ptr[0][j++]; 

    // vort
    SphereData& vort = o_Y->get_vort();
    for (int i = 0; i < vort.sphereDataConfig->physical_array_data_number_of_elements; ++i) {
      vort.physical_space_data[i] = i_flat_data_ptr[0][j++]; 
    }

    // v
    SphereData& div = o_Y->get_div();
    for (int i = 0; i < div.sphereDataConfig->physical_array_data_number_of_elements; ++i) {
      div.physical_space_data[i] = i_flat_data_ptr[0][j++]; 
    }

    // tell sweet that the physical data is up to date
    phi.physical_space_data_valid  = true;
    vort.physical_space_data_valid = true;
    div.physical_space_data_valid  = true;

    // the spectral data needs to be recomputed
    phi.spectral_space_data_valid  = false;
    vort.spectral_space_data_valid = false;
    div.spectral_space_data_valid  = false;

    // make sure that the spectral data is up to date
    phi.request_data_spectral();  
    vort.request_data_spectral();  
    div.request_data_spectral();  
    
  }

  // computes io_Y = i_a * i_X + io_Y
  void c_sweet_data_saxpy(
			  double i_a, 
			  SphereDataVars *i_X, 
			  SphereDataVars *io_Y
			  ) 
  {
    const SphereData& phi_x = i_X->get_phi();
    const SphereData& vort_x = i_X->get_vort();
    const SphereData& div_x = i_X->get_div();

    SphereData&       phi_y = io_Y->get_phi();
    SphereData&       vort_y = io_Y->get_vort();
    SphereData&       div_y = io_Y->get_div();
    
    phi_y  = i_a * phi_x  + phi_y;
    vort_y = i_a * vort_x + vort_y;
    div_y  = i_a * div_x  + div_y;

  }

  // prints the data to the terminal
  void c_sweet_data_eprint(
			   SphereDataVars *i_Y
			   ) 
  {
    const SphereData& phi  = i_Y->get_phi();
    const SphereData& vort = i_Y->get_vort();
    const SphereData& div  = i_Y->get_div();
    
    phi.physical_print();
    vort.physical_print();
    div.physical_print();

  }
}


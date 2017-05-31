#include <iomanip>
#include <cstddef>
#include "cencap.hpp"

extern "C"
{
  /*
    "Encap" functions called from Fortran to manipulate PlaneData
   */
  
  // instantiates and returns the sweet data encapsulated object
  void c_sweet_data_create(
			   PlaneDataCtx   *i_ctx, 
			   int             i_level, 
			   PlaneDataVars **o_Y, 
			   int            *o_size
			   ) 
  {
    PlaneDataConfig *Y_config = i_ctx->get_plane_data_config(i_level);

    // create the PlaneDataVars object
    *o_Y  = new PlaneDataVars(
			      Y_config,
			      i_level
			      );
    
    PlaneData& phi = (*o_Y)->get_phi();
    PlaneData& u   = (*o_Y)->get_u();
    PlaneData& v   = (*o_Y)->get_v();

    // initialize the PlaneData vectors
    phi = 0.0;
    u   = 0.0;
    v   = 0.0;

    // return the size of the number of elements 
    *o_size = phi.planeDataConfig->physical_array_data_number_of_elements 
            +   u.planeDataConfig->physical_array_data_number_of_elements 
            +   v.planeDataConfig->physical_array_data_number_of_elements;
  }

  // calls the destructor of the sweet data encapsulated object
  void c_sweet_data_destroy(
			    PlaneDataVars *i_Y
			    ) 
  {
    delete i_Y; // call the sweet object destructor
  }

  // sets the value of the sweet data encapsulated object
  void c_sweet_data_setval(
			   PlaneDataVars *io_Y, 
			   double i_val
			   ) 
  {
    PlaneData& phi = io_Y->get_phi();
    PlaneData& u   = io_Y->get_u();
    PlaneData& v   = io_Y->get_v();

    // set the PlaneData vectors to i_val
    phi = i_val;
    u   = i_val;
    v   = i_val;
  }

  // copies i_src into o_dst
  void c_sweet_data_copy(PlaneDataVars *i_src, 
			 PlaneDataVars *o_dst) 
  {
    const PlaneData& phi_src = i_src->get_phi();
    const PlaneData& u_src   = i_src->get_u();
    const PlaneData& v_src   = i_src->get_v();

    PlaneData&       phi_dst = o_dst->get_phi();
    PlaneData&       u_dst   = o_dst->get_u();
    PlaneData&       v_dst   = o_dst->get_v();

    phi_dst = phi_src;
    u_dst   = u_src;
    v_dst   = v_src;
  }
  
  // computes the norm of the sweet data encapsulated object
  void c_sweet_data_norm(
			 PlaneDataVars *i_Y, 
			 double *o_val
			 ) 
  {
    const PlaneData& phi = i_Y->get_phi();
    const PlaneData& u   = i_Y->get_u();
    const PlaneData& v   = i_Y->get_v();

    *o_val = phi.reduce_maxAbs();
    const double u_max = u.reduce_maxAbs();
    const double v_max = v.reduce_maxAbs();

    // L-infinity norm
    if (u_max > *o_val) 
      *o_val = u_max;
    if (v_max > *o_val)
      *o_val = v_max;
  }

  // packs all the values contained in the sweet data object into a flat array
  void c_sweet_data_pack(
			 PlaneDataVars *io_Y, 
			 double **o_flat_data_ptr
			 ) 
  {
    PlaneData& phi = io_Y->get_phi();
    PlaneData& u   = io_Y->get_u();
    PlaneData& v   = io_Y->get_v();

    // make sure that the data is up to date
    phi.request_data_physical();  
    u.request_data_physical();  
    v.request_data_physical();  

    // allocate the flat data array
    const int n_elems = phi.planeDataConfig->physical_array_data_number_of_elements 
                      +   u.planeDataConfig->physical_array_data_number_of_elements 
                      +   v.planeDataConfig->physical_array_data_number_of_elements;
    io_Y->allocate_flat_data_array(n_elems);
    double*& flat_data_array = io_Y->get_flat_data_array();
    
    int j = 0;

    // phi
    for (int i = 0; i < phi.planeDataConfig->physical_array_data_number_of_elements; ++i) 
      flat_data_array[j++] = phi.physical_space_data[i];

    // u
    for (int i = 0; i < u.planeDataConfig->physical_array_data_number_of_elements; ++i) {
      flat_data_array[j++] = u.physical_space_data[i];
    }

    // v
    for (int i = 0; i < v.planeDataConfig->physical_array_data_number_of_elements; ++i) {
      flat_data_array[j++] = v.physical_space_data[i]; 
    }

    // return the pointer to the array
    *o_flat_data_ptr = flat_data_array;
  } 

  // unpacks the flat array into the sweet data object 
  void c_sweet_data_unpack(
			   double **i_flat_data_ptr, 
			   PlaneDataVars *o_Y
			   ) 
  {
    int j = 0;

    // copy the values into physical_space_data array

    // phi
    PlaneData& phi = o_Y->get_phi();
    for (int i = 0; i < phi.planeDataConfig->physical_array_data_number_of_elements; ++i) 
      phi.physical_space_data[i] = i_flat_data_ptr[0][j++]; 

    // u
    PlaneData& u = o_Y->get_u();
    for (int i = 0; i < u.planeDataConfig->physical_array_data_number_of_elements; ++i) {
      u.physical_space_data[i] = i_flat_data_ptr[0][j++]; 
    }

    // v
    PlaneData& v = o_Y->get_v();
    for (int i = 0; i < v.planeDataConfig->physical_array_data_number_of_elements; ++i) {
      v.physical_space_data[i] = i_flat_data_ptr[0][j++]; 
    }


    // tell sweet that the physical data is up to date
    phi.physical_space_data_valid = true;
    u.physical_space_data_valid   = true;
    v.physical_space_data_valid   = true;

    // the spectral data needs to be recomputed
    phi.spectral_space_data_valid = false;
    u.spectral_space_data_valid   = false;
    v.spectral_space_data_valid   = false;

  }

  // computes io_Y = i_a * i_X + io_Y
  void c_sweet_data_saxpy(
			  double i_a, 
			  PlaneDataVars *i_X, 
			  PlaneDataVars *io_Y
			  ) 
  {
    const PlaneData& phi_x = i_X->get_phi();
    const PlaneData& u_x   = i_X->get_u();
    const PlaneData& v_x   = i_X->get_v();

    PlaneData&       phi_y = io_Y->get_phi();
    PlaneData&       u_y   = io_Y->get_u();
    PlaneData&       v_y   = io_Y->get_v();
    
    phi_y = i_a * phi_x + phi_y;
    u_y   = i_a * u_x   + u_y;
    v_y   = i_a * v_x   + v_y;

  }

  // prints the data to the terminal
  void c_sweet_data_eprint(
			   PlaneDataVars *i_Y
			   ) 
  {
    const PlaneData& phi = i_Y->get_phi();
    const PlaneData& u   = i_Y->get_u();
    const PlaneData& v   = i_Y->get_v();
    
    phi.print_physicalArrayData();
    u.print_physicalArrayData();
    v.print_physicalArrayData();

  }
}


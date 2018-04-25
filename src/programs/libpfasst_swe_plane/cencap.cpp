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
    
    PlaneData& h    = (*o_Y)->get_h();
    PlaneData& vort = (*o_Y)->get_vort();
    PlaneData& div  = (*o_Y)->get_div();

    // initialize the PlaneData vectors
    h.physical_set_zero();
    vort.physical_set_zero();
    div.physical_set_zero();

    // return the size of the number of elements 
    *o_size = (h.planeDataConfig->physical_array_data_number_of_elements 
  	     + vort.planeDataConfig->physical_array_data_number_of_elements 
	     + div.planeDataConfig->physical_array_data_number_of_elements);
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
    PlaneData& h    = io_Y->get_h();
    PlaneData& vort = io_Y->get_vort();
    PlaneData& div  = io_Y->get_div();

    // set the PlaneData vectors to i_val
    h.physical_set_all(i_val);
    vort.physical_set_all(i_val);
    div.physical_set_all(i_val);
  }

  // copies i_src into o_dst
  void c_sweet_data_copy(PlaneDataVars *i_src, 
			 PlaneDataVars *o_dst) 
  {
    const PlaneData& h_src    = i_src->get_h();
    const PlaneData& vort_src = i_src->get_vort();
    const PlaneData& div_src  = i_src->get_div();

    PlaneData&       h_dst    = o_dst->get_h();
    PlaneData&       vort_dst = o_dst->get_vort();
    PlaneData&       div_dst  = o_dst->get_div();

    h_dst    = h_src;
    vort_dst = vort_src;
    div_dst  = div_src;
  }
  
  // computes the norm of the sweet data encapsulated object
  void c_sweet_data_norm(
			 PlaneDataVars *i_Y, 
			 double *o_val
			 ) 
  {
    const PlaneData& h    = i_Y->get_h();
    const PlaneData& vort = i_Y->get_vort();
    const PlaneData& div  = i_Y->get_div();

    *o_val = PlaneData(h).reduce_maxAbs();
    const double vort_max = PlaneData(vort).reduce_maxAbs();
    const double div_max  = PlaneData(div).reduce_maxAbs();

    // L-infinity norm
    // if (vort_max > *o_val) 
    //   *o_val = vort_max;
    // if (div_max > *o_val)
    //   *o_val = div_max;    
  }

  // packs all the values contained in the sweet data object into a flat array
  void c_sweet_data_pack(
			 PlaneDataVars *io_Y, 
			 double **o_flat_data_ptr
			 ) 
  {
    PlaneData& h    = io_Y->get_h();
    PlaneData& vort = io_Y->get_vort();
    PlaneData& div  = io_Y->get_div();

    // make sure that the physical data is up to date
    h.request_data_physical();  
    vort.request_data_physical();  
    div.request_data_physical();  

    // allocate the flat data array
    const int n_elems = (h.planeDataConfig->physical_array_data_number_of_elements 
		      +  vort.planeDataConfig->physical_array_data_number_of_elements
		      +  div.planeDataConfig->physical_array_data_number_of_elements);
    io_Y->allocate_flat_data_array(n_elems);
    double*& flat_data_array = io_Y->get_flat_data_array();
    
    int j = 0;

    // h
    for (int i = 0; i < h.planeDataConfig->physical_array_data_number_of_elements; ++i) 
      flat_data_array[j++] = h.physical_space_data[i];

    // vort
    for (int i = 0; i < vort.planeDataConfig->physical_array_data_number_of_elements; ++i)
      flat_data_array[j++] = vort.physical_space_data[i];

    // div
    for (int i = 0; i < div.planeDataConfig->physical_array_data_number_of_elements; ++i)
      flat_data_array[j++] = div.physical_space_data[i]; 

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

    // h
    PlaneData& h = o_Y->get_h();
    for (int i = 0; i < h.planeDataConfig->physical_array_data_number_of_elements; ++i) 
      h.physical_space_data[i] = i_flat_data_ptr[0][j++]; 

    // vort
    PlaneData& vort = o_Y->get_vort();
    for (int i = 0; i < vort.planeDataConfig->physical_array_data_number_of_elements; ++i) {
      vort.physical_space_data[i] = i_flat_data_ptr[0][j++]; 
    }

    // v
    PlaneData& div = o_Y->get_div();
    for (int i = 0; i < div.planeDataConfig->physical_array_data_number_of_elements; ++i) {
      div.physical_space_data[i] = i_flat_data_ptr[0][j++]; 
    }

    // tell sweet that the physical data is up to date
    h.physical_space_data_valid    = true;
    vort.physical_space_data_valid = true;
    div.physical_space_data_valid  = true;

    // the spectral data needs to be recomputed
    h.spectral_space_data_valid    = false;
    vort.spectral_space_data_valid = false;
    div.spectral_space_data_valid  = false;

    // make sure that the spectral data is up to date
    //h.request_data_spectral();  
    //vort.request_data_spectral();  
    //div.request_data_spectral();  
    
  }

  // computes io_Y = i_a * i_X + io_Y
  void c_sweet_data_saxpy(
			  double i_a, 
			  PlaneDataVars *i_X, 
			  PlaneDataVars *io_Y
			  ) 
  {
    const PlaneData& h_x    = i_X->get_h();
    const PlaneData& vort_x = i_X->get_vort();
    const PlaneData& div_x  = i_X->get_div();

    PlaneData&       h_y    = io_Y->get_h();
    PlaneData&       vort_y = io_Y->get_vort();
    PlaneData&       div_y  = io_Y->get_div();
    
    h_y  = i_a * h_x  + h_y;
    vort_y = i_a * vort_x + vort_y;
    div_y  = i_a * div_x  + div_y;

  }

  // prints the data to the terminal
  void c_sweet_data_eprint(
			   PlaneDataVars *i_Y
			   ) 
  {
    const PlaneData& h    = i_Y->get_h();
    const PlaneData& vort = i_Y->get_vort();
    const PlaneData& div  = i_Y->get_div();
    
    h.print_physicalArrayData();
    vort.print_physicalArrayData();
    div.print_physicalArrayData();

  }
}


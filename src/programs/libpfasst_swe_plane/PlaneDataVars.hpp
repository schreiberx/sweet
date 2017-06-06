#ifndef _PLANE_DATA_VARS_HPP_
#define _PLANE_DATA_VARS_HPP_

#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataConfig.hpp>

// Class containing the prognotic PlaneData variables h, u, v

class PlaneDataVars {

public:

  // Constructor
  PlaneDataVars(
                PlaneDataConfig *planeDataConfig,
		int i_level
		)

    : prog_h(planeDataConfig),
      prog_u(planeDataConfig),
      prog_v(planeDataConfig),
      flat_data_array(nullptr),
      flat_data_array_size(0),
      level(i_level)
  {}
  
  // Destructor
  ~PlaneDataVars()
  {
    if (flat_data_array != nullptr) 
      {
	// release the memory
	MemBlockAlloc::free(
			    flat_data_array, 
			    flat_data_array_size*sizeof(double)
			    );
      }
  }

  void allocate_flat_data_array(int i_n_elems) 
  {
    if (flat_data_array == nullptr) 
      {
	// allocate the memory
	flat_data_array_size = i_n_elems;
	flat_data_array      = MemBlockAlloc::alloc<double>(
							    i_n_elems*sizeof(double)
							    );
      }						   
  }
  
  // getters for the PlaneData variables
  const PlaneData& get_h() const {return prog_h;}
  PlaneData&       get_h()       {return prog_h;}
  const PlaneData& get_u() const {return prog_u;}
  PlaneData&       get_u()       {return prog_u;}
  const PlaneData& get_v() const {return prog_v;}
  PlaneData&       get_v()       {return prog_v;}

  // getters for the flat data array
  double*&         get_flat_data_array()            {return flat_data_array;}
  const int&       get_flat_data_array_size() const {return flat_data_array_size;}
  int&             get_flat_data_array_size()       {return flat_data_array_size;}

  // getters for the level 
  const int&       get_level() const {return level;};

protected:
  
  // height
  PlaneData prog_h;
  // velocities
  PlaneData prog_u;
  PlaneData prog_v;

  // flat data array vector (currently used to pack and unpack)
  double *flat_data_array;
  int     flat_data_array_size;
  
  // pfasst level
  const int level;

  // default constructor, copy constructor, and operator= are disabled
  PlaneDataVars();
  PlaneDataVars(const PlaneDataVars&);
  PlaneDataVars& operator=(const PlaneDataVars&);

};

#endif

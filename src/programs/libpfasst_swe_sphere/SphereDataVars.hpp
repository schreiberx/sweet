#ifndef _SPHERE_DATA_VARS_HPP_
#define _SPHERE_DATA_VARS_HPP_

#include <sweet/sphere/SphereData.hpp>
#include <sweet/sphere/SphereDataConfig.hpp>

// Class containing the prognotic SphereData variables h, u, v

class SphereDataVars {

public:

  // Constructor
  SphereDataVars(
                SphereDataConfig *sphereDataConfig,
		int i_level
		)

    : prog_h(sphereDataConfig),
      prog_u(sphereDataConfig),
      prog_v(sphereDataConfig),
      flat_data_array(nullptr),
      flat_data_array_size(0),
      level(i_level)
  {}
  
  // Destructor
  ~SphereDataVars()
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
  
  // getters for the SphereData variables
  const SphereData& get_h() const {return prog_h;}
  SphereData&       get_h()       {return prog_h;}
  const SphereData& get_u() const {return prog_u;}
  SphereData&       get_u()       {return prog_u;}
  const SphereData& get_v() const {return prog_v;}
  SphereData&       get_v()       {return prog_v;}

  // getters for the flat data array
  double*&         get_flat_data_array()            {return flat_data_array;}
  const int&       get_flat_data_array_size() const {return flat_data_array_size;}
  int&             get_flat_data_array_size()       {return flat_data_array_size;}

  // getters for the level 
  const int&       get_level() const {return level;};

protected:
  
  // height
  SphereData prog_h;
  // velocities
  SphereData prog_u;
  SphereData prog_v;

  // flat data array vector (currently used to pack and unpack)
  double *flat_data_array;
  int     flat_data_array_size;
  
  // pfasst level
  const int level;

  // default constructor, copy constructor, and operator= are disabled
  SphereDataVars();
  SphereDataVars(const SphereDataVars&);
  SphereDataVars& operator=(const SphereDataVars&);

};

#endif

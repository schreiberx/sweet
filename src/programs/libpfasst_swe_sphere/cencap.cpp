#include <iomanip>
#include <cstddef>
#include "cencap.hpp"

extern "C"
{
/*
    "Encap" functions called from Fortran to manipulate SphereDataSpectral
 */

// instantiates and returns the sweet data encapsulated object
void c_sweet_data_create(
		SphereDataCtx   *i_ctx,
		int             i_level,
		SphereDataVars **o_Y,
		int            *o_size
)
{
	SphereData_Config *Y_config = i_ctx->get_sphere_data_config(i_level);

	// create the SphereDataVars object
	*o_Y  = new SphereDataVars(
			Y_config,
			i_level
	);

	SphereData_Spectral& phi  = (*o_Y)->get_phi();
	SphereData_Spectral& vrt = (*o_Y)->get_vrt();
	SphereData_Spectral& div  = (*o_Y)->get_div();

	// initialize the SphereDataSpectral vectors
	phi.spectral_set_zero();
	vrt.spectral_set_zero();
	div.spectral_set_zero();

	// return the size of the number of elements
	*o_size = 2*(phi.sphereDataConfig->spectral_array_data_number_of_elements
			+ vrt.sphereDataConfig->spectral_array_data_number_of_elements
			+ div.sphereDataConfig->spectral_array_data_number_of_elements);
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
	SphereData_Spectral& phi  = io_Y->get_phi();
	SphereData_Spectral& vrt = io_Y->get_vrt();
	SphereData_Spectral& div  = io_Y->get_div();


	phi.spectral_set_zero();
	vrt.spectral_set_zero();
	div.spectral_set_zero();

	if (i_val == 0)
	{
		// set the SphereDataSpectral vectors to zero in spectral space
	}
	else
	{
		// set the SphereDataSpectral vectors to i_val in physical space
		phi.spectral_add_physical_constant(i_val);
		vrt.spectral_add_physical_constant(i_val);
		div.spectral_add_physical_constant(i_val);
	}
}


// copies i_src into o_dst
void c_sweet_data_copy(SphereDataVars *i_src,
		SphereDataVars *o_dst)
{
	const SphereData_Spectral& phi_src  = i_src->get_phi();
	const SphereData_Spectral& vrt_src  = i_src->get_vrt();
	const SphereData_Spectral& div_src  = i_src->get_div();

	SphereData_Spectral&       phi_dst  = o_dst->get_phi();
	SphereData_Spectral&       vrt_dst  = o_dst->get_vrt();
	SphereData_Spectral&       div_dst  = o_dst->get_div();

	phi_dst  = phi_src;
	vrt_dst = vrt_src;
	div_dst  = div_src;
}

// computes the norm of the sweet data encapsulated object
void c_sweet_data_norm(
		SphereDataVars *i_Y,
		double *o_val
)
{
	const SphereData_Spectral& phi  = i_Y->get_phi();
//	const SphereData_Spectral& vrt = i_Y->get_vrt();
//	const SphereData_Spectral& div  = i_Y->get_div();

	*o_val = phi.toPhys().physical_reduce_max_abs();
//	const double vrt_max = vrt.toPhys().physical_reduce_max_abs();
//	const double div_max  = div.toPhys().physical_reduce_max_abs();

	// L-infinity norm
	// if (vrt_max > *o_val)
	//   *o_val = vrt_max;
	// if (div_max > *o_val)
	//   *o_val = div_max;
}

// packs all the values contained in the sweet data object into a flat array
// void c_sweet_data_pack(
// 			 SphereDataVars *io_Y,
// 			 double **o_flat_data_ptr
// 			 )
// {
//   SphereDataSpectral& phi  = io_Y->get_phi();
//   SphereDataSpectral& vrt = io_Y->get_vrt();
//   SphereDataSpectral& div  = io_Y->get_div();

//   // make sure that the physical data is up to date
//   phi.request_data_physical();
//   vrt.request_data_physical();
//   div.request_data_physical();

//   // allocate the flat data array
//   const int n_elems = (phi.sphereDataConfig->physical_array_data_number_of_elements
// 		      +  vrt.sphereDataConfig->physical_array_data_number_of_elements
// 		      +  div.sphereDataConfig->physical_array_data_number_of_elements);
//   io_Y->allocate_flat_data_array(n_elems);
//   double*& flat_data_array = io_Y->get_flat_data_array();

//   int j = 0;

//   // phi
//   for (int i = 0; i < phi.sphereDataConfig->physical_array_data_number_of_elements; ++i)
//     flat_data_array[j++] = phi.physical_space_data[i];

//   // vrt
//   for (int i = 0; i < vrt.sphereDataConfig->physical_array_data_number_of_elements; ++i)
//     flat_data_array[j++] = vrt.physical_space_data[i];

//   // div
//   for (int i = 0; i < div.sphereDataConfig->physical_array_data_number_of_elements; ++i)
//     flat_data_array[j++] = div.physical_space_data[i];

//   // return the pointer to the array
//   *o_flat_data_ptr = flat_data_array;
// }

// packs all the values contained in the sweet data object into a flat array
void c_sweet_data_pack(
		SphereDataVars *io_Y,
		double **o_flat_data_ptr
)
{
	SphereData_Spectral& phi  = io_Y->get_phi();
	SphereData_Spectral& vrt  = io_Y->get_vrt();
	SphereData_Spectral& div  = io_Y->get_div();


	// allocate the flat data array
	const int n_elems = 2*(phi.sphereDataConfig->spectral_array_data_number_of_elements
			+  vrt.sphereDataConfig->spectral_array_data_number_of_elements
			+   div.sphereDataConfig->spectral_array_data_number_of_elements);
	io_Y->allocate_flat_data_array(n_elems);
	double*& flat_data_array = io_Y->get_flat_data_array();

	int j = 0;

	// real and imaginary parts

	// phi
	for (int i = 0; i < phi.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
	{
		flat_data_array[j++] = phi.spectral_space_data[i].imag();
		flat_data_array[j++] = phi.spectral_space_data[i].real();
	}

	// vrt
	for (int i = 0; i < vrt.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
	{
		flat_data_array[j++] = vrt.spectral_space_data[i].imag();
		flat_data_array[j++] = vrt.spectral_space_data[i].real();
	}

	// div
	for (int i = 0; i < div.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
	{
		flat_data_array[j++] = div.spectral_space_data[i].imag();
		flat_data_array[j++] = div.spectral_space_data[i].real();
	}

	// return the pointer to the array
	*o_flat_data_ptr = flat_data_array;
}

// // unpacks the flat array into the sweet data object
// void c_sweet_data_unpack(
// 			   double **i_flat_data_ptr,
// 			   SphereDataVars *o_Y
// 			   )
// {
//   int j = 0;

//   // copy the values into physical_space_data array

//   // phi
//   SphereDataSpectral& phi = o_Y->get_phi();
//   for (int i = 0; i < phi.sphereDataConfig->physical_array_data_number_of_elements; ++i)
//     phi.physical_space_data[i] = i_flat_data_ptr[0][j++];

//   // vrt
//   SphereDataSpectral& vrt = o_Y->get_vrt();
//   for (int i = 0; i < vrt.sphereDataConfig->physical_array_data_number_of_elements; ++i) {
//     vrt.physical_space_data[i] = i_flat_data_ptr[0][j++];
//   }

//   // v
//   SphereDataSpectral& div = o_Y->get_div();
//   for (int i = 0; i < div.sphereDataConfig->physical_array_data_number_of_elements; ++i) {
//     div.physical_space_data[i] = i_flat_data_ptr[0][j++];
//   }

//   // tell sweet that the physical data is up to date
//   phi.physical_space_data_valid  = true;
//   vrt.physical_space_data_valid = true;
//   div.physical_space_data_valid  = true;

//   // the spectral data needs to be recomputed
//   phi.spectral_space_data_valid  = false;
//   vrt.spectral_space_data_valid = false;
//   div.spectral_space_data_valid  = false;

//   // make sure that the spectral data is up to date
//   //phi.request_data_spectral();
//   //vrt.request_data_spectral();
//   //div.request_data_spectral();

// }

// unpacks the flat array into the sweet data object
void c_sweet_data_unpack(
		double **i_flat_data_ptr,
		SphereDataVars *o_Y
)
{
	int j = 0;

	// copy the values into physical_space_data array

	// phi
	SphereData_Spectral& phi = o_Y->get_phi();
	for (int i = 0; i < phi.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
	{
		phi.spectral_space_data[i] = std::complex<double>(
				i_flat_data_ptr[0][j++],
				i_flat_data_ptr[0][j++]
		);
	}

	// vrt
	SphereData_Spectral& vrt = o_Y->get_vrt();
	for (int i = 0; i < vrt.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
	{
		vrt.spectral_space_data[i] = std::complex<double>(
				i_flat_data_ptr[0][j++],
				i_flat_data_ptr[0][j++]
		);
	}

	// v
	SphereData_Spectral& div = o_Y->get_div();
	for (int i = 0; i < div.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
	{
		div.spectral_space_data[i] = std::complex<double>(
				i_flat_data_ptr[0][j++],
				i_flat_data_ptr[0][j++]
		);

	}
}


// computes io_Y = i_a * i_X + io_Y
void c_sweet_data_saxpy(
		double i_a,
		SphereDataVars *i_X,
		SphereDataVars *io_Y
)
{
	const SphereData_Spectral& phi_x  = i_X->get_phi();
	const SphereData_Spectral& vrt_x = i_X->get_vrt();
	const SphereData_Spectral& div_x  = i_X->get_div();

	SphereData_Spectral&       phi_y  = io_Y->get_phi();
	SphereData_Spectral&       vrt_y = io_Y->get_vrt();
	SphereData_Spectral&       div_y  = io_Y->get_div();

	phi_y  = i_a * phi_x  + phi_y;
	vrt_y = i_a * vrt_x + vrt_y;
	div_y  = i_a * div_x  + div_y;

}

// prints the data to the terminal
void c_sweet_data_eprint(
		SphereDataVars *i_Y
)
{
	const SphereData_Spectral& phi  = i_Y->get_phi();
	const SphereData_Spectral& vrt = i_Y->get_vrt();
	const SphereData_Spectral& div  = i_Y->get_div();

	phi.toPhys().physical_print();
	vrt.toPhys().physical_print();
	div.toPhys().physical_print();

}
}


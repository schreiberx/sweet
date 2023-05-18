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
        SphereDataCtxSDC *i_ctx,
        int i_level,
        SphereDataVars **o_Y,
        int *o_size)
    {
        sweet::SphereData_Config *Y_config = i_ctx->get_sphere_data_config();

        // create the SphereDataVars object
        *o_Y = new SphereDataVars(
            Y_config,
            i_level);

        sweet::SphereData_Spectral &phi_pert = (*o_Y)->get_phi_pert();
        sweet::SphereData_Spectral &vrt = (*o_Y)->get_vrt();
        sweet::SphereData_Spectral &div = (*o_Y)->get_div();

        // initialize the SphereDataSpectral vectors
        phi_pert.spectral_set_zero();
        vrt.spectral_set_zero();
        div.spectral_set_zero();

        // return the size of the number of elements
        *o_size = 2 * (phi_pert.sphereDataConfig->spectral_array_data_number_of_elements + vrt.sphereDataConfig->spectral_array_data_number_of_elements + div.sphereDataConfig->spectral_array_data_number_of_elements);
    }

    // calls the destructor of the sweet data encapsulated object
    void c_sweet_data_destroy(
        SphereDataVars *i_Y)
    {
        delete i_Y; // call the sweet object destructor
    }

    // sets the value of the sweet data encapsulated object
    void c_sweet_data_setval(
        SphereDataVars *io_Y,
        double i_val)
    {
        sweet::SphereData_Spectral &phi_pert = io_Y->get_phi_pert();
        sweet::SphereData_Spectral &vrt = io_Y->get_vrt();
        sweet::SphereData_Spectral &div = io_Y->get_div();

        phi_pert.spectral_set_zero();
        vrt.spectral_set_zero();
        div.spectral_set_zero();

        if (i_val == 0)
        {
            // set the SphereDataSpectral vectors to zero in spectral space
        }
        else
        {
            // set the SphereDataSpectral vectors to i_val in physical space
            phi_pert.spectral_add_physical_constant(i_val);
            vrt.spectral_add_physical_constant(i_val);
            div.spectral_add_physical_constant(i_val);
        }
    }

    // copies i_src into o_dst
    void c_sweet_data_copy(SphereDataVars *i_src,
                           SphereDataVars *o_dst)
    {
        const sweet::SphereData_Spectral &phi_pert_src = i_src->get_phi_pert();
        const sweet::SphereData_Spectral &vrt_src = i_src->get_vrt();
        const sweet::SphereData_Spectral &div_src = i_src->get_div();

        sweet::SphereData_Spectral &phi_pert_dst = o_dst->get_phi_pert();
        sweet::SphereData_Spectral &vrt_dst = o_dst->get_vrt();
        sweet::SphereData_Spectral &div_dst = o_dst->get_div();

        phi_pert_dst = phi_pert_src;
        vrt_dst = vrt_src;
        div_dst = div_src;
    }

    // computes the norm of the sweet data encapsulated object
    void c_sweet_data_norm(
        SphereDataVars *i_Y,
        double *o_val)
    {
        const sweet::SphereData_Spectral &phi_pert = i_Y->get_phi_pert();

        *o_val = phi_pert.toPhys().physical_reduce_max_abs();
    }

    // packs all the values contained in the sweet data object into a flat array
    void c_sweet_data_pack(
        SphereDataVars *io_Y,
        double **o_flat_data_ptr)
    {
        sweet::SphereData_Spectral &phi_pert = io_Y->get_phi_pert();
        sweet::SphereData_Spectral &vrt = io_Y->get_vrt();
        sweet::SphereData_Spectral &div = io_Y->get_div();

        // allocate the flat data array
        const int n_elems = 2 * (phi_pert.sphereDataConfig->spectral_array_data_number_of_elements + vrt.sphereDataConfig->spectral_array_data_number_of_elements + div.sphereDataConfig->spectral_array_data_number_of_elements);
        io_Y->allocate_flat_data_array(n_elems);
        double *&flat_data_array = io_Y->get_flat_data_array();

        int j = 0;

        // real and imaginary parts

        // phi_pert
        for (int i = 0; i < phi_pert.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
        {
            flat_data_array[j++] = phi_pert.spectral_space_data[i].imag();
            flat_data_array[j++] = phi_pert.spectral_space_data[i].real();
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

    // unpacks the flat array into the sweet data object
    void c_sweet_data_unpack(
        double **i_flat_data_ptr,
        SphereDataVars *o_Y)
    {
        int j = 0;

        // copy the values into physical_space_data array

        // phi_pert
        sweet::SphereData_Spectral &phi_pert = o_Y->get_phi_pert();
        for (int i = 0; i < phi_pert.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
        {
            phi_pert.spectral_space_data[i] = std::complex<double>(
                i_flat_data_ptr[0][j],
                i_flat_data_ptr[0][j + 1]);
            j += 2;
        }

        // vrt
        sweet::SphereData_Spectral &vrt = o_Y->get_vrt();
        for (int i = 0; i < vrt.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
        {
            vrt.spectral_space_data[i] = std::complex<double>(
                i_flat_data_ptr[0][j],
                i_flat_data_ptr[0][j + 1]);
            j += 2;
        }

        // div
        sweet::SphereData_Spectral &div = o_Y->get_div();
        for (int i = 0; i < div.sphereDataConfig->spectral_array_data_number_of_elements; ++i)
        {
            div.spectral_space_data[i] = std::complex<double>(
                i_flat_data_ptr[0][j],
                i_flat_data_ptr[0][j + 1]);
            j += 2;
        }
    }

    // computes io_Y = i_a * i_X + io_Y
    void c_sweet_data_saxpy(
        double i_a,
        SphereDataVars *i_X,
        SphereDataVars *io_Y)
    {
        const sweet::SphereData_Spectral &phi_pert_x = i_X->get_phi_pert();
        const sweet::SphereData_Spectral &vrt_x = i_X->get_vrt();
        const sweet::SphereData_Spectral &div_x = i_X->get_div();

        sweet::SphereData_Spectral &phi_pert_y = io_Y->get_phi_pert();
        sweet::SphereData_Spectral &vrt_y = io_Y->get_vrt();
        sweet::SphereData_Spectral &div_y = io_Y->get_div();

        phi_pert_y = i_a * phi_pert_x + phi_pert_y;
        vrt_y = i_a * vrt_x + vrt_y;
        div_y = i_a * div_x + div_y;
    }

    // prints the data to the terminal
    void c_sweet_data_eprint(
        SphereDataVars *i_Y)
    {
        const sweet::SphereData_Spectral &phi_pert = i_Y->get_phi_pert();
        const sweet::SphereData_Spectral &vrt = i_Y->get_vrt();
        const sweet::SphereData_Spectral &div = i_Y->get_div();

        phi_pert.toPhys().physical_print();
        vrt.toPhys().physical_print();
        div.toPhys().physical_print();
    }
}

/*
 * Author: Martin SCHREIBER <schreiberx@gmail.com>
 */

#ifndef SRC_PROGRAMS_PDE_SWE_SPHERE_FILE_OUTPUT_HPP_
#define SRC_PROGRAMS_PDE_SWE_SPHERE_FILE_OUTPUT_HPP_


// This is just for the editor to show code as used within precompiler #if ... directives
#include <sweet/core/defaultPrecompilerValues.hpp>

// Error handling
#include <sweet/core/ErrorBase.hpp>

// Include everything we need for simulations on the plane
#include <sweet/core/sphere/Sphere.hpp>

// Different shacks we need in this file
#include <sweet/core/shacksShared/ShackSphereDataOps.hpp>
#include <sweet/core/shacksShared/ShackIOData.hpp>
#include <sweet/core/shacksShared/ShackTimestepControl.hpp>

#include <sweet/core/Convert_SphereDataSpectral_To_PlaneDataPhysical.hpp>
#include <sweet/core/Convert_SphereDataPhysical_To_PlaneDataPhysical.hpp>

#include "ShackPDESWESphere.hpp"


class PDESWESphereFileOutput
{
public:
	sweet::ErrorBase error;

	sweet::ShackIOData *shackIOData;
	sweet::ShackTimestepControl *shackTimestepControl;
	ShackPDESWESphere *shackPDESWESphere;

	void setup(
			sweet::ShackIOData *i_shackIOData,
			sweet::ShackTimestepControl *i_shackTimestepControl,
			ShackPDESWESphere *i_shackPDESWESphere
	)
	{
		shackIOData = i_shackIOData;
		shackTimestepControl = i_shackTimestepControl;
		shackPDESWESphere = i_shackPDESWESphere;
	}

	void clear()
	{
		shackIOData = nullptr;
		shackTimestepControl = nullptr;
		shackPDESWESphere = nullptr;
	}

	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv_spec_evol(
			const sweet::SphereData_Spectral &i_sphereData,
			const char* i_name		///< name of output variable
	)
	{
		char buffer[1024];
		std::string phase = "_arg";
		std::string ampl = "_amp"; 

		const char* filename_template_ampl = "output_spec_ampl_%s.txt"; //.c_str();
		const char* filename_template_arg = "output_spec_arg_%s.txt"; //.c_str();
		int reduce_mode_factor = 4;

		sprintf(buffer, filename_template_arg, i_name);
		i_sphereData.spectrum_phase_file_write_line(buffer, 
			i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale,
			20, 10e-20, reduce_mode_factor);

		sprintf(buffer, filename_template_ampl, i_name);
		i_sphereData.spectrum_abs_file_write_line(buffer, 
			i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale,
			20, 10e-20, reduce_mode_factor);

		return buffer;
	}


	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_csv(
			const sweet::SphereData_Spectral &i_sphereData,
			const char* i_name,		///< name of output variable
			bool i_phi_shifted = false
	)
	{
		char buffer[1024];

		// create copy
		sweet::SphereData_Physical sphereData = i_sphereData.toPhys();

		const char* filename_template = shackIOData->output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale);

		if (i_phi_shifted)
			sphereData.physical_file_write_lon_pi_shifted(buffer, "vorticity, lon pi shifted");
		else
			sphereData.physical_file_write(buffer);

		return buffer;
	}



	/**
	 * Write file to data and return string of file name
	 */
	std::string write_file_bin(
			const sweet::SphereData_Spectral &i_sphereData,
			const char* i_name
	)
	{
		char buffer[1024];

		sweet::SphereData_Spectral sphereData(i_sphereData);
		const char* filename_template = shackIOData->output_file_name.c_str();
		sprintf(buffer, filename_template, i_name, shackTimestepControl->current_simulation_time*shackIOData->output_time_scale);
		sphereData.file_write_binary_spectral(buffer);

		return buffer;
	}


	std::string output_reference_filenames;

	void write_file_output(
			sweet::SphereOperators &i_ops,
			sweet::SphereData_Spectral &i_prog_phi_pert,
			sweet::SphereData_Spectral &i_prog_div,
			sweet::SphereData_Spectral &i_prog_vrt
	)
	{
#if SWEET_MPI
		if (mpi_rank > 0)
			return;
#endif

		if (shackIOData->output_file_name.length() == 0)
			return;

		std::cout << "Writing output files at simulation time: " << shackTimestepControl->current_simulation_time << " secs" << std::endl;

		if (shackIOData->output_file_mode == "csv")
		{
			std::string output_filename;

			sweet::SphereData_Spectral h = i_prog_phi_pert*(1.0/shackPDESWESphere->gravitation);
			h += shackPDESWESphere->h0;

			output_filename = write_file_csv(h, "prog_h");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << " (min: " << h.toPhys().physical_reduce_min() << ", max: " << h.toPhys().physical_reduce_max() << ")" << std::endl;

			output_filename = write_file_csv(i_prog_phi_pert, "prog_phi_pert");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << i_prog_phi_pert.toPhys().physical_reduce_min() << ", max: " << i_prog_phi_pert.toPhys().physical_reduce_max() << ")" << std::endl;

			sweet::SphereData_Physical phi_phys = h.toPhys() * shackPDESWESphere->gravitation;
			sweet::SphereData_Spectral phi(i_ops.sphereDataConfig);
			phi.loadSphereDataPhysical(phi_phys);
			output_filename = write_file_csv(phi, "prog_phi");
			output_reference_filenames = output_filename;
			std::cout << " + " << output_filename << " (min: " << phi_phys.physical_reduce_min() << ", max: " << phi_phys.physical_reduce_max() << ")" << std::endl;

			sweet::SphereData_Physical u(i_ops.sphereDataConfig);
			sweet::SphereData_Physical v(i_ops.sphereDataConfig);

			i_ops.vrtdiv_to_uv(i_prog_vrt, i_prog_div, u, v);

			output_filename = write_file_csv(u, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(v, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(i_prog_vrt, "prog_vrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			output_filename = write_file_csv(i_prog_div, "prog_div");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;

			sweet::SphereData_Spectral potvrt = (i_prog_phi_pert/shackPDESWESphere->gravitation)*i_prog_vrt;

			output_filename = write_file_csv(potvrt, "prog_potvrt");
			output_reference_filenames += ";"+output_filename;
			std::cout << " + " << output_filename << std::endl;
		}
		else if (shackIOData->output_file_mode == "bin")
		{
			std::string output_filename;

			{
				output_filename = write_file_bin(i_prog_phi_pert, "prog_phi_pert");
				output_reference_filenames = output_filename;
				sweet::SphereData_Physical prog_phys = i_prog_phi_pert.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}

			{
				output_filename = write_file_bin(i_prog_vrt, "prog_vrt");
				output_reference_filenames += ";"+output_filename;
				sweet::SphereData_Physical prog_phys = i_prog_vrt.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}

			{
				output_filename = write_file_bin(i_prog_div, "prog_div");
				output_reference_filenames += ";"+output_filename;
				sweet::SphereData_Physical prog_phys = i_prog_div.toPhys();

				std::cout << " + " << output_filename << " (min: " << prog_phys.physical_reduce_min() << ", max: " << prog_phys.physical_reduce_max() << ")" << std::endl;
			}
		}
		else if (shackIOData->output_file_mode == "csv_spec_evol"){

			std::string output_filename;

			{ 
				/*
				* Spectral kinetic energy and potential enstrophy calculation and output
				*
				* Details in Jakob-Chien, Ruediger, James J. Hack, and David L. Williamson. 
				* "Spectral transform solutions to the shallow water test set." Journal of Computational Physics 119, no. 1 (1995): 164-187.
				*/
				// Kinetic energy is given in spectral space as
				// KE per mode = a^2/((n(n+1)))*(vrt*conj(vrt))+a^2/((n(n+1)))*(div*conj(div))
				// r = a/(sqrt(n(n+1))) (root_laplace)
				// KE per mode = (r*vrt*conj(r*vrt))+(r*div*conj(r*div))
				sweet::SphereData_Spectral rlap_vrt = i_ops.inv_root_laplace(i_prog_vrt);
				sweet::SphereData_Spectral rlap_div = i_ops.inv_root_laplace(i_prog_div);
				sweet::SphereData_Spectral kin_en = rlap_vrt + rlap_div ;

				output_filename = write_file_csv_spec_evol(kin_en, "kin_en"); 
				std::cout << " + " << output_filename << " (Total Kin Energy : " << 0.25*kin_en.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// For Barotropic vort eq: See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere (2009) for eps=0
				// Kinetic energy is given in spectral space as
				// Vortical energy per mode is (0.5 n*(n+1) / a^2) *psi*conj(psi) in spectral space
				//SphereData_Spectral psi = op.inv_laplace(prog_vrt); // 
				// multiply psi by sqrt( n * (n+1))/a (apply root laplacian)
				//SphereData_Spectral psi_root = op.root_laplace(psi);
				//output_filename = write_file_csv_spec_evol(psi_root*std::sqrt(0.5), "spec_energy"); 
				//std::cout << " + " << output_filename << " (Kinetic energy : " << (0.5)*psi_root.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

				// See Schubert Shallow Water Quasi-Geostrophic Theory on the Sphere (2009) for eps=0
				// enstrophy per mode is 0.5 vrt*conj(vrt) in spectral space
				// Total enstrophy is the sum of these (counting twice modes with m>0 and once when m=0)
				output_filename = write_file_csv_spec_evol(i_prog_vrt, "enstrophy");
				std::cout << " + " << output_filename << " (Total Enstrophy : " << i_prog_vrt.spectral_reduce_sum_sqr_quad() << ")" << std::endl;

			}
		}
		else
		{
			SWEETError("Unknown output file mode '"+shackIOData->output_file_mode+"'");
		}
	}
};

#endif

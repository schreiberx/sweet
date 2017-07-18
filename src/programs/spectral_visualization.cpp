
#if !SWEET_USE_PLANE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif


#if SWEET_GUI != 1
	#error "Activate GUI to run this program"
#endif


#include <sweet/plane/PlaneData.hpp>
#include <sweet/plane/PlaneDataComplex.hpp>
#include <sweet/plane/PlaneOperators.hpp>
#include <sweet/VisSweet.hpp>
#include <sweet/SimulationVariables.hpp>
#include <sweet/Stopwatch.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>
#include <stdio.h>

// Plane data config
PlaneDataConfig planeDataConfigInstance;
PlaneDataConfig *planeDataConfig = &planeDataConfigInstance;



SimulationVariables simVars;

class TestSpectral
{
public:
	PlaneData tmp;
	PlaneOperators op;
	char vis_description[1024];

	/**
	 * visualization mode
	 * 0: visualize some spectras
	 * 1: visualize particular distributions
	 */
	int vis_mode = 0;


	Stopwatch stopwatch;

public:
	TestSpectral()	:
		tmp(planeDataConfig),

		op(planeDataConfig, simVars.sim.domain_size, simVars.disc.use_spectral_basis_diffs)
	{
		stopwatch.reset();

		planeDataConfig->printInformation();
		simVars.outputConfig();

#if 0
		//for (std::size_t j = 0; j < planeDataConfig->spectral_complex_data_size[1]; j++)
		std::size_t j = 0;
		{
			for (std::size_t i = 0; i < planeDataConfig->spectral_complex_data_size[0]; i++)
			//std::size_t i = 0;
			{
				{
					std::cout << std::endl;
					std::cout << "complex: " << j << ", " << i << std::endl;
					PlaneDataComplex tmp(planeDataConfig);
					tmp.spectral_set_zero();
					tmp.p_spectral_set(j, i, 1);
					tmp.print_physicalArrayData();
				}

				{
					std::cout << std::endl;
					std::cout << "real: " << j << ", " << i << std::endl;
					PlaneData tmp(planeDataConfig);
					tmp.spectral_set_zero();
					tmp.p_spectral_set(j, i, 1);
					tmp.print_physicalArrayData();
				}

			}
		}
#endif

		exit(1);
	}


	void reset()
	{
		// ugly, but quick
		sprintf(vis_description, "[initialization, press 'v' and 'V' to change the spectrum]");
	}


	void run_timestep()
	{
	}


public:
	void timestep_output(
			std::ostream &o_ostream = std::cout
	)
	{
	}



public:
	bool should_quit()
	{
		return false;
	}


	/**
	 * postprocessing of frame: do time stepping
	 */
	void vis_post_frame_processing(
			int i_num_iterations
	)
	{
	}



	void vis_get_vis_data_array(
			const PlaneData **o_dataArray,
			double *o_aspect_ratio,
			int *o_render_primitive,
			void **o_bogus_data
	)
	{
		*o_dataArray = &tmp;
		*o_aspect_ratio = simVars.sim.domain_size[1] / simVars.sim.domain_size[0];

		if (simVars.timecontrol.run_simulation_timesteps)
		{
			tmp.spectral_set_all(0, 0);

			int spec_array[][4] =
			{
					// j, i, re, im
					{0, 0, 1, 0},
					{0, 0, 0, 1},
					{0, 0, 1, 1},

					{1, 0, 1, 0},
					{1, 0, 0, 1},
					{1, 0, 1, 1},

					{2, 0, 1, 0},
					{2, 0, 0, 1},
					{2, 0, 1, 1},

					{0, 1, 1, 0},
					{0, 1, 0, 1},
					{0, 1, 1, 1},

					{1, 1, 1, 0},
					{1, 1, 0, 1},
					{1, 1, 1, 1},

					{(int)simVars.disc.res_physical[1]-1, 0, 1, 0},
					{(int)simVars.disc.res_physical[1]-1, 0, 0, 1},
					{(int)simVars.disc.res_physical[1]-1, 0, 1, 1},

					{(int)simVars.disc.res_physical[1]-2, 0, 1, 0},
					{(int)simVars.disc.res_physical[1]-2, 0, 0, 1},
					{(int)simVars.disc.res_physical[1]-2, 0, 1, 1},

					{0, (int)simVars.disc.res_physical[0]/2, 1, 0},
					{0, (int)simVars.disc.res_physical[0]/2, 0, 1},
					{0, (int)simVars.disc.res_physical[0]/2, 1, 1},

					{0, (int)simVars.disc.res_physical[0]/2-1, 1, 0},
					{0, (int)simVars.disc.res_physical[0]/2-1, 0, 1},
					{0, (int)simVars.disc.res_physical[0]/2-1, 1, 1},

					{0, (int)simVars.disc.res_physical[0]/2-2, 1, 0},
					{0, (int)simVars.disc.res_physical[0]/2-2, 0, 1},
					{0, (int)simVars.disc.res_physical[0]/2-2, 1, 1},
			};

			int id = simVars.misc.vis_id % (sizeof(spec_array)/sizeof(spec_array[0]));

			sprintf(vis_description, "spec_coord (j, i) = (%i, %i), value = %i + i*%i", spec_array[id][0], spec_array[id][1], spec_array[id][2], spec_array[id][3]);
			tmp.p_spectral_set(spec_array[id][0], spec_array[id][1], {(double)spec_array[id][2], (double)spec_array[id][3]});

			return;
		}

		PlaneData dataArray(planeDataConfig);

		double max = 1.0;
		double shift = stopwatch.getTimeSinceStart()*0.1;

		double c = cos(2.0*M_PIl*shift);
		double s = sin(2.0*M_PIl*shift);
		c *= max;
		s *= max;

		sprintf(vis_description, "test spectrum");
		tmp.spectral_set_all(0, 0);
		int asdf = 2;
		tmp.p_spectral_set(0, asdf, {c, s});
		tmp.p_spectral_set(asdf, 0, {c, s});
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		static char output_string[2048];

		sprintf(output_string, "%i: %s, res x/y: %i/%i, inv: %f", simVars.misc.vis_id, vis_description, simVars.disc.res_physical[0], simVars.disc.res_physical[1], 1.0/(double)(simVars.disc.res_physical[0]*simVars.disc.res_physical[1]));
		return output_string;
	}



	void vis_pause()
	{
		simVars.timecontrol.run_simulation_timesteps = !simVars.timecontrol.run_simulation_timesteps;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			simVars.misc.vis_id++;
			break;

		case 'V':
			simVars.misc.vis_id--;
			break;
		}
	}


};



int main(int i_argc, char *i_argv[])
{
	if (!simVars.setupFromMainParameters(i_argc, i_argv))
	{
		return -1;
	}

	planeDataConfigInstance.setupAuto(simVars.disc.res_physical, simVars.disc.res_spectral);

	TestSpectral *testSpectral = new TestSpectral;

	VisSweet<TestSpectral> visSweet(testSpectral);

	delete testSpectral;

	return 0;
}

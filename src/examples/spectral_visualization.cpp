
#if !SWEET_USE_SPECTRAL_SPACE
	#error "Spectral space not activated"
#endif


#if SWEET_GUI != 1
	#error "Activate GUI to run this program"
#endif


#include <sweet/DataArray.hpp>
#include <sweet/Operators2D.hpp>
#include <sweet/VisSweet.hpp>
#include <sweet/SimulationParameters.hpp>
#include <sweet/Stopwatch.hpp>

#include <ostream>
#include <sstream>
#include <unistd.h>
#include <iomanip>

SimulationParameters parameters;

class TestSpectral
{
public:
	DataArray<2> tmp;
	Operators2D op;
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
		tmp(parameters.res),

		op(parameters.res, parameters.sim_domain_size, parameters.use_spectral_diffs)
	{
		stopwatch.reset();
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
			const DataArray<2> **o_dataArray,
			double *o_aspect_ratio
	)
	{
		*o_dataArray = &tmp;
		*o_aspect_ratio = parameters.sim_domain_size[1] / parameters.sim_domain_size[0];

		if (parameters.run_simulation)
		{
			tmp.spec_setAll(0, 0);

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

					{(int)parameters.res[1]-1, 0, 1, 0},
					{(int)parameters.res[1]-1, 0, 0, 1},
					{(int)parameters.res[1]-1, 0, 1, 1},

					{(int)parameters.res[1]-2, 0, 1, 0},
					{(int)parameters.res[1]-2, 0, 0, 1},
					{(int)parameters.res[1]-2, 0, 1, 1},

					{0, (int)parameters.res[0]/2, 1, 0},
					{0, (int)parameters.res[0]/2, 0, 1},
					{0, (int)parameters.res[0]/2, 1, 1},

					{0, (int)parameters.res[0]/2-1, 1, 0},
					{0, (int)parameters.res[0]/2-1, 0, 1},
					{0, (int)parameters.res[0]/2-1, 1, 1},

					{0, (int)parameters.res[0]/2-2, 1, 0},
					{0, (int)parameters.res[0]/2-2, 0, 1},
					{0, (int)parameters.res[0]/2-2, 1, 1},
			};

			int id = parameters.vis_id % (sizeof(spec_array)/sizeof(spec_array[0]));

			sprintf(vis_description, "spec_coord (j, i) = (%i, %i), value = %i + i*%i", spec_array[id][0], spec_array[id][1], spec_array[id][2], spec_array[id][3]);
			tmp.spec_set(spec_array[id][0], spec_array[id][1], spec_array[id][2], spec_array[id][3]);

			return;
		}

		DataArray<2> dataArray(parameters.res);

		double max = 1.0;
		double shift = stopwatch.getTimeSinceStart();

		double c = cos(2.0*M_PIl*shift);
		double s = sin(2.0*M_PIl*shift);
		c *= max;
		s *= max;

		sprintf(vis_description, "test spectrum");
		tmp.spec_setAll(0, 0);
		int asdf = 10;
		tmp.spec_spectrum_set(0, asdf, c, s);
		tmp.spec_spectrum_set(asdf, 0, c, s);
	}



	/**
	 * return status string for window title
	 */
	const char* vis_get_status_string()
	{
		static char output_string[2048];

		sprintf(output_string, "%i: %s, res x/y: %lu/%lu, inv: %f", parameters.vis_id, vis_description, parameters.res[0], parameters.res[1], 1.0/(double)(parameters.res[0]*parameters.res[1]));
		return output_string;
	}



	void vis_pause()
	{
		parameters.run_simulation = !parameters.run_simulation;
	}



	void vis_keypress(int i_key)
	{
		switch(i_key)
		{
		case 'v':
			parameters.vis_id++;
			break;

		case 'V':
			parameters.vis_id--;
			break;
		}
	}


};



int main(int i_argc, char *i_argv[])
{
	parameters.setup(i_argc, i_argv);

	TestSpectral *testSpectral = new TestSpectral;

	VisSweet<TestSpectral> visSweet(testSpectral);

	delete testSpectral;

	return 1;
}

#include "PlaneDataVars.hpp"
#include "PlaneDataCtx.hpp"
#include "ceval.hpp"

extern "C"
{
  void cecho_error(PlaneData* sd, 
		   int step)
  {
    // not implemented
  }

  void cecho_residual(PlaneData* sd, 
                      int step, 
	              int iter)
  {
    // not implemented
  }

  void cecho_output_solution(
			     PlaneDataCtx *i_ctx,
			     PlaneDataVars *i_Y,
			     int i_current_iter,
			     int i_nnodes,
			     int i_niters
			     )
  {
    const PlaneData& h_Y = i_Y->get_h();
    const PlaneData& u_Y = i_Y->get_u();
    const PlaneData& v_Y = i_Y->get_v();

    // get the SimulationVariables object from context
    SimulationVariables* simVars(i_ctx->get_simulation_variables());

    // write the data to file
    std::string filename = "prog_h_current_iter_"+std::to_string(i_current_iter)
                                +"_nnodes_"      +std::to_string(i_nnodes)
                                +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, h_Y, filename.c_str());

    filename = "prog_u_current_iter_"+std::to_string(i_current_iter)
                    +"_nnodes_"      +std::to_string(i_nnodes)
                    +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, u_Y, filename.c_str());

    filename = "prog_v_current_iter_"+std::to_string(i_current_iter)
                    +"_nnodes_"      +std::to_string(i_nnodes)
                    +"_niters_"      +std::to_string(i_niters);
    write_file(*i_ctx, v_Y, filename.c_str());

  }

}
